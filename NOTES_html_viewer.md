# The interactive HTML tree viewer

Clann can emit a **single self-contained interactive HTML file** that renders its
trees in a web browser — pan/zoom, reroot, collapse, search, and export, with no
server, internet connection, or plug-in required. Four commands produce it via an
`htmlview=` option: `showtrees` (the input gene trees), `hs` (the best
supertree(s)), `nj` (the neighbour-joining tree), and `reconstruct` (all
reconciled gene trees, with duplication/loss events marked).

This note explains why it exists, how it is built and embedded, the data format
it consumes, the public C API, and what is left to do.

---

## 1. Why it exists

Clann's native tree output is ASCII art on the terminal plus Newick/NHX files.
For anything but the smallest tree that is hard to read, and reconciliations —
which layer duplication and loss events onto a gene tree — are worse, because the
events are invisible in plain Newick.

The reconciliation ecosystem does have a standard interchange format,
**recPhyloXML**, and a good static renderer, **thirdkind** (recPhyloXML → SVG,
including the "gene tree inside the species-tree tube" embedding). But those are
static images produced by a separate toolchain. There was no *interactive*
reconciliation viewer that a user could open straight from a Clann run — reroot,
collapse a clade, search for a taxon — without installing anything.

The design goal was therefore a **zero-dependency, offline, single-file**
artefact: something Clann writes at the end of a command that the user can
double-click, email, or archive, and that still works years later with no server
or CDN to go stale.

---

## 2. Architecture: template → generator → embedded strings

The viewer is authored as one ordinary HTML page and **baked into the Clann
binary** as two C string constants that bracket the per-run data.

- **Source of truth:** `tools/clannview.template.html` — the whole viewer (HTML +
  CSS + vanilla JS, no external assets). Edit this to change the viewer.
- **Generator:** `python3 tools/gen_viewer_template.py` splits the template at the
  markers `/*CLANN_DATA_BEGIN*/ … /*CLANN_DATA_END*/` into two C strings,
  `VIEWER_HTML_HEAD` and `VIEWER_HTML_TAIL`, and writes `viewer_template.h`.
  **Re-run the generator after every template edit, then `make`.**
- **Embedding:** `reconcile.c` does `#include "viewer_template.h"` and writes
  `VIEWER_HTML_HEAD` + `<JSON data>` + `VIEWER_HTML_TAIL` to the output file.

So a produced `.html` file is: the fixed viewer chrome, then a JSON blob
describing this run's trees, then the closing tags. All the interactivity lives
in the template's JavaScript; the C side only serialises trees to JSON.

---

## 3. Data format

Between HEAD and TAIL, Clann writes a single JSON object:

```json
{ "type": "tree" | "reconciliation",
  "meta": { "dataset": "…", "criterion": "recon", "lossmodel": "standard" },
  "trees": [ { "name": "tree_0", "score": 4, "dups": 1, "losses": 3,
               "tree": <node> }, … ] }
```

- `type` is `"reconciliation"` for `reconstruct` (events present) and `"tree"`
  for the plain-tree commands.
- `score`/`dups`/`losses` are emitted only for reconciliations.
- A file with a single element in `trees` hides the navigator; more than one
  shows it.

Each `<node>`:

```json
{ "name": "Human.1", "species": "Human", "length": 0.05, "support": 95,
  "event": "speciation" | "duplication" | "loss", "children": [ … ] }
```

- `event` and `species` are present only for reconciliations. `species` is the
  mapped species; a `loss` node is a lost lineage.
- Leaf labels are the **full per-copy taxon names** (e.g. `Human.alpha`), not the
  collapsed species — the viewer shows exactly what was in the source tree.
- The viewer computes a bottom-up `lost` flag (a node all of whose descendant
  leaves are losses) and **dots the branches** of such fully-lost clades.

---

## 4. Public C API

Declared in `reconcile.h`, defined in `reconcile.c`:

```c
FILE *html_view_open(const char *filename, const char *metajson, int recon);
void  html_view_add_tree (FILE *f, struct taxon *tree, const char *name,
                          float score, int recon, int first);
void  html_view_add_newick(FILE *f, const char *newick, const char *name,
                           int treenum, int first);
void  html_view_close(FILE *f, const char *filename);
```

Usage is always: **open once, add each tree, close once.** `first` must be `1`
(TRUE) for the first tree added and `0` afterwards (it controls the JSON comma
between array elements). A `NULL` file (open failed) makes every add/close a
safe no-op.

There are two ways to add a tree:

- **`html_view_add_tree` — from an in-memory annotated tree** (`recon=1`). Used by
  `reconstruct`, which has already built and annotated `tree_top` (duplication =
  `loss==2`, loss = `loss==-1`). It emits the per-tree score/dups/losses. This
  path works for both loss models, because `reconstruct` annotates `tree_top` in
  both (`annotate_standard` for `standard`, `add_losses` for `legacy`) before the
  illustration step.
- **`html_view_add_newick` — from a Newick string** (`recon=0`, plain tree). Two
  input conventions, distinguished by `treenum`:
  - `treenum >= 0`: `newick` is gene tree `treenum`'s numeric form
    (`fundamentals[treenum]`); it is expanded to full taxon names via
    `returntree_fullnames`. Used by `showtrees`.
  - `treenum == -1`: `newick` already carries full names (e.g. a supertree in
    `retained_supers[]`); it is parsed as-is. Used by `hs` and `nj`.

    In both cases the tree is built with `basic_tree_build(…, TRUE)`, which
    stores the literal leaf label in each node's `fullname` — this is why
    per-paralog names such as `Human.alpha` survive into the viewer rather than
    collapsing to the species id.

Static helpers in `reconcile.c` do the serialisation: `hv_json_node`,
`hv_json_tree`, `hv_json_str`, `hv_count_dups`.

### Which commands emit it, and default filenames

| Command | `type` | `htmlview=yes` default name |
|---------|--------|-----------------------------|
| `showtrees` | tree | `<input>.trees.html` |
| `hs` | tree | `<input>.supertree.html` |
| `nj` | tree | `<input>.nj.html` |
| `reconstruct` | reconciliation | `<input>.recon.html` |

(These default output names are gitignored.)

---

## 5. Viewer features

All implemented in the template's JavaScript, driven entirely by the embedded
JSON — nothing is computed by Clann at view time:

- **Navigator** (previous/next, dropdown, `filter trees…` box) when more than one
  tree is present. The filter matches tree **name**, tree **number**, or a
  **contained taxon name**.
- **Cladogram / phylogram** layout; adjustable row spacing and font size.
- Toggles for **support values, branch lengths, internal labels, aligned tip
  labels**.
- **Click-to-reroot** (Reroot mode) and click-a-node-dot to **collapse/expand**.
  Rerooting a *reconciliation* flags its mapping **stale**, because the
  duplication/loss events were computed for the original rooting and are not
  recomputed in the browser.
- A **Find** box that highlights taxa.
- **Export**: the current view as **SVG**, or a **Newick string reflecting the
  current (possibly rerooted) rooting**.
- **Light/dark** theme.
- Fully-lost clades have their branches **dotted** (see §3).

---

## 6. Testing a produced file

```
python3 -m http.server 8777 --bind 127.0.0.1      # in the dir with the .html
```

Then drive it with the browser tools. The JS exposes handles useful for
inspection: `TREES`, `root`, `loadTree(i)`, `applyTreeFilter(q)`. Screenshot
coordinates are scaled and unreliable for clicking — use element refs from the
accessibility tree, or evaluate JS directly.

---

## 7. Future work

- **Support values through the plain-tree path.** `html_view_add_newick` does not
  yet carry bootstrap/support labels into the JSON `support` field; wire them
  through so `hs`/`nj`/`showtrees` trees show support like reconciliations can.
- **Species-tree "tube" embedding.** Render the gene tree nested inside the
  species tree (the thirdkind style), which makes duplications and losses read at
  a glance. This needs the species tree alongside each reconciled gene tree in
  the data.
- **A `recphyloxml=` exporter** for ecosystem interoperability. recPhyloXML is the
  standard reconciliation interchange format; NHX cannot express the species
  embedding. Emitting recPhyloXML would let Clann reconciliations be rendered by
  `thirdkind` and consumed by other reconciliation tools.
