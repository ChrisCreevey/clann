# Gene-tree Decomposition into Ortholog Subtrees (`decomposegenetrees`)

This note is an implementation plan for a new top-level command,
`decomposegenetrees`, plus an `autodecompose=yes|no` auto-run flag. It was
designed in a planning conversation (not yet implemented) and is written so
that a fresh Claude Code session — with no memory of that conversation — can
pick it up and build it. Every referenced function/variable has been read
directly from the current codebase; file:line references are accurate as of
this writing but re-check them before editing, since line numbers drift.

---

## 1. Motivation

Clann can already score a candidate supertree against a set of source trees
using several criteria (`scoring.h`): `dfit`, `sfit`, `qfit`, `mrp`, `avcon`,
`rf`, `ml`, plus the duplication/loss reconciliation criterion `recon`
(`reconcile.c:977`, `get_recon_score`). A separate planning discussion
concluded that `criterion=qfit` (`scoring.c:1373`) is, mathematically, already
close to an ASTRAL-style quartet score: it sums per-source-tree quartet
agreement between the candidate tree (pruned to each source tree's taxon set)
and each source tree's induced quartets. Under the multispecies coalescent
(MSC), maximizing that sum is the same objective ASTRAL is built to optimize,
which is what makes `qfit` usable as a coalescent-consistent criterion without
reimplementing ASTRAL's search machinery.

The problem: the MSC assumes gene tree discordance comes only from incomplete
lineage sorting (ILS) between single-copy orthologs. Multi-copy gene family
trees (paralogs from duplication/loss/HGT) violate that assumption outright —
feeding them to `qfit` (or a future deep-coalescence/MDC criterion) conflates
duplication-driven discordance with ILS-driven discordance and biases the
result. Clann's *existing* whole-tree filter for this, `apply_singlecopy_filter`
(`treecompare2.c:2993`), is all-or-nothing: one in-paralog anywhere in a family
tree discards the entire tree's signal from the search.

`decomposegenetrees` is meant to do better: decompose each multi-copy gene
family tree into maximal single-copy-ish **ortholog subtrees**, salvaging
signal instead of discarding whole trees, using machinery Clann already has
for gene-tree/species-tree reconciliation.

### Terminology (for whoever implements this)

- **Orthologs**: genes related by a speciation event.
- **Paralogs**: genes related by a duplication event.
  - **In-paralogs** (relative to a duplication node D): copies that arose
    *after* D, redundant within one descendant lineage of D.
  - **Out-paralogs** (relative to D): the two (or more) clades *on either side*
    of D — these must never be treated as orthologous to each other.
- Decomposition means: cut at out-paralog splits, collapse in-paralog
  redundancy, and do neither where the evidence is too weak or the resulting
  pieces would be uninformative.

---

## 2. Existing building blocks to reuse (do not reimplement these)

| Purpose | Function | Location | Notes |
|---|---|---|---|
| LCA-map a gene tree onto a species tree, tagging each node with its mapped species-tree node id | `label_gene_tree` | `reconcile.c:518` | Sets `position->tag`. Uses `get_min_node` for the LCA. |
| Detect duplication nodes from the LCA mapping | `reconstruct_map` | `reconcile.c:588` | A node is a duplication iff a child's `tag` equals the node's own `tag` (`reconcile.c:606-611`). Marks it via `position->loss = 2` (`reconcile.c:617`). **This flag is the exact hook for finding cut points.** |
| Full duplication+loss reconciliation cost for one gene tree against one species tree (searches over gene-tree rootings) | `tree_map` | `reconcile.c:131` | Calls `label_gene_tree` + `reconstruct_map` + `add_losses`/`count_losses`. Read this to understand the calling convention, not to reuse its cost formula. |
| Full reconciliation search over gene-tree AND species-tree rootings (used by `criterion=recon`) | `get_recon_score` | `reconcile.c:977` | Shows the pattern for rerooting an unrooted gene tree and picking the best mapping (`reconcile.c:1097-1129`, the `best_mapping`/`best_total` loop). Decomposition should reuse this "find the best rooting of the gene tree against the guide tree" logic before walking for cut points, since duplication calls are rooting-dependent. |
| Detect species-specific (pure in-paralog) clades **without any reference species tree** | `identify_species_specific_clades` | `prune.c:356` | Recursive, checks both below (`list_taxa_in_clade`, `prune.c:496`) and above (`list_taxa_above`, `prune.c:456`) a node to catch monophyly obscured by arbitrary rooting. Tags a "clann ID" per species-specific clade. |
| Whole-command wrapper around the above: prune every multi-copy source tree down to one representative per species-specific clade | `prune_monophylies` | `prune.c:196` | Existing command `prunemonophylies`. Writes `prunedtrees.txt` + `prunedtrees_info.txt`. Drops (excludes) any resulting tree below 4 taxa (`prune.c:331`). Representative choice: random, or longest sequence via `selection=length` (parses embedded length from `fullname`, `extract_length`, `prune.c:531`). **Stage 1 of the new command should call into this same logic** (see §5 on refactoring needed). |
| Whole-tree multi-copy filter (current crude alternative) | `apply_singlecopy_filter` / `restore_singlecopy_filter` | `treecompare2.c:2993` / `treecompare2.c:3029` | Binary: excludes any tree with `presence_of_taxa[i][j] > 1` for any taxon `j`. Toggles `sourcetreetag[i]`. Good reference for how "exclude vs include in search" is represented, but decomposition should *replace* trees, not just toggle tags. |
| Preserve pristine trees when `fundamentals[]` gets modified in place | `original_fundamentals[]` | declared/used in `main.c:2225-2317`, swapped in `reconcile.c:1504-1518` and swapped back at `reconcile.c:2087` | **Important invariant**: `autoprunemono` mutates `fundamentals[i]` in place and stashes the pristine version in `original_fundamentals[i]`. `reconstruct` swaps the originals back in before running, then swaps back after. `decomposegenetrees` must do the same swap-in before Stage 1 runs, so it always decomposes from the pristine tree, never from an already-mutated copy. |
| Auto-run-at-load-time pattern to mirror for `autodecompose=` | `autoprunemono_apply` | `main.c:2215`, dispatched from `main.c:3247` via `do_autoprunemono` flag parsed at `main.c:2811` | Exact template for `autodecompose`. |
| "Species tree source" resolution pattern to mirror for the guide tree | `reconstruct`'s `speciestree=memory\|first\|<file>` option | parsing at `reconcile.c:1467-1477`, resolution at `reconcile.c:1585-1639` | `memory` reads `retained_supers[0]` (`reconcile.c:1620`); explicit file is read raw (`reconcile.c:1626-1630`); with nothing supplied and nothing in memory, `reconstruct` **hard-errors** (`reconcile.c:1592`). `decomposegenetrees` should reuse the `speciestree=` option name/values for consistency, but change the "nothing available" behavior: instead of erroring, fall back to calling `nj()` (see §4). |
| Auto-generate a guide tree with no criterion/search needed | `nj()` | `reconcile.c:42` | Builds an NJ tree from average-consensus distances over whatever's currently loaded, writes to `retained_supers[0]`, sets `trees_in_memory = 1` (`reconcile.c:97-110`). This is the fallback tier of guide-tree resolution. |
| Per-source-tree weight already threaded through every criterion | `tree_weights[]` | used in `scoring.c:978,1145,1246,1349,1448` | Multiplied into every criterion's per-tree contribution. This is the hook for down-weighting fragments so one duplication-heavy family doesn't out-vote families that never duplicated. |
| The source-tree pool itself | `fundamentals[]`, `Total_fund_trees`, `presence_of_taxa[][]`, `sourcetreetag[]`, `tree_names[]` | throughout `treecompare2.c`/`scoring.c` | Decomposed subtrees are just new entries in this pool — nothing about `hs`, `qfit`, `rf`, etc. needs to change once the pool is populated correctly. |
| Species identity from a multi-copy taxon name | delimiter parsing | `tree_io.c:817-830` | When `delimiter` mode is on, only the substring before `delimiter_char` is used to resolve a leaf to its `taxa_names[]` index (`tree_io.c:820-825,837`). **This means `Human.alpha` and `Human.beta` already resolve to the same numeric taxon id as `Human`** — `label_gene_tree`/`reconstruct_map`'s duplication test (`tag` equality) works correctly on multi-copy input with zero extra species-lookup code. Per-copy distinctness survives only in `fullname` (used for `select_longest`). |
| Existing example dataset with exactly the edge cases discussed | `examples/tutorial_multicopy.ph` | repo root `examples/` | See §7 for a line-by-line breakdown — this file already contains pure in-paralog cases and scattered/out-paralog cases suitable as test fixtures. |

---

## 3. Command surface

### `decomposegenetrees [options]`

A new manual command, dispatched the same way `prunemonophylies` is
(`main.c:1232-1240`: `if(strcmp(parsed_command[0], "prunemonophylies") == 0) prune_monophylies();` — add an equivalent branch calling the new function; add `"decomposegenetrees"` to the command-name array near `main.c:55`; add help text near `main.c:1563` and `main.c:2013`).

| Option | Values | Default | Notes |
|---|---|---|---|
| `speciestree` | `memory` \| `<file>` | unset → resolved per §4 | Reuse `reconstruct`'s option name/semantics for consistency. |
| `dupsupport` | float | TBD, conservative (see §6) | Minimum confidence to trust a duplication call before cutting on it. Units/source of "confidence" need deciding during implementation — see open question in §6. |
| `minfragtaxa` | int | 4 | Reuses the existing quartet floor already hardcoded in `scoring.c:1409` (`ntaxa_i < 4`) — expose it as a shared constant/option rather than a second hardcoded literal. |
| `minfragspecies` | int | 2 | A fragment mapping to fewer than this many distinct species is dropped/merged, not kept as its own tree. |
| `selection` | `length` \| `random` | `random` | Same semantics as `prune_monophylies`'s `selection=length` (`prune.c:219`) for in-paralog representative choice — reuse directly, don't reimplement. |
| `filename` | `<name>` | `decomposedtrees.txt` | Mirrors `prune_monophylies`'s `filename=` (`prune.c:221`). Output: `<filename>` (the decomposed subtrees) and `<filename>_info.txt` (log), same pattern as `prunedtrees.txt`/`prunedtrees_info.txt`. |

### `autodecompose=yes|no`

Global option alongside `autoprunemono` (`main.c:65,152`), applied at load
time the same way (`main.c:3247`), or optionally triggered automatically the
first time a coalescent criterion (`qfit`, and later a prospective `mdc`) is
invoked on data that still contains multi-copy trees. Exact trigger point
(load time vs. first search) is an implementation decision — load time is
simpler and matches `autoprunemono`'s existing precedent, so default to that
unless it proves awkward.

---

## 4. Guide-tree resolution order

Agreed order, in priority:

1. **`speciestree=<file>`** supplied explicitly → use it, but **hard-error** if
   its taxon set does not exactly match the full registry (`taxa_names[]` /
   `number_of_taxa`). Check both directions: every name in `taxa_names[]` must
   appear in the guide tree, and the guide tree must contain no taxon absent
   from `taxa_names[]`. Fail loudly and list what's missing/extra — do not
   silently prune to the intersection. (`reconstruct`'s file-reading code at
   `reconcile.c:1623-1631` is the base to extend with this validation; it
   currently only checks the file opens, not that taxa match.)
2. **`speciestree=memory`**, or nothing supplied but `trees_in_memory > 0` →
   use `retained_supers[0]` (`reconcile.c:1620` shows the exact read).
3. **Nothing available** → call `nj()` (`reconcile.c:42`) to generate one on
   the fly. This is the one deliberate deviation from `reconstruct`'s existing
   behavior, which hard-errors in this case (`reconcile.c:1592`) — for
   `decomposegenetrees` a missing guide tree should never be fatal, since the
   whole point is that this command can run standalone with zero setup.

---

## 5. Algorithm

### Stage 0 — restore pristine input

Before anything else, if `autoprunemono_active` and `original_fundamentals[]`
is populated, swap the pristine trees back into `fundamentals[]` (same pattern
as `reconcile.c:1504-1518`), run the whole pipeline, then swap back at the end
(`reconcile.c:2087` shows the swap-back). Decomposition must always start from
the pristine gene tree, never from something already collapsed by
`autoprunemono`/`prunemonophylies` in a previous step of the same session.

### Stage 1 — monophyly collapse (no guide tree needed)

For every multi-copy source tree, run the same logic as `prune_monophylies`
(`prune.c:196`). **Refactoring note**: `prune_monophylies` today is a single
monolithic command function that loops over `fundamentals[]` itself, builds
each tree, calls `identify_species_specific_clades`, and writes output files
directly. Extract the per-tree "collapse species-specific clades in this tree,
return the collapsed tree + a log of what was removed" logic
(`prune.c:238-343`) into a reusable function that both `prune_monophylies` and
`decomposegenetrees` call, rather than duplicating the loop body. Trees that
collapse fully to single-copy go straight into the candidate pool with weight
1. Trees still multi-copy after this stage proceed to Stage 2.

### Stage 2 — LCA-mapping decomposition (needs the guide tree)

For each tree still multi-copy after Stage 1:

1. Find the best rooting of the gene tree against the guide tree, reusing the
   rerooting-and-remap-cost pattern from `get_recon_score`
   (`reconcile.c:1097-1129`) — map via `label_gene_tree`, score via whatever
   cost is cheapest to compute (duplication count alone is sufficient here,
   losses aren't needed for the decomposition decision).
2. Run `reconstruct_map` (`reconcile.c:588`) to get the `loss == 2`
   duplication tags.
3. Walk the tagged tree top-down. At each duplication node D, apply, in order:
   - **Support gate**: if the call at D isn't well-supported (see open
     question in §6 on what "supported" means here), treat D as *not* a
     duplication — do not cut, continue the recursion through D as if it were
     an ordinary node.
   - **Pure in-paralog check**: if every leaf under D maps to a single
     species, collapse to one representative (reuse the `selection=`
     representative-choice logic from Stage 1/`prune_monophylies` — this case
     should be rare here since Stage 1 already removed most of it, but keep
     the check for defense-in-depth against cases Stage 1's topology-only test
     can miss, e.g. gene-tree-error-scattered in-paralogs that later map back
     together under the guide tree).
   - **Informativeness floor**: evaluate each child clade of D independently
     against `minfragtaxa`/`minfragspecies`. If one side fails the floor,
     prune that side down to a representative tip and merge it back into its
     sibling — do **not** spin off an uninformative fragment. If both sides
     fail, collapse the whole clade under D to one representative.
   - **Otherwise (both sides pass, well-supported, genuinely cross-species)**:
     cut. Each child of D becomes the root of an independent subtree; recurse
     into each separately, repeating this same decision procedure.
4. Base case: a subtree whose root and every retained descendant node is a
   speciation node (no further duplication tags) is a finished ortholog
   subtree — stop recursing.

### Stage 3 — weighting and pool integration

For every original family tree that produced `k` surviving fragments (`k >=
1`), set each fragment's `tree_weights[i] = 1/k` (simple equal split for a
first implementation — do not over-engineer a size-proportional scheme yet,
see §6). Add each surviving fragment (Stage 1 pass-throughs included, with
weight 1) as a new entry in `fundamentals[]`/`Total_fund_trees`, with its own
row in `presence_of_taxa[][]` computed the normal way trees are loaded
elsewhere (look at how `tree_build`/taxon registration populates
`presence_of_taxa` for a freshly-parsed tree — do not hand-roll this).

### Output

- `<filename>` — the decomposed subtrees, one Newick per line, same physical
  format as `prunedtrees.txt`.
- `<filename>_info.txt` — a log of every decision made: which nodes were
  collapsed (in-paralog), which were pruned-and-merged (thin side), which were
  cut (with the resulting fragment count and assigned weight), and which
  duplication calls were suppressed by the support gate. This is what makes
  the command useful standalone, not just as internal search plumbing — a
  user should be able to read this file and understand exactly what happened
  to their gene family trees.

---

## 6. Lifecycle: applying the output, backing up originals, and restoring

This section resolves a question that came up during planning: once Stage 3
produces the decomposed fragment set, what actually happens to the live
`fundamentals[]`, and how do we get back to the pre-decomposition state later?

### 6.1 Two invocation modes, two destructiveness defaults

Clann already has precedent for exactly this fork, in the
`prune_monophylies`/`autoprunemono` pair, and the two halves behave
differently on purpose:

- **`prune_monophylies`** (standalone command) is **non-destructive**: it
  reads `fundamentals[]`, computes the pruning, and only writes
  `prunedtrees.txt`/`prunedtrees_info.txt` to disk. It never reassigns
  `fundamentals[i]`. The user must explicitly `exe prunedtrees.txt` to adopt
  the result.
- **`autoprunemono=yes`** (load-time flag) **is destructive/immediate**: as
  part of the same `exe`, it mutates `fundamentals[]` in place and stashes the
  pristine copy in `original_fundamentals[]` (`main.c:2214` comment,
  `main.c:2225-2317`) so `reconstruct` can still get at the originals within
  the same session.

`decomposegenetrees`/`autodecompose` should mirror this split exactly:

- **`decomposegenetrees` (manual command) defaults to non-destructive.**
  Write `<filename>` and `<filename>_info.txt` only; do not touch live
  `fundamentals[]`; print a message telling the user to `exe <filename>` if
  they want to adopt it. This needs none of the backup machinery in §6.3,
  because nothing in memory changes. (Optionally add an `apply=yes` option
  later for one-step adoption, which would then go through the same path as
  `autodecompose` below — but don't build this in the first pass; it's not
  needed for the core feature to work.)
- **`autodecompose=yes` commits immediately**, in-memory, as part of the same
  `exe` — this is the path that needs the backup/restore design below.

### 6.2 Why the reload must go through the *same* registration path as loading a new file — not a hand-rolled patch

`Total_fund_trees` changes after decomposition (one family tree can become
`k` fragments). That invalidates every piece of per-tree derived state
computed at load time for the old count: `presence_of_taxa[][]`,
`fund_bipart_sets` (`scoring.c:994`, `rf_precompute_fund_biparts`), the `dfit`
distance matrices (`cal_fund_scores`), `sourcetree_scores[]`, `tree_weights[]`
sizing, taxon hash values, etc. Re-deriving all of that by hand in a second
code path, parallel to what `execute_command` already does for a freshly
loaded file (`main.c:2750` onward), is exactly the kind of duplicated logic
this plan has tried to avoid elsewhere.

Instead: **write the decomposed fragments to a file, then call the existing
full load/registration path on that file** — `execute_command(filename, TRUE)`
directly, or the public wrapper `clann_load_trees(filename, parse_opts)`
(`clann_api.c:303`, confirmed to be a thin wrapper: it just builds `"exe
<file>"` and calls `execute_command`). This is not a workaround, it's the
correct way to get a fully-consistent re-registration for free, using the
one code path Clann already trusts to do this correctly.

This also directly answers "do the preprocessing steps need to be re-run" —
yes, genuinely, twice, for an unavoidable structural reason, not just for
convenience: Stage 2's guide-tree resolution (§4) needs the *original* trees
already registered (to fall back to `nj()`, which itself needs
`presence_of_taxa`/average-consensus distances already computed) before
decomposition can even run — so the first registration pass (of the raw
input) has to happen before decomposition, and a second pass (of the
decomposed output) has to happen after, because the tree count changed.

### 6.3 The backup problem, and why a new mechanism is needed (not `original_fundamentals[]`)

`execute_command` unconditionally tears down and frees `fundamentals[]` *and*
`original_fundamentals[]` at the start of every load (`main.c:2883-2907`,
including `autoprunemono_active = 0` at `main.c:2907`). That means if
`autodecompose` reloads through the standard path per §6.2, the pre-decompose
pristine trees are gone from memory the instant the reload happens — unless
something snapshots them first.

`original_fundamentals[]` can't be reused for this snapshot: it's sized and
indexed 1:1 against the *old* tree count, and it's about to be freed by the
very reload it would need to survive.

**Required new state** (names indicative, pick better ones during
implementation):

- `pre_decompose_fundamentals[]` / `pre_decompose_tree_names[]` /
  `pre_decompose_total_trees` — a snapshot of the complete pristine gene-tree
  set as it stood immediately before decomposition (i.e. after Stage 0's
  swap-back, before Stage 1 runs), taken *before* calling into the reload
  path. Independent of `original_fundamentals[]`, own lifetime.
- `decompose_active` — a session flag mirroring `autoprunemono_active`
  (`reconcile.c:1506`, `main.c:2907`).
- `reconstruct` and `criterion=recon` (`get_recon_score`) must check
  `decompose_active` the same way they already check `autoprunemono_active`
  (`reconcile.c:1504-1520`, swap-back at `reconcile.c:2087-2095`) and use the
  pristine multi-copy trees instead of the decomposed fragments — reconciling
  duplication/loss against already-decomposed single-copy fragments would
  silently report near-zero duplications everywhere, which is wrong, not just
  imprecise. **Difference from the existing mechanism**: `autoprunemono`'s
  swap is element-by-element (same tree count, same indices, just two strings
  swapped per index). This can't work here because the tree count differs.
  Treat it instead as swapping two *whole dataset snapshots*: stash the
  current (decomposed) `fundamentals`/`Total_fund_trees`/`presence_of_taxa`/etc.
  aside, install the `pre_decompose_*` snapshot in their place, run
  `reconstruct`/`recon`, then swap back. Print an equivalent message to the
  existing `"(Using original unpruned trees for reconstruct)"` one so the
  behavior is visible, not silent.
- Internal callers should treat `decompose_active` as an idempotency guard
  too: if it's already set for the current in-memory dataset, a second
  request to decompose (e.g. from a future criterion's own auto-trigger) is a
  no-op, not a re-run.

### 6.4 Restoring the true original — no new "restore" command needed

Decomposition never modifies the file the user originally loaded from
(`inputfilename`, set at `main.c:2858`) — it only reads it and writes *new*
output files. So the full, session-independent way back to the
pre-decomposition (and pre-`autoprunemono`) state is simply:

```
exe <inputfilename>
```

which re-triggers the exact same full load/registration path from scratch.
This works even in a brand new process, doesn't depend on
`pre_decompose_fundamentals[]` still being alive, and needs no new command.
Print this as guidance (e.g. "To restore the original gene trees, run: exe
<inputfilename>") after `autodecompose` commits a change, rather than
building a bespoke `restoregenetrees` command that would just duplicate `exe`.
The in-memory `pre_decompose_*` snapshot from §6.3 exists for a different,
narrower purpose — letting `reconstruct`/`recon` work correctly *within the
same session* without forcing a disk round-trip — not as the primary restore
mechanism.

### 6.5 What a future coalescent criterion's internal auto-trigger can skip

If a future `criterion=mdc` (or `qfit`, detecting multi-copy input at search
start) wants to invoke this automatically rather than relying on
`autodecompose` having already run at load time:

- **Must run once per analysis setup, never per candidate tree.** This is a
  preprocessing step that changes the input pool, not a per-SPR-move rescore
  like `criterion=recon`. Once `fundamentals[]` holds the decomposed
  fragments, `hs`/`qfit`/`rf`/etc. run exactly as they do today — nothing
  about the search loop changes.
- **Skip re-running if `decompose_active` is already set** (§6.3) — check
  first, don't decompose twice.
- **Skip the CLI/option-parsing layer entirely.** An internal caller should
  call the Stage 1/2/3 functions (and `execute_command()`/
  `clann_load_trees()` for the reload) directly with explicit arguments, not
  build a command string and route it through `parsed_command[]`. This is
  exactly why the Stage 1 refactor in §5 (extracting `prune_monophylies`'s
  per-tree logic out of its command wrapper) matters — the reusable logic,
  not the CLI command, is the thing a future criterion links against.
- **The taxon-set hard-error validation for `speciestree=<file>` is simply
  not exercised** in the internal path, since there's normally no
  user-supplied file — nothing to explicitly skip, it's just inactive.
- **The `nj()` guide-tree fallback will rarely trigger** in this context,
  since by the time a coalescent criterion runs mid-session there's usually
  already something in `retained_supers[]` from an earlier command. It must
  still be implemented (fresh single-command sessions do hit it), just don't
  expect it to be the common path.
- **Do not skip writing `<filename>_info.txt`.** Keep the decision log even
  for automatic invocation — a user running `hs criterion=qfit` on multi-copy
  data should still be able to see what happened to their gene trees, same
  rationale as §5's "Output" subsection.

---

## 7. Open questions — decide during implementation, don't guess silently

- **What exactly gates "well-supported duplication call"?** Candidates: (a)
  reuse the existing `bsweight`/bootstrap-support infrastructure
  (`scoring.h:34`) if bootstrap values are present on the gene tree at node D;
  (b) a minimum branch-length separating D from its parent (near-zero branch
  → likely gene-tree noise); (c) both, with (a) preferred when available and
  (b) as fallback. Pick a conservative default that biases toward *not*
  cutting when evidence is ambiguous — under-splitting loses less information
  than over-fragmenting on noise.
- **`dupsupport` default value and units** depend on the answer above (a
  bootstrap percentage if (a), a branch-length threshold if (b)). Do not
  invent a number without at least running it against
  `examples/tutorial_multicopy.ph` (§7) and a larger real dataset to sanity
  check it doesn't either cut everything or nothing.
- **Fragment weighting scheme**: start with equal-split (`1/k`) as specified
  in Stage 3. A size- or informativeness-proportional scheme (e.g. weight by
  number of quartets a fragment can generate, normalized so the family still
  sums to 1) is a plausible future refinement but out of scope for the first
  implementation — don't build it speculatively.
- **Where `autodecompose` triggers**: at load time (simplest, matches
  `autoprunemono` precedent) vs. lazily the first time `hs criterion=qfit` (or
  a future `criterion=mdc`) runs on data still containing multi-copy trees.
  Default to load time unless that proves awkward during implementation.

---

## 8. Test fixtures already in the repo

`examples/tutorial_multicopy.ph` (8 trees) already contains the relevant edge
cases without needing new test data:

- Lines 1–3: ordinary single-copy trees — should pass through both stages
  completely unchanged, weight 1.
- Line 4 `(((Human.alpha,Human.beta),(Chimp.a,Chimp.b)),Gorilla),Orangutan),Macaque)`:
  clean monophyletic in-paralog pairs for both Human and Chimp — Stage 1
  (`prune_monophylies` logic) should collapse both pairs and the tree should
  become fully single-copy without ever reaching Stage 2.
- Line 5 `((((Human.x,Human.y),(Chimp.x,Chimp.y)),(Gorilla.a,Gorilla.b)),((Mouse.a,Mouse.b),(Rat.a,Rat.b)))`:
  same pattern across more taxa — same expected outcome as line 4.
- Line 6 `(((Human.1,Human.2),(Chimp.1,Chimp.2)),((Gorilla.1,Gorilla.2),(Orangutan.1,Orangutan.2)))`:
  same again.
- Line 7 `((((Human.alpha,Chimp),(Human.beta,Gorilla)),Orangutan),Macaque)`:
  **the important one** — `Human.alpha` and `Human.beta` are *not* sisters
  (Human.alpha pairs with Chimp, Human.beta pairs with Gorilla). Stage 1's
  monophyly test will not fire (neither subtree is purely Human). This is a
  genuine test of Stage 2's out-paralog cut logic against a guide tree — use
  it to verify the support gate and the cut/merge decision actually engage.
- Line 8 `(((Human.a,(Chimp,Gorilla)),(Human.b,(Orangutan,Macaque))),(Mouse,Rat))`:
  same category as line 7, different topology — a second Stage 2 test case.

`TUTORIAL.md`'s existing "Part 5: Multicopy gene trees and autoprunemono"
section (`TUTORIAL.md:232`) documents `autoprunemono` against this same file
and is the template to extend with a new "Part 6" (or subsection) documenting
`decomposegenetrees`/`autodecompose` once built, including a worked example
against line 7 or 8 above showing the `_info.txt` output.

---

## 9. Suggested implementation order

1. Refactor `prune_monophylies`'s per-tree collapse logic (`prune.c:238-343`)
   into a reusable function taking a tree and returning a collapsed tree + log
   entries, called by both the existing command and the new one. Verify
   `prunemonophylies` still behaves identically after the refactor (run it
   against `examples/tutorial_multicopy.ph`, diff output against the current
   `prunedtrees.txt`/`prunedtrees_info.txt`).
2. Implement guide-tree resolution (§4) as its own function, including the
   taxon-set validation hard-error for explicit files. Unit-test this in
   isolation (missing taxon, extra taxon, exact match, `memory` with/without
   `trees_in_memory`, nothing supplied and nothing in memory → falls through
   to `nj()`).
3. Implement Stage 2's cut/collapse/merge decision procedure as a function
   operating on an already-LCA-mapped-and-tagged gene tree (i.e. build on top
   of existing `label_gene_tree`/`reconstruct_map`, don't touch those
   functions). Get this correct against lines 7–8 of
   `examples/tutorial_multicopy.ph` before wiring up anything else.
4. Implement Stage 3 weighting and pool integration (adding fragments into
   `fundamentals[]`/`presence_of_taxa[][]`/`tree_weights[]`).
5. Wire up the `decomposegenetrees` command (option parsing, output files,
   help text) following the `prunemonophylies` pattern exactly
   (`main.c:1232-1240`, `1563`, `2013`).
6. Wire up `autodecompose=yes|no` following the `autoprunemono` pattern
   exactly (`main.c:65,152,2215,2811,3247`), including the §6.3
   `pre_decompose_*`/`decompose_active` snapshot taken *before* the §6.2
   reload, and the §6.3 changes to `reconstruct`/`get_recon_score` to swap to
   that snapshot instead of `original_fundamentals[]`-style element swapping.
7. End-to-end test: `clann exe examples/tutorial_multicopy.ph
   autodecompose=yes` then `hs criterion=qfit`, confirm all 8 trees (or their
   decomposed fragments) now contribute instead of only the 3 that are
   single-copy today. Also verify: `reconstruct` still reports the original
   duplication/loss counts (not near-zero) after `autodecompose=yes`, and
   `exe examples/tutorial_multicopy.ph` (no options) fully restores the
   pre-decomposition state.
8. Documentation: extend `TUTORIAL.md` Part 5 or add Part 6; add
   `decomposegenetrees`/`autodecompose` to the option-listing help text blocks
   in `main.c` (search for where `autoprunemono` appears in help text and add
   the new option alongside each occurrence, e.g. `main.c:1585,1702,1811`).

---

## 10. Explicitly out of scope for this work

- Implementing the deep-coalescence/MDC criterion itself (`criterion=mdc`) —
  this note is only about producing clean ortholog-subtree input for
  *existing* criteria (`qfit` today, `mdc` if it's built later).
- Any change to `qfit`/`rf`/`recon`/`ml` scoring code — decomposition only
  changes what goes into `fundamentals[]`, nothing about how those functions
  score it.
- Size- or informativeness-proportional fragment weighting (see §6) — ship
  equal-split first.
- Iterative refinement (decompose → build supertree → re-decompose against
  the improved guide tree → repeat) — worth doing later once the single-pass
  version works and its output has been sanity-checked on real data.
