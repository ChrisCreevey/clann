# Design brief — Clann web client

A brief for a designer (e.g. Claude) to **restyle the web interface** of Clann.
The goal is purely visual/UX polish of an existing, working single-page app — not
new features and not changes to behaviour or the server.

---

## 1. What Clann is

**Clann** is a command-line tool used by biologists (phylogeneticists) to build
**supertrees** — single evolutionary trees assembled from many smaller,
conflicting gene trees — and to **reconcile** gene trees against a species tree to
infer gene duplication and loss events. It is a serious research instrument: the
people using it are scientists analysing real data, often iterating (load trees →
pick a scoring criterion → run a search → inspect the resulting tree → try again).

The output that matters most is a **tree diagram**, shown in an interactive viewer
(pan/zoom, reroot, cladogram vs phylogram, export). Trees can also be expressed as
**Newick** text, e.g. `((Human,Chimp),Gorilla);`.

## 2. What the web client is

Historically Clann is terminal-only. This web client lets a researcher drive it
**entirely from a browser** on their own machine:

1. **Upload** a source-tree file (`.ph`).
2. **Load** it — the session panel shows taxa count, number of source trees.
3. **Set** a scoring criterion (e.g. `recon` for reconciliation).
4. **Run** a command (e.g. a heuristic search `hs`) via a typed option form, with a
   **live streaming log** and a **Stop** button.
5. **View** the resulting tree in the embedded interactive viewer.

It is a **local, single-user, offline tool**. The server binds to `127.0.0.1`
(loopback) with no login; the browser *is* the whole UI. Tone should feel like a
**precise, calm scientific instrument** — trustworthy, legible, uncluttered — not a
consumer marketing site. Think "well-made lab tool / good developer tool," not
"landing page."

## 3. What we want from you

Finesse the look and feel of the interface: **typography, spacing, colour,
hierarchy, card and control styling, states (loading / streaming / empty / error /
cancelled), micro-interactions, and responsiveness.** Make it feel considered and
professional in **both light and dark mode**. Improve clarity of the run flow and
the results presentation. You have wide latitude on aesthetics as long as the hard
constraints in §5 hold.

Concrete things worth improving (suggestions, not a spec):
- A clearer visual through-line for the **three-step flow** (session/data → run →
  results) so a first-time user knows where to start.
- Better **states**: connecting vs connected; a run that is streaming; an empty
  "no results yet"; an error; a cancelled run.
- A more refined **results area** — the scores summary line, the Newick table, the
  collapsible log, and the framing around the embedded tree viewer `<iframe>`.
- Considered **buttons and form controls** (the option form is generated
  dynamically from a schema — see §6).
- Tasteful use of the `Clann` wordmark / header.

## 4. Current interface at a glance

Single page, served at `/`. Layout:

- **Header** — `Clann` wordmark + a connection status chip (`connecting…` →
  `connected`).
- **Two cards side by side** (a responsive 2-column grid, ~1000px max width):
  - **Session & data** — a key/value readout (Input file, Taxa, Source trees,
    Criterion, In memory), a file upload row, a "loaded file" select + Load, a
    "New session" button.
  - **Run a command** — a command `<select>`, a one-line description, a
    **dynamically generated typed option form**, a collapsible "Advanced (raw
    options)" text field, and **Run** / **Stop** buttons.
- **Results card** (full width, hidden until the first run) — a scores/summary
  line, the **embedded interactive tree viewer** in an `<iframe>`, and two
  collapsible `<details>`: "Trees (Newick)" (a table) and "Log" (a dark streaming
  console).
- **Footer** — a one-line "local session · loopback · per-session worker" note.

The current styling is a competent-but-plain baseline: system font, CSS custom
properties for colour, light/dark via `prefers-color-scheme`, rounded cards. It's
a fine starting point to elevate.

## 5. Hard constraints (do not break these)

1. **Single self-contained file.** Everything lives in one `index.html` with
   **inline `<style>` and inline `<script>`**. Keep it that way: **no external
   stylesheets, no CDN scripts, no web-font URLs, no remote images.** It runs
   offline on a researcher's laptop. System font stacks and inline SVG/data-URIs
   are fine.
2. **Preserve every element `id` the JavaScript uses** (see §7). The script does
   `getElementById(...)` and `querySelectorAll('[data-name]')` extensively. You may
   restructure the markup and restyle freely, but these IDs must survive and keep
   their roles (input/select/button/container as noted).
3. **Preserve the class names and CSS-variable names the code references**
   (see §7). Some are read from JavaScript and inline styles, so renaming them
   silently breaks behaviour.
4. **Keep both light and dark themes.** Dark mode is not optional; many users work
   in it. Drive it the way it's done today (CSS variables + `prefers-color-scheme`),
   or better, but both must look intentional.
5. **Behaviour, DOM-population, and the server stay untouched.** The JS builds the
   command form, streams the log, toggles the viewer, etc. Don't change what the
   script does — only how it looks. (If a change genuinely needs a markup hook, add
   a class or `data-` attribute rather than removing an existing one.)
6. **Accessible & responsive.** Keep labels associated with inputs, maintain colour
   contrast in both themes, and make the two-column layout collapse gracefully on
   narrow screens. The page body must not scroll horizontally; wide content (the
   Newick table, the log) should scroll within its own container.

## 6. Behaviour you should understand (so the states you design are real)

- **Runs are asynchronous.** Clicking **Run** posts the command, then the UI opens
  a **Server-Sent Events** stream and appends log lines **live** while the search
  runs (this can take seconds to minutes). Design for a genuinely *streaming*,
  growing log — not a spinner-then-dump.
- **Stop** appears only while a run is active. Cancelling **resets the session**
  (loaded trees cleared; uploaded files survive) and shows an informational
  message — design that recovery state.
- **The option form is generated from a schema.** Each command
  (`hs`, `nj`, `reconstruct`, `set`, …) returns typed options — enum
  (`<select>`), int/float (`number`), text — each with a default and a one-line
  help hint. The form is a 2-column grid of `.fitem`s. Some commands have ~12
  options; some have none. See `clann_web/command_schema.json` for the real
  content, and `clann_web/commands.py` for the command list + blurbs.
- **Results shapes.** A run returns `scores` and `trees`. Example realistic data:
  - session `state`: `{input_file:"tutorial_multicopy.ph", num_taxa:9,
    num_source_trees:8, criterion:"recon", trees_in_memory:1}`
  - a supertree result: `scores:[17.0]`, one tree row with a Newick string.
  - a reconciliation result: 8 trees, each with a `score` and duplication/loss
    counts shown as `(d2/l5)`.
- **The tree viewer is in scope.** It is a separate self-contained HTML document
  (`tools/clannview.template.html`) loaded into the `#viewer` `<iframe>` from
  `/api/viewer` — it renders the interactive tree (pan/zoom, reroot, cladogram vs
  phylogram, row-spacing/font sliders, SVG/Newick export, reconciliation
  duplication/loss events). Please finesse **both** the frame around the iframe
  **and** the viewer document itself so the two read as one designed system
  (shared palette, type, control styling, light/dark). The viewer is the emotional
  centre of the tool — it's what the researcher actually came to see — so it's
  worth real attention. It has its own inline `<style>`/`<script>` and its own
  controls; keep it self-contained and don't change what it computes, only how it
  looks. It also has a data-injection marker pair
  (`/*CLANN_DATA_BEGIN*/ … /*CLANN_DATA_END*/`) that the server fills at runtime —
  leave that mechanism intact.

## 7. The contract: names the code depends on

**Element IDs (must exist, keep their role):**

| id | role |
|----|------|
| `conn` | status chip; text set to `connecting…`/`connected`, colour set via `var(--ok)` |
| `s-file`,`s-taxa`,`s-trees`,`s-crit`,`s-mem` | session key/value readouts |
| `file` | `<input type=file>` · `upload` | Upload button |
| `files` | `<select>` of uploaded files · `load` | Load button |
| `new-session` | button · `s-id` | session-id text |
| `data-msg` | status/error line for the session card |
| `cmd` | command `<select>` · `cmd-blurb` | its description line |
| `form` | container the typed option form is generated into |
| `opts` | "Advanced (raw options)" `<input type=text>` |
| `run` | primary Run button · `stop` | Stop button (hidden unless running) |
| `run-cmd` | shows the built `clann> …` command · `run-msg` | run status line |
| `results-card` | full-width results section (hidden until a result) |
| `scores` | scores/summary line |
| `viewer` | `<iframe>` for the interactive tree viewer |
| `trees-wrap` | container for the generated Newick `<table>` |
| `logbox` | `<details>` wrapping the log · `log` | `<pre>` the log streams into |

**Class names referenced from CSS/JS (keep or restyle, don't drop):**
`card`, `full`, `kv`/`dt`/`dd`, `fld`, `blurb`, `btn`, `btn.primary`, `chip`,
`muted`, `row`, `msg` + `msg.ok`/`msg.err`, `form-grid`, `fitem`, `fitem .name`,
`fitem .hint`, `log`, `trees` (+ `th`/`td`/`td.newick`).

**CSS custom properties referenced from JS/inline styles (keep these token names;
add more freely):** `--ok` (used for the connected chip), `--bg` and `--line`
(used inline on the viewer iframe), plus the palette in general: `--panel`,
`--ink`, `--muted`, `--accent`, `--accent-ink`, `--chip`, `--warn`, `--err`,
`--logbg`, `--logink`.

## 8. How to preview

**Easiest — the static mock.** `clann_web/static/index.mock.html` is a
**data-populated standalone copy** of the app: open it directly in a browser (no
server needed) and every state is visible. It has a **mock-only state switcher**
(top-right: Idle / Running / Results / Error / Cancelled) so you can see the run
flow, the streaming-log and Stop states, the results table, and a representative
tree in the viewer frame. Use it as your working canvas for the layout/CSS, then
port the styling into the real `index.html` (and the viewer template).

> Everything in the mock marked `data-mock` (the switcher and its script) plus the
> `.mock-viewer` stand-in are **preview-only** — they must not end up in the real
> `index.html`. In the real app the viewer is a live `<iframe id="viewer">`.

**The tree viewer, populated.** `tools/clannview.mock.html` is a data-populated
standalone render of the viewer (`tools/clannview.template.html`) — open it
directly to design the viewer with real content: a 4-tree reconciliation dataset
exercising every event type (speciation dots, red duplication squares, dashed
loss branches), plus support values, branch lengths, the tree navigator, and all
the controls. It is a faithful render; the only baked-in difference is the
`const DATA = {…}` object the server fills at runtime. **Do your styling in
`clannview.template.html` and re-preview via the mock** — don't ship edits made
only to the mock.

**Full fidelity — the real app.** To see live data and the real interactive
viewer, run the server:

```bash
# from the repo root (see clann_web/README.md for the one-time lib build)
./venv/bin/clann-web            # or: PYTHONPATH="$PWD" arch -x86_64 /usr/bin/python3 -m clann_web
# open http://127.0.0.1:8765, upload examples/tutorial_multicopy.ph, Load,
# run `set criterion=recon`, then `hs` to see a full results state + live viewer.
```

## 9. Files to hand the designer

**Essential (the files being restyled):**
- `clann_web/static/index.html` — the single-page app: markup, inline CSS, inline
  JS. **The main file to edit.**
- `tools/clannview.template.html` — the interactive tree viewer shown in the
  results `<iframe>` (in scope — see §6). Its own self-contained page to finesse.

**Preview / canvas (open directly, no server):**
- `clann_web/static/index.mock.html` — data-populated standalone mock of the SPA
  with a state switcher. Preview canvas for `index.html` (see §8).
- `tools/clannview.mock.html` — data-populated standalone render of the tree
  viewer. Preview canvas for `clannview.template.html` (see §8).

**For understanding (read-only context):**
- `clann_web/DESIGN_BRIEF.md` — this document.
- `clann_web/README.md` — what the tool does and its endpoints.
- `clann_web/command_schema.json` — the real option content the form renders.
- `clann_web/commands.py` — the command list + one-line blurbs shown in the UI.
- `clann_web/server.py` — optional: how `/` and `/api/viewer` are served and the
  API response shapes.

**Do NOT need to touch:** the Python server logic, the engine/worker code, tests,
or the build. This is a front-end styling pass on the two HTML files.
