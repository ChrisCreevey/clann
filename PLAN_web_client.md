# Plan — a web client + headless Clann server

**Goal.** Let a user drive Clann entirely from a browser: pick commands from
menus, upload/choose input files, set every command option through
command-specific forms, run the analysis, and see the results — including the
**interactive tree/reconciliation visualisation** already shipped in the HTML
viewer — without touching a terminal. Clann runs as a local **headless server**
on a port; the browser is the whole UI.

This document is written for **successive Claude instances**. Each step below is
scoped to be picked up cold, is small enough for one working session, and ends in
a **verifiable artefact** (a passing test, a `curl` transcript, a screenshot).
Read §1–§4 for orientation before starting any step in §6.

> Status legend in §6: `[ ]` not started · `[~]` in progress · `[x]` done.
> When you finish a step, tick it and note the verifying command/artefact.

---

## 1. The good news: how much already exists

Clann is **much closer to this than a standalone CLI would be.** Before planning
new work, know what is already in the tree:

| Asset | Where | What it gives us |
|-------|-------|------------------|
| **In-process library API** | `clann_api.c` / `clann_api.h` | `clann_init()`, `clann_load_trees()`, `clann_run_command()`, `clann_reset()`, `clann_set_output_fn()`. Runs any REPL command in-process, no subprocess. |
| **Output capture** | `printf2()` in `utils.c` (guarded by `CLANN_LIBRARY_MODE`) | Every line of Clann output can be redirected to a callback instead of stdout — the basis for returning results over HTTP. |
| **Non-fatal error unwinding** | `clann_jmp_exit` / `clann_exit_code`, `clean_exit()` | In library mode `clean_exit()` `longjmp`s back to the API frame instead of `_exit()`, so a Clann error does not kill the server. |
| **Shared library build** | `make libclann.so` (Makefile ~1015) | Builds `libclann.so` (a `.dylib` on macOS) with `-DCLANN_LIBRARY_MODE -fPIC`. Verified building today. |
| **Machine-readable option lists** | `opts_hs[]`, `opts_nj[]`, `opts_showtrees[]`, `opts_reconstruct[]`, … in `main.c` (~73–150) | The exact set of options each command accepts — the seed for auto-generated menus. |
| **Per-command help text** | `<command> ?` blocks in `main.c` | Human-readable option descriptions, values, and defaults — parseable into menu metadata. |
| **Interactive tree viewer** | `tools/clannview.template.html` → `viewer_template.h`; `html_view_*` in `reconcile.c` | The exact visualisation the client needs (navigator, reroot, collapse, search, phylogram/cladogram, SVG/Newick export, reconciliation events). Already emits a clean **JSON tree format** (see `NOTES_html_viewer.md` §3). |
| **Python binding** | `pyclann/` (`_core.py`, `_commands.py`, `_parser.py`) | A Python package that already models Clann commands, parses results, and returns typed `ClannResult` objects. Currently drives the **binary via subprocess**, but the command/result modelling is reusable. |

**Implication for architecture:** we do **not** need to write an HTTP server in
C, re-implement command parsing, or invent a tree-drawing UI. We need to (a) wrap
`libclann.so` behind a thin web server, (b) turn Clann's text output into
structured JSON, (c) expose option metadata, and (d) reuse the viewer. Most steps
are integration + hardening.

---

## 2. Recommended architecture

```
 ┌─────────────────────────── browser (the client) ──────────────────────────┐
 │  Single-page app:                                                          │
 │   • command palette / dropdown (populated from /api/commands)              │
 │   • per-command option form (populated from /api/commands/<name>/schema)   │
 │   • file picker / upload                                                   │
 │   • run button → job; live log stream; results panel                       │
 │   • tree visualisation = the existing clannview viewer, embedded           │
 └───────────────────────────────┬────────────────────────────────────────────┘
                                  │ HTTP + JSON (+ SSE/WebSocket for live logs)
 ┌────────────────────────────────▼───────────────────────────────────────────┐
 │  Web server  (Python: FastAPI/uvicorn — recommended)                        │
 │   • REST: /api/session, /api/files, /api/commands, /api/run, /api/jobs/<id> │
 │   • session registry (one Clann worker per session)                         │
 │   • job queue (async; long hs runs don't block the event loop)             │
 │   • serves the static client bundle                                         │
 └───────────────┬─────────────────────────────────────────────────────────────┘
                 │  in-process (ctypes/CFFI) OR one worker subprocess per session
 ┌───────────────▼─────────────────────────────────────────────────────────────┐
 │  libclann.so   (CLANN_LIBRARY_MODE)                                          │
 │   clann_init → clann_load_trees → clann_run_command … → clann_reset          │
 │   printf2 output → callback → captured per job                               │
 └──────────────────────────────────────────────────────────────────────────────┘
```

**Why a Python server wrapping `libclann.so`, not a C HTTP server:**

- Reuses `clann_api` + the `pyclann` command/result modelling directly.
- HTTP, JSON, SSE/WebSockets, file upload, sessions, and async job queues are
  trivial in FastAPI and painful in C.
- Keeps the concurrency/security-sensitive glue in a memory-safe language; the C
  core stays a single-threaded analysis engine.
- The server stays **local-first** (bind `127.0.0.1`), so it is a personal tool,
  not a public service — matching how a researcher would use it.

**The one hard constraint that shapes everything: Clann is built on massive
mutable global state and is not re-entrant.** `retained_supers`, `tree_top`,
`criterion`, `parsed_command`, the `threadprivate` OpenMP scratch — all global.
**Two analyses cannot run in one process at once.** This forces the concurrency
model (see §3.1) and is the single most important design fact in this plan.

### 2.1 Two concurrency models (pick in Step 1)

- **Model A — one worker *subprocess* per session.** The Python server spawns a
  small worker process per session; each worker loads `libclann.so` and owns one
  independent copy of Clann's global state. Requests for a session are routed to
  its worker and serialised there. **Pro:** true isolation, a crash or a runaway
  `hs` kills only that session, natural parallelism across sessions. **Con:** IPC
  plumbing (pipe/socket protocol between server and worker). **Recommended.**
- **Model B — single in-process engine, global lock.** The server holds one
  `libclann.so` instance and a global mutex; every analysis serialises. **Pro:**
  simplest to build first. **Con:** one slow `hs` blocks all users; a Clann
  crash takes the whole server down; no cross-session parallelism.

Start with **B for the first working slice** (Steps 1–4, single user), then
migrate to **A** at Step 8 when multi-session and robustness matter. The REST API
is identical for both, so the client never changes.

---

## 3. Major challenges & risk areas

Read these before estimating any step — several steps exist *only* to defuse
these.

### 3.1 Global state / re-entrancy (highest impact)
As above. Consequences: no concurrent analyses per process; `clann_reset()` must
reliably return the engine to a clean baseline between runs (it exists — verify
it in Step 2 with a load→run→reset→load→run loop that produces identical results
to two fresh processes).

### 3.2 Long-running, blocking commands
`hs` can run for minutes. The HTTP request cannot block that long. Need an
**async job model**: `POST /api/run` returns a `job_id` immediately; the client
polls `GET /api/jobs/<id>` or subscribes to a log stream. Requires:
- **Progress streaming.** `hs` already prints per-rep progress lines via
  `printf2`; capture them and forward as Server-Sent Events. (`hs progress=`
  controls cadence in the parallel path.)
- **Cancellation.** Clann has a `controlc` SIGINT handler that cleanly stops a
  search. In Model A, cancel = send SIGINT to the worker. In Model B, this is
  much harder (can't interrupt a blocking in-process call cleanly) — another
  reason A wins for production.

### 3.3 Text output, not structured data
`printf2` emits human ASCII (tables, tree drawings, "Supertree 1 of 1 score =
17.0000"). The client needs Newick/JSON + numbers. Two complementary tactics:
- **Reuse the viewer's JSON emitter.** `html_view_*` already serialise trees (and
  reconciliations) to the clean JSON in `NOTES_html_viewer.md` §3. Add a mode
  that writes that JSON to a **string/file the API can read**, rather than a full
  HTML file — this is the cleanest path to structured tree results.
- **Parse the rest from captured text**, reusing `pyclann/_commands.py` which
  already extracts trees/scores from Clann output. Treat this as the fallback for
  scalars until first-class JSON exists.
Risk: text parsing is brittle across commands. Prefer extending the JSON emitter
over growing the parser.

### 3.4 File input & sandboxing
The client uploads or selects tree files. Server must: write them to a per-session
**sandbox dir**, pass only that path to `clann_load_trees`, and **never** let a
command reference arbitrary filesystem paths. Clann commands take file paths
freely (`exe /any/path`, `savetrees=/any/path`, `htmlview=/any/path`), so the
server must validate/rewrite paths to stay inside the session sandbox.

### 3.5 Security surface (must-fix before binding anything but loopback)
- **Shell escapes in Clann.** `!` runs `system("bash")` (`main.c:1011`); several
  commands call PAUP\* via `system(system_call)`. `html_view_launch` calls the
  platform opener. **A server build must hard-disable the `!` command and any
  `system()` path**, ideally behind a `CLANN_SERVER_MODE`/`CLANN_NO_SHELL`
  compile guard or a runtime flag. Verify no `system()` is reachable from
  `clann_run_command` in server mode.
- **Path traversal** (see §3.4).
- **Bind to `127.0.0.1` by default.** Exposing to `0.0.0.0` must be an explicit,
  documented opt-in with a warning — there is no auth model.
- **Auto-open browser** already gates on `isatty`; in server mode disable
  `html_view_launch` entirely (the browser is the client, not spawned by Clann).

### 3.6 Option metadata is thin
`opts_*[]` gives option **names** only. Rich forms need **type** (bool / int /
float / enum / filename / string), **default**, **allowed values**, and **help
text**. Sources: parse the `<command> ?` help blocks, or hand-author a
`command_schema.json`. Hand-authoring is more reliable and is itself a verifiable
artefact; auto-parsing help is a nice-to-have. Keep the schema in the repo so the
CLI help and the web forms cannot drift (a test can assert every `opts_*[]` name
appears in the schema and vice-versa).

### 3.7 Statefulness
Clann is a stateful REPL: `set criterion=recon` then `hs`. The web client must
carry a **session** whose state = loaded trees + `set` options + trees in memory.
Map one browser session → one Clann engine/worker. Surface current state
(`/api/session` returns loaded file, taxa count, criterion, trees in memory) so
the UI can show it.

### 3.8 Build/packaging
`libclann.so` builds today, but: OpenMP linkage, macOS `.dylib` vs Linux `.so`
naming, and shipping the client bundle + Python deps. Package as a `pip install`
that vendors or builds the lib, plus a `clann-web` entry point that starts the
server and prints the URL.

---

## 4. What "done" looks like (acceptance)

A user runs `clann-web` (or `python -m clann_web`), a browser opens at
`http://127.0.0.1:PORT`, and they can: upload a `.ph` file; see it loaded (taxa,
tree count); choose `hs`, pick `criterion=recon`, set options via a form; run it
with a live progress log; and view the resulting supertree(s) in the interactive
viewer — then run `reconstruct` and view the reconciliations with duplication/loss
events. No terminal involved.

---

## 5. Testing philosophy

Every step ends in something runnable:
- Server steps → a `pytest` test hitting the endpoint with `httpx`, and a pasted
  `curl` transcript in the step's completion note.
- Engine/C steps → a small C or Python harness asserting behaviour (e.g. the
  reset-loop determinism test), plus the `17.0000` recon regression must still
  hold (`NOTES`/handover anchor).
- Client steps → a screenshot via the browser tools and a note on what was
  clicked. Reuse the pattern already used to verify the HTML viewer.

Keep a golden dataset: `examples/tutorial_multicopy.ph` (recon best = 17 legacy).

---

## 6. Step-by-step plan

Each step: **Goal · Work · Verify · Done-when.** Steps are ordered so each builds
on a green predecessor. A single Claude instance should take one step.

### Phase 0 — Foundations & de-risking

**Step 0.1 — Prove one live engine runs a persistent multi-command session.** `[x]`
- *Goal:* confirm the recommended design's load-bearing assumption — that a single
  in-process `libclann.so` engine holds state across many commands (the REPL feel
  over an API), **unlike** today's one-shot `pyclann` subprocess model (§1, §3.7).
- *Work / Verify (DONE):* a ctypes harness (`tools/session_persistence_check.py`,
  committed) `dlopen`s `libclann.so` **once** (after `make libclann.so`),
  registers an output callback via `clann_set_output_fn`,
  then in a **single process** runs the full user scenario with **no reload
  between commands**:
  `clann_load_trees(tutorial_multicopy.ph)` → `set criterion=recon` → `nj` →
  `hs nreps=5 seed=42` → `showtrees` → `reconstruct speciestree=memory`.
  Result (captured): `hs` reaches `17.000000`; **`showtrees` still reports
  "8 source trees" after `nj`+`hs`** (original trees survived in memory);
  `reconstruct speciestree=memory` consumed the best tree `hs` left behind
  (`total dup+loss score: 17.0000`). Harness prints `PASS`.
  - *Key facts confirmed:* state persists in Clann's globals between
    `clann_run_command` calls; `quit` is a no-op in library mode (session stays
    live until the host ends it); one session must **serialise** its commands
    (single engine, single global state).
  - *Build gotcha for future agents:* `make libclann.so` here produces an
    **x86_64** dylib (Homebrew `gcc-15`). Load it from a matching interpreter —
    e.g. `arch -x86_64 /usr/bin/python3` on Apple Silicon — or an arch mismatch
    (`incompatible architecture`) will fail the `dlopen`.
- *Done-when:* ✅ demonstrated. The persistent-session engine the web server needs
  already exists and works.

**Step 0.1b — Confirm `clann_reset()` returns to a clean baseline.** `[ ]`
- *Goal:* the *other* half of §3.1/§3.7 — that ending a session and starting a new
  one in the same process is deterministic (needed for worker reuse / "start
  fresh", not for the persistent path proven in 0.1).
- *Work:* extend the 0.1 harness: run the scenario, `clann_reset()`, run it again;
  assert the second run reproduces the first exactly and matches a fresh-process
  baseline (same `17.0000`, same best tree). Watch for state that `clean_exit()`
  frees but does not re-null (the API already zeroes a known set — verify no
  others leak).
- *Verify:* harness prints `PASS`; two post-reset runs are byte-identical on the
  score/best-tree lines.
- *Done-when:* reset determinism is demonstrated, or a specific leak is documented
  as a bug to fix before Model A worker reuse.

**Step 0.2 — Server-safe build guard (`CLANN_SERVER_MODE`).** `[ ]`
- *Goal:* a compile mode with no shell/`system()` reachable and no browser
  auto-open (defuses §3.5).
- *Work:* add `-DCLANN_SERVER_MODE`; guard `!`/`system("bash")`, the PAUP\*
  `system()` calls, and `html_view_launch` so they are compiled out or refuse at
  runtime with a clear message. Add a `libclann-server.so` make target.
- *Verify:* `grep`-driven test (or `nm`) shows no `system` call is reachable from
  `clann_run_command` in the server build; running `!` returns an error string,
  not a shell. Standalone binary behaviour unchanged.
- *Done-when:* server lib builds and the shell surface is provably gone.

**Step 0.3 — Structured tree output from the engine.** `[ ]`
- *Goal:* get the viewer's JSON (a tree/reconciliation as in `NOTES_html_viewer.md`
  §3) as a **string the API can capture**, not only as a full `.html` file.
- *Work:* factor the JSON-emitting core out of `html_view_*` so it can write to a
  memory buffer / temp file; add a tiny API (e.g. `clann_last_trees_json(buf)` or
  a `resultjson=<file>` option on `hs`/`nj`/`reconstruct`/`showtrees`).
- *Verify:* a harness runs `hs …` then reads back valid JSON whose leaf set equals
  the input taxa; reconciliation JSON contains `event` fields.
- *Done-when:* structured results exist without scraping ASCII.

### Phase 1 — Minimal vertical slice (single user, Model B)

**Step 1.1 — Skeleton FastAPI server wrapping the engine.** `[ ]`
- *Goal:* one process, global-lock engine, three endpoints end-to-end.
- *Work:* `clann_web/` Python package. Load `libclann.so` via ctypes (reuse
  `pyclann` patterns). Endpoints: `POST /api/session` (init/reset engine),
  `POST /api/load` (path to a bundled example), `POST /api/run` (a command
  string, **synchronous** for now), returning captured text. Bind `127.0.0.1`.
- *Verify:* `pytest` + `curl`: create session → load tutorial → `run "set
  criterion=recon"` → `run "hs nreps=5 seed=42"` → response text contains
  `17.0000`.
- *Done-when:* a browserless HTTP round-trip drives a real analysis.

**Step 1.2 — File upload into a session sandbox.** `[ ]`
- *Goal:* replace "bundled example" with real uploads, safely (§3.4).
- *Work:* `POST /api/files` (multipart) writes to a per-session temp dir; `load`
  accepts only a sandbox-relative filename; reject/relocate absolute or `..`
  paths. Same rewrite for any output path option.
- *Verify:* `pytest`: upload `tutorial_multicopy.ph`, load it, `showtrees` reports
  8 trees; an upload named `../evil` is rejected; a `savetrees=/etc/x` is confined
  to the sandbox.
- *Done-when:* arbitrary paths cannot escape the session dir.

**Step 1.3 — Structured run result.** `[ ]`
- *Goal:* `/api/run` returns JSON `{ ok, log, trees[], scores[], session_state }`,
  not raw text.
- *Work:* use Step 0.3's JSON for trees; reuse `pyclann/_commands.py` parsing for
  scalars; add `/api/session` returning loaded file, taxa count, criterion, trees
  in memory.
- *Verify:* `pytest`: after `hs`, `trees[0]` is valid Newick over the input taxa
  and `scores[0]==17.0` (recon legacy); `session_state.criterion=="recon"`.
- *Done-when:* the client can consume results without parsing prose.

### Phase 2 — The browser client

**Step 2.1 — Static SPA shell + command palette.** `[ ]`
- *Goal:* a served single page that lists commands and shows a session panel.
- *Work:* `GET /api/commands` returns command names (+ one-line blurbs) derived
  from the `opts_*[]`/help inventory. Client renders a dropdown and a session
  status area. Serve the bundle from the server.
- *Verify:* browser screenshot: page loads at `127.0.0.1:PORT`, dropdown lists
  `hs`, `nj`, `showtrees`, `reconstruct`, …; session panel shows "no file loaded".
- *Done-when:* the shell renders and talks to `/api/commands`.

**Step 2.2 — Command schema → dynamic option form.** `[ ]`
- *Goal:* per-command forms with typed inputs, defaults, and help (§3.6).
- *Work:* author `command_schema.json` (name, options[type,default,values,help]);
  `GET /api/commands/<name>/schema`; client renders bool→checkbox, enum→select,
  int/float→number, filename→file picker, with defaults + tooltips. Add a test
  asserting schema ⇄ `opts_*[]` consistency.
- *Verify:* screenshot: selecting `hs criterion=recon` shows `lossmodel`
  (legacy/standard), `numspeciesrootings`, `htmlview`, etc. with correct widgets;
  consistency test passes.
- *Done-when:* every command's options are editable from a generated form.

**Step 2.3 — Upload + load + run from the UI.** `[ ]`
- *Goal:* full loop without a terminal (synchronous run is fine here).
- *Work:* wire file picker → `/api/files`; "Load" → `/api/load`; "Run" → `/api/run`
  with form values; render the returned `log` and scores.
- *Verify:* screenshot sequence: upload tutorial → load (panel shows 8 trees,
  16 taxa) → run `hs` → score `17` shown.
- *Done-when:* a user can complete an analysis by clicking.

**Step 2.4 — Embed the interactive tree viewer.** `[ ]`
- *Goal:* results show the existing clannview visualisation, fed by API JSON.
- *Work:* refactor `tools/clannview.template.html` so its JS can be mounted in the
  SPA and **loaded with JSON from the API** (not only from a baked `.html`). Route
  `trees[]` from `/api/run` into it. Reuse navigator/reroot/collapse/search as-is.
- *Verify:* screenshot: after `hs`, the supertree renders in the embedded viewer;
  after `reconstruct`, reconciliations render with duplication/loss events and the
  tree navigator (matches the standalone viewer already verified).
- *Done-when:* the visualisation the user asked for is live in the client.

### Phase 3 — Long-running jobs & robustness

**Step 3.1 — Async job model + live log stream.** `[ ]`
- *Goal:* non-blocking `hs`; progress in the UI (§3.2).
- *Work:* `/api/run` returns `job_id`; run the engine in a worker thread/process;
  stream `printf2` lines over SSE (`/api/jobs/<id>/stream`); `GET /api/jobs/<id>`
  for status/result. Client shows a live log and a spinner.
- *Verify:* browser + `curl`: a long `hs nreps=50` streams per-rep lines and
  finishes with a result; the UI stays responsive.
- *Done-when:* long runs report progress and don't block the server.

**Step 3.2 — Cancellation.** `[ ]`
- *Goal:* stop a running search from the UI.
- *Work:* `POST /api/jobs/<id>/cancel`. (Model A: SIGINT to the worker via the
  existing `controlc` handler.) A "Stop" button in the UI.
- *Verify:* start `hs nreps=100`, cancel mid-run; server returns partial/best-so-
  far cleanly and accepts the next command.
- *Done-when:* runaway searches are interruptible.

**Step 3.3 — Migrate to per-session worker processes (Model A).** `[ ]`
- *Goal:* isolation, cross-session parallelism, crash containment (§2.1, §3.1).
- *Work:* replace the in-process global-lock engine with a worker-process pool,
  one per session, speaking a small line/JSON protocol (load/run/reset/cancel).
  Server routes by session id. API unchanged.
- *Verify:* two sessions run `hs` **concurrently** with independent state; killing
  one worker leaves the other healthy; the Phase 1–3 tests still pass unchanged.
- *Done-when:* multi-user isolation holds and nothing above regressed.

### Phase 4 — Hardening & packaging

**Step 4.1 — Security pass.** `[ ]`
- *Goal:* safe defaults before anyone binds beyond loopback (§3.5).
- *Work:* confirm server build has no shell reachability; enforce sandbox on every
  path option; loopback-only default with an explicit `--host 0.0.0.0` warning;
  request size limits; per-session temp cleanup.
- *Verify:* a `pytest` security suite: `!` refused; path traversal blocked; output
  paths confined; default bind is `127.0.0.1`.
- *Done-when:* the security suite is green.

**Step 4.2 — Packaging & one-command launch.** `[ ]`
- *Goal:* `pip install` + `clann-web` opens the app.
- *Work:* build/vendor the server lib per-platform (`.so`/`.dylib`), ship the
  client bundle + schema, add a `clann-web` entry point that starts uvicorn and
  prints/opens the URL. Document in `USER_MANUAL.md` and a new
  `NOTES_web_client.md`.
- *Verify:* in a clean venv, `pip install .` then `clann-web` serves the app;
  smoke test drives one full analysis.
- *Done-when:* a non-developer can install and use it.

**Step 4.3 — Command coverage & polish.** `[ ]`
- *Goal:* extend beyond the core four commands.
- *Work:* add schema + result handling for `bootstrap`, `consensus`, `alltrees`,
  `usertrees`, `rfdists`, `mlscores`, `excludetrees`/`includetrees`, etc.; surface
  output files (histograms, landscape TSV) for download; show `set` state.
- *Verify:* each added command has a passing round-trip test and a screenshot.
- *Done-when:* the common workflows are all clickable.

---

## 7. Open decisions for the user (resolve before Phase 1)

1. **Server language** — Python/FastAPI (recommended, reuses `pyclann`) vs Node vs
   embedded C HTTP server. Affects every server step.
2. **Concurrency model start point** — begin at Model B then migrate (recommended)
   vs build Model A up front.
3. **Structured tree output** — extend the viewer's JSON emitter in C (cleaner,
   recommended) vs parse ASCII in Python (faster to start, brittle).
4. **Scope of v1 commands** — just `hs`/`nj`/`showtrees`/`reconstruct`, or the full
   command set from the start.
5. **Distribution** — local-only personal tool (recommended) vs a shared/hosted
   deployment (needs auth, real sandboxing, resource limits — out of scope here).

---

## 8. Suggested next three sessions

**Step 0.1 is already done** (§6) — the persistent in-process engine is proven, so
the core premise of the whole design is sound. Next:

1. **Step 0.1b** — reset-determinism (small; extends the existing 0.1 harness;
   confirms the "start fresh" / worker-reuse path).
2. **Step 0.2** — `CLANN_SERVER_MODE` shell lockdown (unblocks binding a port).
3. **Step 1.1** — skeleton FastAPI server driving a real, *persistent* session
   over HTTP (load once, then `nj`/`hs`/`reconstruct` against the same engine).

After those, the architecture is proven end-to-end (engine ↔ HTTP ↔ real
multi-command session) and the remaining steps are incremental UI + robustness.
