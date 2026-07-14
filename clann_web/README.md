# clann_web — a local browser client for Clann

Run Clann from a web browser: upload source-tree files, pick commands from
menus, set options, run analyses, and explore the results — including the
interactive tree/reconciliation viewer — without touching a terminal.

It is a **local, single-user tool**: the server binds to loopback (`127.0.0.1`)
and drives one in-process Clann engine. See
[`../PLAN_web_client.md`](../PLAN_web_client.md) for the design and roadmap.

---

## 1. One-time setup: build the server library

The web server talks to a hardened build of Clann — `libclann-server.so` — that
has **no shell/`system()` surface** (safe to sit behind a port). Build it once
from the repository root:

```bash
cd /path/to/clann
make libclann-server.so
```

## 2. Start the server

From the repository root:

```bash
PYTHONPATH="$PWD" arch -x86_64 /usr/bin/python3 -m clann_web
```

You'll see:

```
clann-web serving on http://127.0.0.1:8765  (session …)
```

Open **http://127.0.0.1:8765** in your browser. Stop the server with `Ctrl+C`.

Options:

```bash
# choose a different port
PYTHONPATH="$PWD" arch -x86_64 /usr/bin/python3 -m clann_web --port 9000
```

### Why the `arch -x86_64 /usr/bin/python3` prefix?

On Apple Silicon, the Homebrew `gcc-15` build produces an **x86_64**
`libclann-server.so`, so it must be loaded by a matching (x86_64) Python.
`/usr/bin/python3` ships an x86_64 slice that runs under Rosetta;
`arch -x86_64` selects it. A native-arm64 `python3` will fail to `dlopen` the
library with an *"incompatible architecture"* error.

On Linux (or once a native/universal library is built), just use your normal
Python:

```bash
PYTHONPATH="$PWD" python3 -m clann_web
```

If the library is somewhere non-default, point at it with `CLANN_LIB`:

```bash
CLANN_LIB=/path/to/libclann-server.so PYTHONPATH="$PWD" python3 -m clann_web
```

---

## 3. Using it in the browser

1. **Upload** — *Choose file*, pick a source-tree file, then **Upload**. Bundled
   examples live in the repo's `examples/` folder. Good starting points:
   - `examples/tutorial_multicopy.ph` — 8 gene trees, several with duplicated
     copies (best for seeing reconciliations)
   - `examples/tutorial_single.ph` — single-copy trees
   - `examples/tutorial_trees.ph`
2. **Load** — click **Load**; the session panel updates (e.g. *8 trees, 9 taxa*).
3. **Set the criterion** — choose `set`, put `criterion=recon` in options, **Run**.
   Criteria: `dfit` (default), `sfit`, `qfit`, `rf`, `ml`, `recon`.
4. *(optional, for reproducibility)* run `set` with `seed=42`.
5. **Search** — choose `hs`, options e.g. `nreps=10 nthreads=4`, **Run**. The best
   supertree renders in the interactive viewer (pan/zoom, reroot, cladogram /
   phylogram, SVG / Newick export).
6. **Reconcile** — choose `reconstruct`, options `speciestree=memory`, **Run**. The
   viewer's tree navigator now shows every reconciled gene tree with
   duplication/loss events.

**New session** clears everything (loaded trees, criterion, results).

### Quick check without a browser

```bash
curl -s -X POST http://127.0.0.1:8765/api/session | python3 -m json.tool
```

Returns a `session_id` and an empty `state`.

---

## 4. Known limitations (being worked on)

- **`hs seed=` is silently ignored** — only the `set seed=` command seeds the RNG.
  For reproducible searches, run `set seed=42` *before* `hs` (not `hs seed=42`).
- **Long searches block the page** — a large `hs` run holds the HTTP request until
  it finishes; there is no progress bar or cancel yet (planned: async jobs + live
  log streaming). Keep `nreps` modest while experimenting.
- **One analysis at a time** — the engine keeps all state in process globals, so a
  session runs commands one after another. (Per-session worker processes for
  concurrent users are planned.)
- **Loopback only** — there is no authentication. Do not expose the server beyond
  `127.0.0.1`. Binding another host prints a warning.

---

## 5. HTTP API (for scripting)

| Method & path            | Body / query                 | Returns |
|--------------------------|------------------------------|---------|
| `POST /api/session`      | —                            | `{session_id, state, files}` (resets to a clean session) |
| `GET  /api/session`      | —                            | `{session_id, state, files}` |
| `POST /api/files?name=F` | raw file bytes               | `{ok, stored, files}` (name must be a bare filename) |
| `POST /api/load`         | `{"file":"F"}`               | `{ok, log, state}` (F must be uploaded first) |
| `POST /api/run`          | `{"command":"hs nreps=10"}`  | `{ok, log, command, state, trees?, scores?, result_type?, has_viewer?}` |
| `GET  /api/commands`     | —                            | `{commands:[{name, blurb}]}` |
| `GET  /api/viewer`       | —                            | interactive viewer HTML for the latest tree result |
| `GET  /`                 | —                            | the single-page app |

`state` is `{input_file, num_taxa, num_source_trees, criterion, trees_in_memory}`.
For tree-producing commands (`hs`, `nj`, `showtrees`, `reconstruct`, `alltrees`),
`/api/run` also returns `trees` (each `{name, newick, tree, score?, dups?,
losses?}`) and parallel `scores`.

Example, end to end:

```bash
B=http://127.0.0.1:8765
curl -s -X POST $B/api/session >/dev/null
curl -s -X POST "$B/api/files?name=tutorial_multicopy.ph" \
     --data-binary @examples/tutorial_multicopy.ph >/dev/null
curl -s -X POST $B/api/load -d '{"file":"tutorial_multicopy.ph"}'
curl -s -X POST $B/api/run  -d '{"command":"set seed=42"}'      >/dev/null
curl -s -X POST $B/api/run  -d '{"command":"set criterion=recon"}' >/dev/null
curl -s -X POST $B/api/run  -d '{"command":"hs nreps=10 nthreads=4"}' | python3 -m json.tool
```

All paths in commands are confined to the session's private sandbox directory,
so `savetrees=/etc/x` is relocated to the sandbox and `exe /etc/passwd` is
rejected.

---

## 6. Running the tests

```bash
PYTHONPATH="$PWD" arch -x86_64 /usr/bin/python3 -m pytest clann_web/tests/ -q
```

(or run any test file directly with the same interpreter). They start the real
server on an ephemeral loopback port and drive it end to end.
