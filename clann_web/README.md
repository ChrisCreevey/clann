# clann_web — a local browser client for Clann

Run Clann from a web browser: upload source-tree files, pick commands from
menus, set options, run analyses, and explore the results — including the
interactive tree/reconciliation viewer — without touching a terminal.

It is a **local, single-user tool**: the server binds to loopback (`127.0.0.1`)
and drives Clann through a per-session worker process. See
[`../PLAN_web_client.md`](../PLAN_web_client.md) for the design and roadmap.

---

## 1. Build the server library

The web server talks to a hardened build of Clann — `libclann-server.so` — that
has **no shell/`system()` surface** (safe to sit behind a port). Build it once
from the repository root:

```bash
cd /path/to/clann
make libclann-server.so
```

## 2. Install and run

`clann-web` is a plain-stdlib package (no third-party dependencies). Install it
and launch with a single command:

```bash
# On Apple Silicon: create an x86_64 venv so the interpreter matches the
# x86_64 library (see the note below). On Linux, a normal venv is fine.
arch -x86_64 /usr/bin/python3 -m venv venv
./venv/bin/pip install .

./venv/bin/clann-web            # starts the server and opens your browser
```

You'll see `clann-web serving on http://127.0.0.1:8765` and a browser tab open
there. Stop with `Ctrl+C`. Options:

```bash
./venv/bin/clann-web --port 9000     # different port
./venv/bin/clann-web --no-browser    # don't auto-open a browser
```

`pip install .` copies the freshly built `libclann-server.so` into the package,
so the installed app is self-contained. To develop against the source tree
without reinstalling, use an editable install (`pip install -e .`) — it finds the
library at the repository root.

Prefer not to install? Run it straight from the checkout:

```bash
PYTHONPATH="$PWD" arch -x86_64 /usr/bin/python3 -m clann_web
```

### Why the x86_64 interpreter on Apple Silicon?

The Homebrew `gcc-15` build produces an **x86_64** `libclann-server.so`, so it
must be loaded by a matching (x86_64) Python — a native-arm64 `python3` fails to
`dlopen` it with an *"incompatible architecture"* error. `/usr/bin/python3` ships
an x86_64 slice that runs under Rosetta; building the venv with
`arch -x86_64 /usr/bin/python3` makes the installed `clann-web` script (and the
worker processes it spawns) x86_64 too. On Linux, or once a native/universal
library is built, use your normal Python. If the library lives somewhere
non-default, point at it with `CLANN_LIB=/path/to/libclann-server.so`.

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

## 4. Notes & limitations

- **Long searches stream live** — `POST /api/run` returns a `job_id` immediately;
  the page streams the log over Server-Sent Events and shows a **Stop** button.
- **Cancelling resets the session** — Clann's own interrupt handlers are
  interactive, so *Stop* recovers a runaway/wedged search by killing the session's
  worker process and starting a fresh one. Uploaded files survive; loaded trees
  and `set` options are cleared, so reload the file to continue.
- **One analysis at a time per session** — each session runs its commands in
  order (single engine). Different sessions have independent worker processes and
  can run concurrently.
- **`hs seed=` is silently ignored** — only the `set seed=` command seeds the RNG.
  For reproducible searches, run `set seed=42` *before* `hs` (not `hs seed=42`).
- **Loopback only** — there is no authentication. Do not expose the server beyond
  `127.0.0.1`. Binding another host prints a warning. Uploads are capped (64 MiB)
  and every file path in a command is confined to the session's sandbox.

---

## 5. HTTP API (for scripting)

| Method & path            | Body / query                 | Returns |
|--------------------------|------------------------------|---------|
| `POST /api/session`      | —                            | `{session_id, state, files}` (resets to a clean session) |
| `GET  /api/session`      | —                            | `{session_id, state, files}` |
| `POST /api/files?name=F` | raw file bytes               | `{ok, stored, files}` (name must be a bare filename) |
| `POST /api/load`         | `{"file":"F"}`               | `{ok, log, state}` (F must be uploaded first) |
| `POST /api/run`          | `{"command":"hs nreps=10"}`  | `202 {ok, job_id, command, status}` (async; `409` if a job is already running) |
| `GET  /api/jobs/<id>`    | —                            | `{status, log, trees?, scores?, result_type?, has_viewer?, error?}` |
| `GET  /api/jobs/<id>/stream` | —                        | Server-Sent Events: live `log` chunks, then an `event: done` with the result |
| `POST /api/jobs/<id>/cancel` | —                        | `202 {ok, job_id, status}` (kills + respawns the session worker) |
| `GET  /api/commands`     | —                            | `{commands:[{name, blurb}]}` |
| `GET  /api/commands/<name>/schema` | —                  | typed option schema for the command form |
| `GET  /api/viewer`       | —                            | interactive viewer HTML for the latest tree result |
| `GET  /`                 | —                            | the single-page app |

`status` is one of `running`, `done`, `error`, `cancelled`.

`state` is `{input_file, num_taxa, num_source_trees, criterion, trees_in_memory}`.
For tree-producing commands (`hs`, `nj`, `showtrees`, `reconstruct`, `alltrees`),
`/api/run` also returns `trees` (each `{name, newick, tree, score?, dups?,
losses?}`) and parallel `scores`.

Example, end to end:

Runs are asynchronous — `POST /api/run` returns a `job_id`; poll
`GET /api/jobs/<id>` (or stream `…/stream`) until `status` leaves `running`:

```bash
B=http://127.0.0.1:8765
run() {   # run a command and wait for the job to finish, printing the result
  id=$(curl -s -X POST $B/api/run -d "{\"command\":\"$1\"}" \
       | python3 -c 'import sys,json;print(json.load(sys.stdin)["job_id"])')
  while [ "$(curl -s $B/api/jobs/$id | python3 -c 'import sys,json;print(json.load(sys.stdin)["status"])')" = running ]; do sleep 0.2; done
  curl -s $B/api/jobs/$id | python3 -m json.tool
}
curl -s -X POST $B/api/session >/dev/null
curl -s -X POST "$B/api/files?name=tutorial_multicopy.ph" \
     --data-binary @examples/tutorial_multicopy.ph >/dev/null
curl -s -X POST $B/api/load -d '{"file":"tutorial_multicopy.ph"}'
run "set seed=42"       >/dev/null
run "set criterion=recon" >/dev/null
run "hs nreps=10 nthreads=4"
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
