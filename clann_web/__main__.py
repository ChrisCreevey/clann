"""Run the Clann web server: ``clann-web`` (or ``python -m clann_web``).

Starts the local server, prints the URL, and opens it in a browser. Options:
``--host`` / ``--port`` to change the bind address, ``--no-browser`` to skip the
auto-open.

On Apple Silicon the Homebrew-gcc ``libclann-server.so`` is x86_64, so run under
a matching interpreter — either install into an x86_64 venv
(``arch -x86_64 /usr/bin/python3 -m venv venv``) so the ``clann-web`` script is
x86_64, or launch with ``arch -x86_64 /usr/bin/python3 -m clann_web``.
"""

import argparse
import threading
import webbrowser

from .server import serve


def _open_browser_later(url: str, delay: float = 1.0) -> None:
    threading.Timer(delay, lambda: webbrowser.open(url)).start()


def main() -> None:
    ap = argparse.ArgumentParser(prog="clann-web")
    ap.add_argument("--host", default="127.0.0.1",
                    help="bind address (default 127.0.0.1; loopback only)")
    ap.add_argument("--port", type=int, default=8765,
                    help="port to serve on (default 8765)")
    ap.add_argument("--no-browser", action="store_true",
                    help="do not open a browser automatically")
    args = ap.parse_args()

    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"WARNING: binding to {args.host} exposes the engine beyond "
              f"loopback; there is no authentication.")

    if not args.no_browser:
        # localhost is reachable even when bound to 0.0.0.0.
        display_host = "localhost" if args.host in ("0.0.0.0", "::") else args.host
        _open_browser_later(f"http://{display_host}:{args.port}")

    serve(args.host, args.port)


if __name__ == "__main__":
    main()
