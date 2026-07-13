"""Run the Clann web server: `python -m clann_web [--host H] [--port P]`.

On Apple Silicon, libclann-server.so is x86_64, so launch under a matching
interpreter, e.g.:  arch -x86_64 /usr/bin/python3 -m clann_web
"""

import argparse

from .server import serve


def main() -> None:
    ap = argparse.ArgumentParser(prog="clann_web")
    ap.add_argument("--host", default="127.0.0.1",
                    help="bind address (default 127.0.0.1; loopback only)")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"WARNING: binding to {args.host} exposes the engine beyond "
              f"loopback; there is no authentication.")
    serve(args.host, args.port)


if __name__ == "__main__":
    main()
