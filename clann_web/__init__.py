"""clann_web — a local browser client for Clann backed by libclann-server.so.

See PLAN_web_client.md. Step 1.1: engine wrapper + minimal HTTP server.
"""
from .engine import ClannEngine, ClannError
from .server import make_server, serve

__all__ = ["ClannEngine", "ClannError", "make_server", "serve"]
