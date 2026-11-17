"""python -m tlfit  → opens the fitting interface in the browser."""
import argparse
import os

from .server import serve


def main():
    ap = argparse.ArgumentParser(prog="tlfit", description="Fitting of thermal lens and thermal mirror transients")
    ap.add_argument("--host", default=os.environ.get("TLFIT_HOST", "127.0.0.1"),
                    help="listening address (0.0.0.0 = accept connections from other machines)")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", 8765)), help="local port (default 8765; 0 = any)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the browser automatically")
    a = ap.parse_args()
    serve(a.port, not a.no_browser, a.host)


if __name__ == "__main__":
    main()
