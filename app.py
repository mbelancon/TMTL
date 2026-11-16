"""TMTL2 — desktop application.

Entry point of the application. Starts the local server on a free port and
opens the interface:

* in its own window, if pywebview is available (Edge WebView2 on
  Windows 10/11);
* otherwise, in the default browser. In that case the program exits by itself
  about 30 s after the tab is closed.

Run it with:  python app.py   [--browser]
"""
import multiprocessing
import os
import sys
import threading
import time
import webbrowser

# few cores for BLAS: the fit is parallelized by the swarm/sweep threads
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from tlfit import server  # noqa: E402

TITLE = "TMTL2 — TL and TM fitting"
IDLE_EXIT = 30  # s without a ping from the page → exit (browser mode)


def log(msg):
    # in an executable without a console, sys.stdout is None
    if sys.stdout is not None:
        print(msg, flush=True)


def main():
    httpd = server.make_server(0)
    url = f"http://127.0.0.1:{httpd.server_address[1]}/"
    threading.Thread(target=httpd.serve_forever, daemon=True).start()

    if "--browser" not in sys.argv:
        try:
            import webview  # pywebview
        except ImportError:
            webview = None
        if webview is not None:
            try:
                webview.create_window(TITLE, url, width=1440, height=900, min_size=(900, 600))
                webview.start()
                httpd.shutdown()
                return
            except Exception as e:  # no WebView2, for example
                log(f"Own window unavailable ({e}); opening in the browser.")

    log(f"{TITLE}\nInterface at {url}\nClose this window (or the browser tab) to quit.")
    webbrowser.open(url)
    server.LAST_PING[0] = time.time() + 20  # time for the browser to open
    try:
        while time.time() - server.LAST_PING[0] < IDLE_EXIT:
            time.sleep(2)
    except KeyboardInterrupt:
        pass
    httpd.shutdown()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
