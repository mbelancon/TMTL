"""TMTL2 — aplicativo de desktop.

Ponto de entrada do executável (PyInstaller). Sobe o servidor local numa porta
livre e abre a interface:

* numa janela própria, se o pywebview estiver disponível (Edge WebView2 no
  Windows 10/11);
* senão, no navegador padrão. Nesse caso o programa encerra sozinho cerca de
  30 s depois que a aba é fechada.

Também funciona sem empacotar:  python app.py   [--browser]
"""
import multiprocessing
import os
import sys
import threading
import time
import webbrowser

# poucos núcleos para o BLAS: o paralelismo do ajuste é feito pelo enxame
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from tlfit import server  # noqa: E402

TITLE = "TMTL2 — Ajuste de TL e TM"
IDLE_EXIT = 30  # s sem ping da página → encerra (modo navegador)


def log(msg):
    # no executável sem console, sys.stdout é None
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
            except Exception as e:  # sem WebView2, por exemplo
                log(f"Janela própria indisponível ({e}); abrindo no navegador.")

    log(f"{TITLE}\nInterface em {url}\nFeche esta janela (ou a aba do navegador) para encerrar.")
    webbrowser.open(url)
    server.LAST_PING[0] = time.time() + 20  # tempo para o navegador abrir
    try:
        while time.time() - server.LAST_PING[0] < IDLE_EXIT:
            time.sleep(2)
    except KeyboardInterrupt:
        pass
    httpd.shutdown()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
