"""python -m tlfit  → abre a interface de ajuste no navegador."""
import argparse
import os

from .server import serve


def main():
    ap = argparse.ArgumentParser(prog="tlfit", description="Ajuste de transientes de lente térmica")
    ap.add_argument("--host", default=os.environ.get("TLFIT_HOST", "127.0.0.1"),
                    help="endereço de escuta (0.0.0.0 = aceitar conexões de outras máquinas)")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", 8765)), help="porta local (padrão 8765; 0 = qualquer)")
    ap.add_argument("--no-browser", action="store_true", help="não abrir o navegador automaticamente")
    a = ap.parse_args()
    serve(a.port, not a.no_browser, a.host)


if __name__ == "__main__":
    main()
