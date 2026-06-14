# Sky Mapper — scans FITS/XISF files, extracts RA/DEC, target names, exposure and
# capture metadata (plate-solved WCS where available), clusters nearby frames, and
# generates an interactive Aladin Lite sky map. A small local web server serves the
# map plus an API for live FITS/XISF headers, SIMBAD lookups, and rescans.
# Calibration frames, mosaics, and coordinate-less test frames are handled
# automatically; UI state persists in the browser.
#
# This entry point is intentionally tiny; the implementation lives in the
# skymap_* modules and config.py, with the page markup in template.html.
import sys

# Make console output UTF-8 safe everywhere. On Windows, when stdout is a pipe or
# redirected to a file, Python defaults to the legacy locale codepage (e.g. cp1252)
# and the emoji/arrow status lines would raise UnicodeEncodeError mid-scan.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

from config import WEB_SERVER_PORT
from skymap_scan import rebuild_catalog
from skymap_server import run_server, ensure_web_server


def main():
    # Background server mode (spawned by ensure_web_server): just serve and block.
    if '--serve' in sys.argv:
        idx = sys.argv.index('--serve')
        port = int(sys.argv[idx + 1]) if len(sys.argv) > idx + 1 else WEB_SERVER_PORT
        run_server(port)
        return

    rebuild_catalog()
    ensure_web_server()


if __name__ == "__main__":
    main()
