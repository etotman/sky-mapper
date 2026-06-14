"""Local web server: static map files + /api endpoints, and server lifecycle."""
import os, sys, json, time, socket, functools, subprocess, webbrowser, threading
import urllib.request
import xml.etree.ElementTree as ET
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, unquote, quote
from astropy.io import fits
from astropy.coordinates import SkyCoord, get_constellation
import astropy.units as u
from config import *
from skymap_io import read_xisf_header
from skymap_simbad import simbad_lookup
from skymap_scan import RESCAN_LOCK, rebuild_catalog

# =============================================================================
# WEB SERVER
# =============================================================================

def read_fits_header_text(filepath: str) -> str:
    out = []
    with fits.open(filepath, ignore_missing_end=True) as hdul:
        for i, hdu in enumerate(hdul):
            out.append(f"===== HDU {i} =====")
            out.append(repr(hdu.header).strip())
    return "\n".join(out)


def read_xisf_header_text(filepath: str) -> str:
    """Render an XISF header to look like a FITS header: the embedded FITS
    keywords as classic `KEY = value / comment` cards, plus the XISF properties
    and image geometry in matching aligned sections."""
    raw = read_xisf_header(filepath)
    try:
        root = ET.fromstring(raw)
    except Exception:
        return raw   # malformed/truncated XML — show as-is

    def local(tag):
        return tag.rsplit('}', 1)[-1]

    fits_cards, props, images = [], [], []
    for el in root.iter():
        tag = local(el.tag)
        if tag == 'FITSKeyword':
            name = (el.get('name') or '').strip()
            val = (el.get('value') or '').strip()
            com = (el.get('comment') or '').strip()
            card = f"{name[:8]:<8}= {val}"
            if com:
                card += f" / {com}"
            fits_cards.append(card)
        elif tag == 'Property':
            pid = (el.get('id') or '').strip()
            val = el.get('value')
            if val is None:
                val = (el.text or '').strip()
            if pid:
                props.append(f"{pid:<34}= {(val or '').strip()}")
        elif tag == 'Image':
            for k in ('geometry', 'sampleFormat', 'colorSpace', 'imageType', 'bounds'):
                if el.get(k):
                    images.append(f"{k:<34}= {el.get(k)}")

    out = []
    if images:
        out += ["===== Image =====", *images, ""]
    out += ["===== FITS Keywords =====", *(fits_cards or ["(none)"])]
    if props:
        out += ["", "===== XISF Properties =====", *props]
    return "\n".join(out)


def _path_is_allowed(real: str) -> bool:
    if not real.lower().endswith(('.fits', '.fit', '.xisf')) or not os.path.isfile(real):
        return False
    return any(real.startswith(os.path.realpath(sd)) for sd, _, _ in SEARCH_DIRS)


class MapRequestHandler(SimpleHTTPRequestHandler):
    """Serves the map folder, plus an /api/header endpoint that reads a FITS/XISF
    header live (so full headers don't have to be embedded in the page)."""

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == '/api/ping':
            self._send_text(f'sky_mapper:{SERVER_VERSION}')
            return
        if parsed.path == '/api/header':
            self._serve_header(parse_qs(parsed.query))
            return
        if parsed.path == '/api/rescan':
            self._serve_rescan(parse_qs(parsed.query))
            return
        if parsed.path == '/api/objectinfo':
            self._serve_objectinfo(parse_qs(parsed.query))
            return
        if parsed.path == '/api/shutdown':
            # Lets a new run stop an older instance cross-platform (no kill/lsof).
            self._send_text('bye')
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        return super().do_GET()

    def _serve_objectinfo(self, query):
        name = unquote((query.get('name') or [''])[0]).strip()
        result = {'ok': True, 'name': name, 'resolved': False}

        def fnum(key):
            try:
                return float((query.get(key) or [''])[0])
            except ValueError:
                return None
        ra, dec = fnum('ra'), fnum('dec')
        if ra is not None and dec is not None:
            try:
                result['constellation'] = get_constellation(SkyCoord(ra, dec, unit='deg'))
            except Exception:
                pass

        if name and name.lower() != 'unknown target':
            try:
                with _OBJINFO_LOCK:
                    sim = simbad_lookup(name)
            except Exception:
                sim = None
            if sim:
                result.update(sim)
                result['resolved'] = True
            result['simbad_url'] = 'https://simbad.cds.unistra.fr/simbad/sim-id?Ident=' + quote(name)
            result['ned_url'] = 'https://ned.ipac.caltech.edu/byname?objname=' + quote(name)
        elif ra is not None and dec is not None:
            # Sexagesimal h/m/s + d/m/s so NED (and SIMBAD) read RA as hours, not
            # degrees — a bare decimal in NED's RA field is interpreted as hours.
            c = SkyCoord(ra, dec, unit='deg')
            ra_hms = c.ra.to_string(unit=u.hour, sep=('h ', 'm ', 's'), precision=0, pad=True)
            dec_dms = c.dec.to_string(unit=u.deg, sep=('d ', 'm ', 's'),
                                      precision=0, alwayssign=True, pad=True)
            result['simbad_url'] = ('https://simbad.cds.unistra.fr/simbad/sim-coo?Coord='
                                    + quote(ra_hms + ' ' + dec_dms) + '&Radius=5&Radius.unit=arcmin')
            result['ned_url'] = ('https://ned.ipac.caltech.edu/conesearch?in_csys=Equatorial'
                                 '&in_equinox=J2000&ra=' + quote(ra_hms)
                                 + '&dec=' + quote(dec_dms) + '&radius=5')
        self._send_json(result)

    def _serve_rescan(self, query):
        rig = unquote((query.get('rig') or ['all'])[0])
        valid = {lbl for _, _, lbl in SEARCH_DIRS}
        only = None if rig == 'all' else {rig}
        if only is not None and not (only & valid):
            self.send_error(400, 'Unknown rig')
            return
        try:
            with RESCAN_LOCK:                       # one rescan at a time
                summary = rebuild_catalog(only_rigs=only)
            self._send_json({'ok': True, **summary})
        except Exception as e:
            self._send_json({'ok': False, 'error': str(e)}, status=500)

    def _send_json(self, obj, status=200):
        body = json.dumps(obj).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _serve_header(self, query):
        path = unquote((query.get('path') or [''])[0])
        real = os.path.realpath(path)
        if not _path_is_allowed(real):
            self.send_error(403, 'Forbidden')
            return
        try:
            if real.lower().endswith('.xisf'):
                text = read_xisf_header_text(real)
            else:
                text = read_fits_header_text(real)
        except Exception as e:
            text = f"Error reading header:\n{e}"
        self._send_text(text)

    def _send_text(self, text):
        body = text.encode('utf-8', errors='replace')
        self.send_response(200)
        self.send_header('Content-Type', 'text/plain; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass   # quiet


def run_server(port: int) -> None:
    handler = functools.partial(MapRequestHandler, directory=SCRIPT_DIR)
    ThreadingHTTPServer(('127.0.0.1', port), handler).serve_forever()


def _server_is_ours(port: int) -> bool:
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/ping', timeout=0.5) as r:
            return r.read().strip() == f'sky_mapper:{SERVER_VERSION}'.encode()
    except Exception:
        return False


def _port_in_use(port: int) -> bool:
    try:
        socket.create_connection(('127.0.0.1', port), timeout=0.3).close()
        return True
    except OSError:
        return False


def _force_free_port(port: int) -> None:
    """Best-effort, cross-platform: free a port held by a process that didn't
    respond to /api/shutdown (e.g. a foreign program or a very old instance)."""
    try:
        if os.name == 'nt':
            out = subprocess.run(['netstat', '-ano'], capture_output=True, text=True).stdout
            pids = {line.split()[-1] for line in out.splitlines()
                    if f':{port}' in line and 'LISTENING' in line.upper()}
            for pid in pids:
                subprocess.run(['taskkill', '/F', '/PID', pid], capture_output=True)
        else:
            subprocess.run(f"kill $(lsof -ti :{port}) 2>/dev/null", shell=True, capture_output=True)
    except Exception:
        pass


def _spawn_server() -> None:
    """Start the background server, detached so it survives this process exiting."""
    kwargs = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if os.name == 'nt':
        # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
        kwargs['creationflags'] = 0x00000008 | 0x00000200
    else:
        kwargs['start_new_session'] = True
    subprocess.Popen(
        [sys.executable, os.path.join(SCRIPT_DIR, "sky-mapper.py"), '--serve', str(WEB_SERVER_PORT)],
        **kwargs,
    )


def ensure_web_server() -> None:
    url = f"http://localhost:{WEB_SERVER_PORT}/{os.path.basename(OUTPUT_HTML)}"
    stop_hint = (f"taskkill /F /PID <pid>" if os.name == 'nt'
                 else f"kill $(lsof -ti :{WEB_SERVER_PORT})")

    if _port_in_use(WEB_SERVER_PORT):
        if _server_is_ours(WEB_SERVER_PORT):
            print(f"\n🌐 Web server already running → {url}  (just refresh your browser)")
            return
        # An older instance (or another program) holds our port. Ask any of our
        # own servers to stop gracefully (cross-platform), then fall back to a
        # best-effort OS kill if something stubborn is still listening.
        print(f"\n♻  Replacing existing server on port {WEB_SERVER_PORT} …")
        try:
            urllib.request.urlopen(f'http://127.0.0.1:{WEB_SERVER_PORT}/api/shutdown', timeout=1).read()
        except Exception:
            pass
        for _ in range(15):
            if not _port_in_use(WEB_SERVER_PORT):
                break
            time.sleep(0.2)
        if _port_in_use(WEB_SERVER_PORT):
            _force_free_port(WEB_SERVER_PORT)
            time.sleep(0.6)
        if _port_in_use(WEB_SERVER_PORT):
            print(f"   ⚠ Port {WEB_SERVER_PORT} is still in use by another program.")
            print(f"     Free it (e.g. {stop_hint}) or set WEB_SERVER_PORT in local_config.py.")
            print(f"     The map file is ready at: {OUTPUT_HTML}")
            return

    _spawn_server()
    print(f"\n🌐 Started web server → {url}")
    print(f"   (Keeps running in the background; stop it with:  {stop_hint})")
    webbrowser.open(url)

