import os
import sys
import json
import glob
import re
import math
import time
import socket
import threading
import functools
import subprocess
import webbrowser
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, unquote, quote
import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, get_constellation
import astropy.units as u

# sky_mapper11 — scans FITS/XISF files, extracts RA/DEC + target names + exposure
# metadata, clusters nearby targets, and generates an interactive Aladin Lite map.
#
# New in v11:
#   • Cluster data is embedded as JSON; the page builds overlays dynamically.
#   • Click a footprint → details panel (frames, integration time per filter, dates).
#   • Per-rig show/hide checkboxes in the legend.
#   • "Dim targets not imageable tonight" mode (uses OBSERVER_LAT/LON below).
#   • Export target list as CSV.
#   • Labels auto-hide when zoomed out past LABEL_HIDE_FOV_DEG.
#   • Optional constellation lines overlay (checkbox; data auto-downloaded once).
#   • Calibration frames (dark/flat/bias) detected via IMAGETYP header (or filename
#     tokens when no header) and skipped.
#   • Files without coordinates (test images etc.) are skipped and remembered.
#   • Mosaic frames (filename/target contains "mosaic"/"panel") can be hidden
#     with a checkbox.
#   • Local web server starts automatically (only if not already running).
#   • UI state (survey, checkboxes, rig toggles) persists in the browser.
#   • Cluster centers averaged with unit vectors (fixes RA 0/360 wraparound).
#   • XISF header length read from the file header instead of a fixed 16 KB.
#
# All output files (cache, HTML, constellation data) live next to this script,
# regardless of the directory you run it from.

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# One entry per equipment series / folder: (absolute_path, hex_color, label).
# These are placeholder defaults — set your real image folders in local_config.py
# (git-ignored). See local_config.example.py.
SEARCH_DIRS = [
    (r"/path/to/rig_a/images/", "#00ff00", "Rig A"),   # green
    (r"/path/to/rig_b/images/", "#ff4444", "Rig B"),   # red
    # Add more rigs here:
    # (r"/Volumes/NAS/rig_c/", "#44aaff", "Rig C"),
]

# Target keywords to exclude from map rendering (case-insensitive)
EXCLUDE_KEYWORDS = [
    "flatwizard",
    "snapshot",
    "tsuchinshan",
]

# Filenames made only of these tokens are treated as calibration frames when the
# file has no IMAGETYP header. (Header check happens first, so a *target* named
# "Dark Shark" shot with NINA/SGP is still kept — its IMAGETYP says LIGHT.)
CALIBRATION_TOKENS = {
    "dark", "darks", "flat", "flats", "bias", "biases",
    "darkflat", "darkflats", "masterdark", "masterflat", "masterbias",
}

# A frame counts as part of a mosaic if its filename or target contains one of these.
MOSAIC_KEYWORDS = ["mosaic", "panel"]

# A target is also treated as a mosaic when its name appears in this many (or more)
# separate pointings — same object, multiple panels. Note: a moving target (comet)
# imaged on several nights also trips this.
MOSAIC_MIN_PANELS = 3

CACHE_FILE  = os.path.join(SCRIPT_DIR, "astrophoto_cache_v2.json")
OUTPUT_HTML = os.path.join(SCRIPT_DIR, "aladin_map.html")
CACHE_SCHEMA = 7    # v7: also capture conditions (altitude, airmass, temps)

GROUPING_TOLERANCE_DEG = 0.5

# Fallback footprint size, used only when a file's headers lack the data needed
# to compute the true field of view (NAXIS / XPIXSZ / FOCALLEN).
FRAME_WIDTH_DEG  = 1.0
FRAME_HEIGHT_DEG = 1.0

# Manual footprint nudge. Normally leave at 0: plate-solved WCS centers (Rig A)
# and post-centering mount coordinates (Rig B) are accurate. The old -0.22/+0.10
# values were compensating for mount-reported positions, which are no longer used
# when a solved WCS is present.
RA_OFFSET_DEG  = 0.0
DEC_OFFSET_DEG = 0.0

# Observing site for "Dim targets too low tonight". Leave as None to auto-detect
# from the most common SITELAT/SITELONG in your image headers; set explicitly to
# override (degrees; longitude EAST positive, so US longitudes are negative).
OBSERVER_LAT     = None
OBSERVER_LON     = None
MIN_ALTITUDE_DEG = 30.0    # a target must reach this altitude during darkness

LABEL_HIDE_FOV_DEG = 90.0  # hide labels when zoomed out wider than this FOV

# Default control-panel position in pixels from the top-left of the window.
# The panel can also be dragged by its title bar; a dragged position is
# remembered by the browser and overrides these defaults.
PANEL_X = 10
PANEL_Y = 90

# Folder-browser panel default position (also draggable; remembered per browser).
# Defaults to the right of the main control panel so they don't overlap.
FOLDER_PANEL_X = 360
FOLDER_PANEL_Y = 90

# Default outline color for footprints highlighted via the folder browser.
# Customize here, or change it live with the color swatch in the folder panel.
HIGHLIGHT_COLOR = "#ff00ff"   # magenta

WEB_SERVER_PORT = 8001         # v12 uses its own port so it can run alongside v11
SERVER_VERSION  = 12           # bump when the server's API code changes, to force a
                               # running background server to be replaced on next run

CONSTELLATION_FILE = os.path.join(SCRIPT_DIR, "constellations.lines.json")
CONSTELLATION_URL  = "https://cdn.jsdelivr.net/gh/ofrohn/d3-celestial@master/data/constellations.lines.json"

# Cache of SIMBAD reference data per object name (so the details-panel info
# section doesn't re-query SIMBAD on every click).
OBJECT_INFO_CACHE = os.path.join(SCRIPT_DIR, "object_info_cache.json")

# Messier/Caldwell catalog (generated once by generate_catalog.py) for the
# completion overlay. The page fetches it when the Catalog panel is opened.
CATALOG_FILE = os.path.join(SCRIPT_DIR, "messier_caldwell.json")

# ---------------------------------------------------------------------------
# Local overrides. Put your real rig folders, observing site, and any preference
# tweaks in local_config.py (git-ignored) — it may redefine any constant above.
# Copy local_config.example.py to local_config.py to get started.
# ---------------------------------------------------------------------------
try:
    from local_config import *  # noqa: F401,F403
except ImportError:
    pass


# =============================================================================
# FILENAME-BASED TARGET EXTRACTION
# =============================================================================

_CATALOGUE_RE = re.compile(
    r'\b('
    r'ic\s*\d+'           # IC objects
    r'|ngc\s*\d+'         # NGC objects
    r'|m\s*\d{1,3}'       # Messier objects  (M1–M110)
    r'|sh2[-_\s]\d+'      # Sharpless
    r'|lbn[-_\s]\d+'      # LBN
    r'|ldn[-_\s]\d+'      # LDN  (e.g. ldn1366b → LDN 1366)
    r'|vdb\s*\d+'         # vdB
    r'|b\s*\d+'           # Barnard
    r'|abell\s*\d+'       # Abell
    r')',
    re.IGNORECASE,
)

def target_from_filename(filepath: str) -> str | None:
    stem = os.path.splitext(os.path.basename(filepath))[0]
    tokens = stem.split('_')

    candidates = list(tokens)
    for i in range(len(tokens) - 1):
        candidates.append(tokens[i] + tokens[i + 1])

    for candidate in candidates:
        m = _CATALOGUE_RE.match(candidate.strip())
        if m:
            raw = m.group(0)
            name = re.sub(r'([A-Za-z]+)[-_\s]*(\d+)', r'\1 \2', raw).upper()
            return name.strip()

    return None


def should_exclude_target(target_name: str) -> bool:
    if not target_name:
        return False
    target_lower = target_name.lower()
    return any(keyword.lower() in target_lower for keyword in EXCLUDE_KEYWORDS)


def filename_is_calibration(filepath: str) -> bool:
    stem = os.path.splitext(os.path.basename(filepath))[0].lower()
    tokens = re.split(r'[\s_\-.]+', stem)
    return any(t in CALIBRATION_TOKENS for t in tokens)


def is_mosaic(filepath: str, target: str) -> bool:
    text = (os.path.basename(filepath) + ' ' + (target or '')).lower()
    return any(k in text for k in MOSAIC_KEYWORDS)


def _is_calibration_imagetyp(imagetyp: str) -> bool:
    return any(k in imagetyp for k in ('dark', 'flat', 'bias'))


# =============================================================================
# FILE READERS
# =============================================================================

def parse_radec(ra_val, dec_val):
    """Decimal degrees if both values are numeric, otherwise sexagesimal (RA in hours)."""
    try:
        return float(ra_val), float(dec_val)
    except (TypeError, ValueError):
        coord = SkyCoord(f"{ra_val} {dec_val}", unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg


def _first(header, *keys):
    for key in keys:
        val = header.get(key)
        if val is not None and str(val).strip() != '':
            return val
    return None


def _fov_deg(naxis, pixsz_um, focallen_mm):
    """Field of view of one axis: pixels × pixel size (µm) / focal length (mm)."""
    try:
        naxis, pixsz_um, focallen_mm = float(naxis), float(pixsz_um), float(focallen_mm)
        if naxis > 0 and pixsz_um > 0 and focallen_mm > 0:
            return round(math.degrees(naxis * pixsz_um / 1000.0 / focallen_mm), 4)
    except (TypeError, ValueError):
        pass
    return None


def _num_kw(getter, *keys):
    """First parseable float among keys, using a getter(key)->value callable."""
    for k in keys:
        v = getter(k)
        if v is not None and str(v).strip() not in ('', '~'):
            try:
                return round(float(v), 3)
            except (TypeError, ValueError):
                pass
    return None


def _wcs_center_and_rotation(header):
    """Field-center sky position and camera position angle from a solved WCS.

    Evaluated at the reference pixel CRPIX (where the solution is exact and SIP
    distortion is zero), NOT at NAXIS/2. Some ASIAIR Live/Preview stacks solve a
    downsampled image, so CRPIX sits near the downsampled-grid center — far from
    the full-frame NAXIS/2. Evaluating at NAXIS/2 there extrapolates the SIP
    polynomial and lands the footprint degrees away (the M51 bug); CRPIX yields
    the true field center for normal and downsampled solves alike."""
    cx, cy = header.get('CRPIX1'), header.get('CRPIX2')
    nx, ny = header.get('NAXIS1'), header.get('NAXIS2')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wcs = WCS(header).celestial
        if not wcs.has_celestial:
            return None
        if cx is not None and cy is not None:
            px, py = cx - 1.0, cy - 1.0      # FITS CRPIX is 1-based; astropy is 0-based
        elif nx and ny:
            px, py = nx / 2.0, ny / 2.0
        else:
            return None
        center = wcs.pixel_to_world(px, py)
        up = wcs.pixel_to_world(px, py + 10.0)
        return (float(center.ra.deg), float(center.dec.deg),
                float(center.position_angle(up).deg))


def get_fits_data(filepath: str) -> dict:
    try:
        with fits.open(filepath, ignore_missing_end=True) as hdul:
            header = hdul[0].header

            imagetyp = str(_first(header, 'IMAGETYP', 'FRAME') or '').strip().lower()
            if imagetyp:
                if _is_calibration_imagetyp(imagetyp):
                    return {'skip': 'calibration'}
            elif filename_is_calibration(filepath):
                return {'skip': 'calibration'}

            ra_val  = _first(header, 'RA', 'OBJCTRA')
            dec_val = _first(header, 'DEC', 'OBJCTDEC')
            solved = None
            if header.get('CRVAL1') is not None and header.get('CRVAL2') is not None:
                try:
                    solved = _wcs_center_and_rotation(header)
                except Exception:
                    solved = None
            if solved is None and (ra_val is None or dec_val is None):
                return {'skip': 'no_coords'}

            rot = 0.0
            if solved:
                # plate-solved WCS beats the mount-reported position
                ra, dec, rot = solved
            else:
                ra, dec = parse_radec(ra_val, dec_val)
                try:
                    rot = float(_first(header, 'OBJCTROT', 'CROTA2', 'CROTA1') or 0)
                except (TypeError, ValueError):
                    rot = 0.0

            fov_w = _fov_deg(header.get('NAXIS1'),
                             header.get('XPIXSZ'), header.get('FOCALLEN'))
            fov_h = _fov_deg(header.get('NAXIS2'),
                             header.get('YPIXSZ') or header.get('XPIXSZ'),
                             header.get('FOCALLEN'))

            target = str(_first(header, 'OBJECT', 'OBJNAME') or '').strip()
            if not target or target.lower() == 'unknown':
                target = target_from_filename(filepath) or "Unknown Target"

            try:
                exptime = float(_first(header, 'EXPTIME', 'EXPOSURE') or 0)
            except (TypeError, ValueError):
                exptime = 0.0
            filt     = str(header.get('FILTER') or '').strip() or None
            date_obs = str(header.get('DATE-OBS') or '').strip()[:10] or None

            site = None
            try:
                site_lat, site_lon = header.get('SITELAT'), header.get('SITELONG')
                if site_lat is not None and site_lon is not None:
                    site = [round(float(site_lat), 4), round(float(site_lon), 4)]
            except (TypeError, ValueError):
                pass

            return {'ra': ra, 'dec': dec, 'target': target, 'exptime': exptime,
                    'filter': filt, 'date_obs': date_obs, 'site': site,
                    'fov_w': fov_w, 'fov_h': fov_h, 'rot': round(rot, 2),
                    'alt': _num_kw(header.get, 'CENTALT', 'OBJCTALT', 'ALTITUDE'),
                    'airmass': _num_kw(header.get, 'AIRMASS'),
                    'ccdtemp': _num_kw(header.get, 'CCD-TEMP', 'CCDTEMP'),
                    'foctemp': _num_kw(header.get, 'FOCTEMP', 'FOCUSTEM', 'AMBTEMP')}
    except Exception as e:
        print(f"  ⚠ Error reading FITS {filepath}: {e}")
        return {'skip': 'error'}


def read_xisf_header(filepath: str) -> str:
    """Read the XML header block, sized from the XISF signature (fallback: 64 KB)."""
    with open(filepath, 'rb') as f:
        sig = f.read(8)
        if sig == b'XISF0100':
            header_len = int.from_bytes(f.read(4), 'little')
            f.read(4)  # reserved
            return f.read(min(header_len, 8 * 1024 * 1024)).decode('utf-8', errors='ignore')
        f.seek(0)
        return f.read(65536).decode('ascii', errors='ignore')


def _xisf_property(header_xml: str, prop_id: str) -> str | None:
    m = re.search(r'<Property id="' + re.escape(prop_id) + r'"[^>]*value="([^"]*)"', header_xml)
    if m:
        return m.group(1).strip()
    m = re.search(r'<Property id="' + re.escape(prop_id) + r'"[^>]*>([^<]+)</Property>', header_xml)
    return m.group(1).strip() if m else None


def _xisf_fitskw(header_xml: str, name: str) -> str | None:
    m = re.search(r'<FITSKeyword[^>]*name="' + re.escape(name) + r'"[^>]*value="\'?\s*([^"\']*?)\s*\'?"',
                  header_xml)
    return m.group(1).strip() if m else None


def _xisf_solved_center(header_xml: str):
    """If the XISF embeds a plate-solved WCS (CRVAL present as FITS keywords),
    return (ra, dec, rotation) computed at CRPIX. Processed XISF exports often
    store a wrong `Observation:Center:RA` (computed at NAXIS/2 over a downsampled
    solve — the M51 bug), but a correct CRVAL, so prefer the WCS when present."""
    if _xisf_fitskw(header_xml, 'CRVAL1') is None or _xisf_fitskw(header_xml, 'CRVAL2') is None:
        return None
    h = fits.Header()
    for k in ('CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
              'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2'):
        v = _xisf_fitskw(header_xml, k)
        if v is not None:
            try:
                h[k] = float(v)
            except ValueError:
                pass
    # Force a plain TAN projection: we only need center + linear rotation, and any
    # SIP coefficients aren't reconstructed from the XISF keyword list.
    h['CTYPE1'], h['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    try:
        return _wcs_center_and_rotation(h)
    except Exception:
        return None


def get_xisf_data(filepath: str) -> dict:
    try:
        hdr = read_xisf_header(filepath)

        imagetyp = (_xisf_fitskw(hdr, 'IMAGETYP')
                    or _xisf_property(hdr, 'Observation:ImageType') or '').lower()
        if imagetyp:
            if _is_calibration_imagetyp(imagetyp):
                return {'skip': 'calibration'}
        elif filename_is_calibration(filepath):
            return {'skip': 'calibration'}

        ra_val  = (_xisf_property(hdr, 'Observation:Center:RA')
                   or _xisf_fitskw(hdr, 'RA') or _xisf_fitskw(hdr, 'OBJCTRA'))
        dec_val = (_xisf_property(hdr, 'Observation:Center:Dec')
                   or _xisf_fitskw(hdr, 'DEC') or _xisf_fitskw(hdr, 'OBJCTDEC'))

        solved = _xisf_solved_center(hdr)   # embedded WCS wins over header center
        if solved is not None:
            ra, dec, solved_rot = solved
        elif ra_val and dec_val:
            ra, dec = parse_radec(ra_val, dec_val)
            solved_rot = None
        else:
            return {'skip': 'no_coords'}

        target = (_xisf_property(hdr, 'Observation:Object:Name')
                  or _xisf_fitskw(hdr, 'OBJECT') or '').strip()
        if not target or target.lower() == 'unknown':
            target = target_from_filename(filepath) or "Unknown Target"

        try:
            exptime = float(_xisf_property(hdr, 'Instrument:ExposureTime')
                            or _xisf_fitskw(hdr, 'EXPTIME')
                            or _xisf_fitskw(hdr, 'EXPOSURE') or 0)
        except (TypeError, ValueError):
            exptime = 0.0
        filt = (_xisf_property(hdr, 'Instrument:Filter:Name')
                or _xisf_fitskw(hdr, 'FILTER')) or None
        date_obs = (_xisf_property(hdr, 'Observation:Time:Start')
                    or _xisf_fitskw(hdr, 'DATE-OBS') or '')[:10] or None

        site = None
        try:
            site_lat = (_xisf_property(hdr, 'Observation:Location:Latitude')
                        or _xisf_fitskw(hdr, 'SITELAT'))
            site_lon = (_xisf_property(hdr, 'Observation:Location:Longitude')
                        or _xisf_fitskw(hdr, 'SITELONG'))
            if site_lat is not None and site_lon is not None:
                site = [round(float(site_lat), 4), round(float(site_lon), 4)]
        except (TypeError, ValueError):
            pass

        if solved_rot is not None:
            rot = solved_rot
        else:
            try:
                rot = float(_xisf_fitskw(hdr, 'OBJCTROT') or 0)
            except (TypeError, ValueError):
                rot = 0.0

        geom = re.search(r'<Image[^>]*geometry="(\d+):(\d+):', hdr)
        nx, ny = (int(geom.group(1)), int(geom.group(2))) if geom else (None, None)
        pixsz_x = _xisf_fitskw(hdr, 'XPIXSZ') or _xisf_property(hdr, 'Instrument:Sensor:XPixelSize')
        pixsz_y = _xisf_fitskw(hdr, 'YPIXSZ') or _xisf_property(hdr, 'Instrument:Sensor:YPixelSize') or pixsz_x
        focal = _xisf_fitskw(hdr, 'FOCALLEN')
        if not focal:
            focal_m = _xisf_property(hdr, 'Instrument:Telescope:FocalLength')  # metres
            try:
                focal = float(focal_m) * 1000.0 if focal_m else None
            except (TypeError, ValueError):
                focal = None
        fov_w = _fov_deg(nx, pixsz_x, focal)
        fov_h = _fov_deg(ny, pixsz_y, focal)

        kw = lambda k: _xisf_fitskw(hdr, k)
        return {'ra': ra, 'dec': dec, 'target': target, 'exptime': exptime,
                'filter': filt, 'date_obs': date_obs, 'site': site,
                'fov_w': fov_w, 'fov_h': fov_h, 'rot': round(rot, 2),
                'alt': _num_kw(kw, 'CENTALT', 'OBJCTALT', 'ALTITUDE'),
                'airmass': _num_kw(kw, 'AIRMASS'),
                'ccdtemp': _num_kw(kw, 'CCD-TEMP', 'CCDTEMP'),
                'foctemp': _num_kw(kw, 'FOCTEMP', 'FOCUSTEM', 'AMBTEMP')}
    except Exception as e:
        print(f"  ⚠ Error reading XISF {filepath}: {e}")
        return {'skip': 'error'}


# =============================================================================
# CACHE
# =============================================================================

def load_cache() -> dict:
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            return json.load(f)
    return {}


def save_cache(cache: dict) -> None:
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=1)


# =============================================================================
# CLUSTERING  (unit-vector centroids — safe across the RA 0/360 boundary)
# =============================================================================

def _radec_to_vector(ra_deg, dec_deg):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return (math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec))


def _vector_to_radec(v):
    x, y, z = v
    norm = math.sqrt(x * x + y * y + z * z) or 1.0
    ra = math.degrees(math.atan2(y, x)) % 360
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z / norm))))
    return ra, dec


def cluster_by_distance(entries: list, tolerance_deg: float) -> list:
    clusters = []
    cos_tol = math.cos(math.radians(tolerance_deg))

    for entry in entries:
        v = _radec_to_vector(entry['ra'], entry['dec'])

        best_idx, best_dot = None, -2.0
        for i, cl in enumerate(clusters):
            cx, cy, cz = cl['vsum']
            norm = math.sqrt(cx * cx + cy * cy + cz * cz) or 1.0
            dot = (v[0] * cx + v[1] * cy + v[2] * cz) / norm
            if dot > best_dot:
                best_dot, best_idx = dot, i

        if best_idx is not None and best_dot >= cos_tol:
            cl = clusters[best_idx]
            cl['vsum'] = (cl['vsum'][0] + v[0], cl['vsum'][1] + v[1], cl['vsum'][2] + v[2])
            cl['members'].append(entry)
        else:
            clusters.append({'vsum': v, 'members': [entry]})

    return clusters


def _median(values, fallback):
    values = sorted(v for v in values if v)
    return values[len(values) // 2] if values else fallback


def _circular_mean_deg(angles):
    if not angles:
        return 0.0
    s = sum(math.sin(math.radians(a)) for a in angles)
    c = sum(math.cos(math.radians(a)) for a in angles)
    if abs(s) < 1e-9 and abs(c) < 1e-9:
        return 0.0
    return math.degrees(math.atan2(s, c))


def footprint_corners(ra, dec, fov_w, fov_h, rot, ndigits=6):
    """Four sky corners of a frame centered at ra/dec, rotated by position angle
    `rot` (degrees east of north), on a local tangent plane."""
    cos_dec = max(math.cos(math.radians(dec)), 0.01)
    th = math.radians(rot)
    corners = []
    for u, v in ((-fov_w / 2, -fov_h / 2), (fov_w / 2, -fov_h / 2),
                 (fov_w / 2, fov_h / 2), (-fov_w / 2, fov_h / 2)):
        east  = u * math.cos(th) + v * math.sin(th)
        north = -u * math.sin(th) + v * math.cos(th)
        corners.append([round(ra + east / cos_dec, ndigits), round(dec + north, ndigits)])
    return corners


def build_cluster_payload(clusters: list) -> list:
    payload = []
    for i, cl in enumerate(clusters):
        ra, dec = _vector_to_radec(cl['vsum'])
        members = cl['members']

        names = [m.get('target') or 'Unknown Target' for m in members]
        valid = [n for n in names if n != 'Unknown Target']
        name = max(set(valid), key=valid.count) if valid else 'Unknown Target'

        colors = [m['color'] for m in members]
        color = max(set(colors), key=colors.count)

        filters = {}
        total_s = 0.0
        dates = []
        for m in members:
            fname = m.get('filter') or 'None/OSC'
            try:
                sec = float(m.get('exptime') or 0)
            except (TypeError, ValueError):
                sec = 0.0
            slot = filters.setdefault(fname, {'frames': 0, 'seconds': 0.0})
            slot['frames'] += 1
            slot['seconds'] = round(slot['seconds'] + sec, 1)
            total_s += sec
            if m.get('date_obs'):
                dates.append(m['date_obs'])

        ra_c  = ra  + RA_OFFSET_DEG
        dec_c = dec + DEC_OFFSET_DEG

        # True frame size and camera angle from the member files' headers.
        fov_w = _median([m.get('fov_w') for m in members], FRAME_WIDTH_DEG)
        fov_h = _median([m.get('fov_h') for m in members], FRAME_HEIGHT_DEG)
        rot   = _circular_mean_deg([m.get('rot') or 0.0 for m in members])

        corners = footprint_corners(ra_c, dec_c, fov_w, fov_h, rot)

        ra_lo,  ra_hi  = min(p[0] for p in corners), max(p[0] for p in corners)
        dec_lo, dec_hi = min(p[1] for p in corners), max(p[1] for p in corners)

        payload.append({
            'id': i,
            'ra': round(ra_c, 6),
            'dec': round(dec_c, 6),
            'corners': corners,
            'rect': [round(ra_lo, 6), round(ra_hi, 6),
                     round(dec_lo, 6), round(dec_hi, 6)],
            'label_dec': round(dec_hi + 0.18, 6),
            'name': name,
            'label': f"{name} ({len(members)})",
            'alt_names': sorted({n for n in valid if n != name})[:8],
            'count': len(members),
            'color': color,
            'mosaic': any(m.get('mosaic') for m in members),
            'rigs': sorted({m.get('rig') for m in members if m.get('rig')}),
            'filters': filters,
            'total_seconds': round(total_s, 1),
            'date_min': min(dates) if dates else None,
            'date_max': max(dates) if dates else None,
            'alt': _median([m.get('alt') for m in members if m.get('alt') is not None], None),
            'airmass': _median([m.get('airmass') for m in members if m.get('airmass') is not None], None),
        })

    # Same target name spread across several clusters = mosaic panels.
    name_spread = Counter(c['name'] for c in payload if c['name'] != 'Unknown Target')
    for c in payload:
        if name_spread.get(c['name'], 0) >= MOSAIC_MIN_PANELS:
            c['mosaic'] = True

    return payload


def build_folder_data(cache: dict) -> dict:
    """Folder tree (one root per rig) plus the unique footprints in each folder,
    for the click-a-folder-to-highlight browser. Footprints are deduplicated
    globally (coarsely, so dithered subs collapse) and referenced by index."""
    rig_of_dir = {sd: lbl for sd, _, lbl in SEARCH_DIRS}

    fp_index, footprints = {}, []          # dedup key -> idx,  list of corner sets
    node_ids, nodes = {}, []               # path key -> int id, list of node dicts

    def get_node(key, label, parent_id, rig):
        nid = node_ids.get(key)
        if nid is None:
            nid = len(nodes)
            node_ids[key] = nid
            nodes.append({'label': label, 'parent': parent_id,
                          'rig': rig, 'fp': [], 'count': 0, 'total': 0})
        return nid

    for filepath, entry in cache.items():
        if entry.get('skip') or 'ra' not in entry:
            continue
        if should_exclude_target(entry.get('target', '')):
            continue
        sd = entry.get('search_dir', '')
        if not sd:
            continue
        rig = entry.get('rig') or rig_of_dir.get(sd, sd)

        root_key = 'rig:' + (rig or sd)
        leaf_id = get_node(root_key, rig or sd, None, rig)

        try:
            rel = os.path.relpath(filepath, sd)
        except ValueError:
            continue
        folder = os.path.dirname(rel)
        if folder and folder != '.':
            accum = root_key
            parent_id = leaf_id
            for part in folder.split(os.sep):
                accum += '/' + part
                parent_id = get_node(accum, part, parent_id, rig)
            leaf_id = parent_id

        ra, dec = entry['ra'], entry['dec']
        fov_w = entry.get('fov_w') or FRAME_WIDTH_DEG
        fov_h = entry.get('fov_h') or FRAME_HEIGHT_DEG
        rot   = entry.get('rot') or 0.0
        # coarse key (0.01° ≈ 36") so dithered/near-identical subs share a footprint
        key = (round(ra, 2), round(dec, 2), round(fov_w, 2), round(fov_h, 2), round(rot, 0))
        idx = fp_index.get(key)
        if idx is None:
            idx = len(footprints)
            fp_index[key] = idx
            footprints.append(footprint_corners(ra, dec, fov_w, fov_h, rot, ndigits=5))

        node = nodes[leaf_id]
        node['count'] += 1
        if idx not in node['fp']:
            node['fp'].append(idx)
        node.setdefault('files', []).append([os.path.basename(filepath), idx])
        if 'adir' not in node:
            node['adir'] = os.path.dirname(filepath)

    for n in nodes:
        if n.get('files'):
            n['files'].sort(key=lambda f: f[0])

    # recursive file totals (for display on parent folders)
    children = {}
    for nid, n in enumerate(nodes):
        children.setdefault(n['parent'], []).append(nid)

    def compute_total(nid):
        total = nodes[nid]['count']
        for child in children.get(nid, []):
            total += compute_total(child)
        nodes[nid]['total'] = total
        return total

    for nid, n in enumerate(nodes):
        if n['parent'] is None:
            compute_total(nid)

    return {'footprints': footprints, 'tree': nodes}


def build_stats(cache: dict) -> dict:
    """Archive-wide statistics for the dashboard, computed from per-file cache."""
    total_frames = 0
    total_sec = 0.0
    by_filter, by_rig, by_month = {}, {}, {}
    target_sec, target_frames = {}, {}
    nights = set()
    alts, airmasses, foctemps = [], [], []

    for entry in cache.values():
        if entry.get('skip') or 'ra' not in entry:
            continue
        if should_exclude_target(entry.get('target', '')):
            continue
        try:
            sec = float(entry.get('exptime') or 0)
        except (TypeError, ValueError):
            sec = 0.0
        total_frames += 1
        total_sec += sec
        by_filter[entry.get('filter') or 'None/OSC'] = by_filter.get(entry.get('filter') or 'None/OSC', 0.0) + sec
        rig = entry.get('rig') or '?'
        by_rig[rig] = by_rig.get(rig, 0.0) + sec
        d = entry.get('date_obs')
        if d:
            by_month[d[:7]] = by_month.get(d[:7], 0.0) + sec
            nights.add(d)
        t = entry.get('target', 'Unknown Target')
        if t and t != 'Unknown Target':
            target_sec[t] = target_sec.get(t, 0.0) + sec
            target_frames[t] = target_frames.get(t, 0) + 1
        if entry.get('alt') is not None:
            alts.append(entry['alt'])
        if entry.get('airmass') is not None:
            airmasses.append(entry['airmass'])
        if entry.get('foctemp') is not None:
            foctemps.append(entry['foctemp'])

    def med(a):
        return round(sorted(a)[len(a) // 2], 2) if a else None

    top = sorted(target_sec.items(), key=lambda kv: -kv[1])[:15]
    return {
        'total_frames': total_frames,
        'total_seconds': round(total_sec),
        'targets': len(target_sec),
        'nights': len(nights),
        'date_min': min(nights) if nights else None,
        'date_max': max(nights) if nights else None,
        'by_filter': {k: round(v) for k, v in sorted(by_filter.items(), key=lambda kv: -kv[1])},
        'by_rig': {k: round(v) for k, v in sorted(by_rig.items(), key=lambda kv: -kv[1])},
        'by_month': {k: round(v) for k, v in sorted(by_month.items())},
        'top_targets': [{'name': n, 'seconds': round(s), 'frames': target_frames[n]} for n, s in top],
        'median_alt': med(alts),
        'median_airmass': med(airmasses),
        'median_foctemp': med(foctemps),
        'cond_count': len(alts),
    }


# =============================================================================
# HTML GENERATION
# =============================================================================

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>Astrophotography Sky Map</title>
    <script type="text/javascript" src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>
    <style>
        body { margin: 0; padding: 0; background-color: #111; color: white; font-family: sans-serif; }
        #aladin-lite-div { width: 100vw; height: 100vh; }
        #controls {
            position: absolute; top: 90px; left: 10px; z-index: 100;
            background: rgba(0,0,0,0.8); padding: 12px; border-radius: 6px;
            min-width: 220px; max-height: calc(100vh - 110px); overflow-y: auto;
        }
        #controls h3 { margin: 0; font-size: 14px; }
        #controls-header { display: flex; align-items: center; justify-content: space-between; gap: 8px; cursor: move; touch-action: none; user-select: none; }
        #controls-collapse { background: none; border: none; color: #aaa; cursor: pointer; font-size: 14px; padding: 0 4px; }
        #controls-collapse:hover { color: white; }
        #controls-body { margin-top: 6px; }
        #controls p  { margin: 0 0 12px 0; font-size: 12px; color: #aaa; }
        #legend { border-top: 1px solid #444; padding-top: 8px; font-size: 13px; }
        .control-group { margin-bottom: 10px; }
        .control-group label { display: block; font-size: 11px; color: #aaa; margin-bottom: 3px; }
        .survey-select { width: 100%; background: #222; color: white; border: 1px solid #444; padding: 5px; border-radius: 4px; font-size: 12px; }
        .search-container { display: flex; gap: 4px; }
        .search-input { flex-grow: 1; background: #222; color: white; border: 1px solid #444; padding: 5px; border-radius: 4px; font-size: 12px; }
        .search-btn { background: #444; color: white; border: none; padding: 5px 10px; border-radius: 4px; cursor: pointer; font-size: 12px; }
        .search-btn:hover { background: #555; }
        .checkbox-container { display: flex; align-items: center; gap: 6px; font-size: 12px; margin-top: 4px; cursor: pointer; color: white !important; }
        .checkbox-container input { cursor: pointer; }
        .legend-row { display: flex; align-items: center; gap: 6px; margin: 4px 0; font-size: 13px; cursor: pointer; }
        .legend-row .swatch { width: 14px; height: 14px; border-radius: 3px; display: inline-block; flex-shrink: 0; }

        #details-panel {
            position: absolute; top: 10px; right: 10px; z-index: 150;
            background: rgba(0,0,0,0.85); padding: 14px; border-radius: 6px;
            width: 270px; max-height: 80vh; overflow-y: auto;
            display: none; font-size: 13px; box-shadow: 0 4px 12px rgba(0,0,0,0.6);
        }
        #details-panel h3 { margin: 0 0 6px 0; font-size: 16px; }
        .details-close { float: right; cursor: pointer; color: #aaa; font-size: 14px; padding: 0 2px; }
        .details-close:hover { color: white; }
        .muted { color: #aaa; font-size: 12px; margin: 4px 0; }
        .stat { margin: 8px 0; }
        .badge { display: inline-block; background: #553300; color: #ffaa44; font-size: 10px; font-weight: bold; padding: 2px 6px; border-radius: 3px; margin-bottom: 4px; }
        table.filters { width: 100%; border-collapse: collapse; font-size: 12px; margin: 6px 0; }
        table.filters th { text-align: left; color: #888; font-weight: normal; border-bottom: 1px solid #444; padding: 2px 4px; }
        table.filters td { padding: 2px 4px; border-bottom: 1px solid #2a2a2a; }

        .ref-data { margin-top: 10px; border-top: 1px solid #444; padding-top: 8px; }
        .ref-head { font-size: 12px; font-weight: bold; color: #cdd; margin-bottom: 4px; }
        table.ref-table { width: 100%; border-collapse: collapse; font-size: 12px; }
        table.ref-table td { padding: 2px 4px; border-bottom: 1px solid #2a2a2a; vertical-align: top; color: #ddd; }
        table.ref-table td.rk { color: #8aa; white-space: nowrap; width: 38%; }
        .ref-links { margin-top: 6px; font-size: 12px; }
        .ref-links a { color: #6cf; text-decoration: none; }
        .ref-links a:hover { text-decoration: underline; }

        #toast {
            position: absolute; bottom: 20px; left: 50%; transform: translateX(-50%);
            background: rgba(0, 150, 255, 0.9); color: white; padding: 10px 20px;
            border-radius: 20px; font-size: 13px; font-weight: bold; z-index: 200;
            opacity: 0; transition: opacity 0.3s ease; pointer-events: none;
            box-shadow: 0 4px 10px rgba(0,0,0,0.5);
        }

        #folders {
            position: absolute; z-index: 100;
            background: rgba(0,0,0,0.8); padding: 12px; border-radius: 6px;
            width: 300px; max-height: calc(100vh - 110px); overflow: hidden;
            display: flex; flex-direction: column;
        }
        #folders-header { display: flex; align-items: center; justify-content: space-between; gap: 8px; cursor: move; touch-action: none; user-select: none; flex-shrink: 0; }
        #folders-header h3 { margin: 0; font-size: 14px; }
        #folders-collapse { background: none; border: none; color: #aaa; cursor: pointer; font-size: 14px; padding: 0 4px; }
        #folders-collapse:hover { color: white; }
        #folders-body { display: flex; flex-direction: column; min-height: 0; flex: 1 1 auto; overflow: hidden; }
        #folders-tools { display: flex; align-items: center; gap: 8px; margin: 8px 0; font-size: 11px; color: #aaa; flex-shrink: 0; }
        #folders-tools input[type=color] { width: 26px; height: 22px; padding: 0; border: 1px solid #444; background: #222; border-radius: 4px; cursor: pointer; }
        #rescan-bar { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; margin-bottom: 8px; font-size: 11px; color: #aaa; flex-shrink: 0; }
        #rescan-bar .rescan-btn { background: #2a3a2a; border: 1px solid #3c5c3c; color: #cfe; padding: 3px 8px; border-radius: 4px; cursor: pointer; font-size: 11px; }
        #rescan-bar .rescan-btn:hover { background: #345234; }
        #rescan-bar .rescan-btn:disabled { opacity: 0.5; cursor: default; }
        #folder-tree { flex: 1 1 auto; min-height: 0; overflow-y: auto; font-size: 12px; line-height: 1.7; }
        .folder-row { display: flex; align-items: center; gap: 4px; white-space: nowrap; cursor: pointer; border-radius: 3px; padding: 0 2px; }
        .folder-row:hover { background: rgba(255,255,255,0.08); }
        .folder-row.selected { background: rgba(255,0,255,0.18); }
        .folder-row.flash { animation: flashRow 2s ease-out; }
        @keyframes flashRow { 0% { background: rgba(255,230,0,0.5); } 100% { background: transparent; } }
        .folder-twirl { width: 12px; flex-shrink: 0; color: #888; text-align: center; }
        .folder-name { overflow: hidden; text-overflow: ellipsis; }
        .folder-count { color: #777; font-size: 10px; margin-left: 4px; flex-shrink: 0; }
        .file-row .folder-name { color: #bcd; }
        .file-icon { width: 12px; flex-shrink: 0; color: #668; text-align: center; font-size: 10px; }

        #header-panel {
            position: absolute; top: 70px; left: 50%; transform: translateX(-50%);
            z-index: 250; background: rgba(10,10,14,0.96); border: 1px solid #444;
            border-radius: 6px; width: 540px; max-width: 92vw; max-height: 78vh;
            display: none; flex-direction: column; box-shadow: 0 6px 24px rgba(0,0,0,0.7);
        }
        #header-panel-bar { display: flex; align-items: center; justify-content: space-between; gap: 8px; padding: 8px 12px; border-bottom: 1px solid #333; cursor: move; touch-action: none; user-select: none; flex-shrink: 0; }
        #header-panel-bar .title { font-size: 13px; font-weight: bold; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
        #header-panel-close { background: none; border: none; color: #aaa; cursor: pointer; font-size: 16px; padding: 0 4px; flex-shrink: 0; }
        #header-panel-close:hover { color: white; }
        #header-panel pre { margin: 0; padding: 10px 12px; overflow: auto; font-size: 11px; line-height: 1.45; white-space: pre; color: #cdd; flex: 1 1 auto; min-height: 0; }

        /* menu bar */
        #menubar {
            position: absolute; top: 0; left: 0; right: 0; z-index: 300; height: 36px;
            display: flex; align-items: center; gap: 6px; padding: 0 10px;
            background: rgba(0,0,0,0.85); border-bottom: 1px solid #333;
        }
        .menubar-title { font-size: 13px; font-weight: bold; margin-right: 8px; }
        .menu-btn { background: #2b2b2b; color: #ccc; border: 1px solid #444; padding: 4px 10px;
                    border-radius: 4px; cursor: pointer; font-size: 12px; }
        .menu-btn:hover { background: #3a3a3a; color: #fff; }
        .menu-btn.active { background: #2a4a6a; border-color: #3a6a9a; color: #fff; }

        /* generic draggable panel (stats/catalog/plan) */
        .panel {
            position: absolute; z-index: 140; display: none; flex-direction: column;
            background: rgba(0,0,0,0.85); border: 1px solid #444; border-radius: 6px;
            width: 320px; max-height: calc(100vh - 90px); box-shadow: 0 4px 12px rgba(0,0,0,0.6);
        }
        .panel.open { display: flex; }
        .panel-head { display: flex; align-items: center; justify-content: space-between;
                      padding: 8px 12px; border-bottom: 1px solid #333; cursor: move;
                      touch-action: none; user-select: none; flex-shrink: 0; }
        .panel-head h3 { margin: 0; font-size: 14px; }
        .panel-x { background: none; border: none; color: #aaa; cursor: pointer; font-size: 14px; }
        .panel-x:hover { color: #fff; }
        .panel-body { padding: 10px 12px; overflow-y: auto; font-size: 12px; }

        /* statistics */
        .stat-big { display: flex; gap: 14px; flex-wrap: wrap; margin-bottom: 10px; }
        .stat-big .num { font-size: 20px; font-weight: bold; color: #cfe; }
        .stat-big .lbl { font-size: 10px; color: #999; text-transform: uppercase; }
        .stat-sec { font-size: 12px; font-weight: bold; color: #cdd; margin: 12px 0 4px; }
        table.kv { width: 100%; border-collapse: collapse; }
        table.kv td { padding: 2px 4px; border-bottom: 1px solid #2a2a2a; }
        table.kv td.r { text-align: right; color: #9bd; white-space: nowrap; }
        .bar { background: #2a3a55; height: 12px; border-radius: 2px; display: inline-block; vertical-align: middle; }
        .barwrap { display: flex; align-items: center; gap: 6px; }

        /* catalog */
        .cat-toggles { display: flex; gap: 12px; margin-bottom: 8px; }
        .cat-prog { font-size: 12px; margin: 4px 0; }
        .cat-prog .pct { color: #6f6; font-weight: bold; }
        #catalog-list { max-height: 46vh; overflow-y: auto; }
        .cat-row { display: flex; align-items: center; gap: 6px; padding: 1px 2px; cursor: pointer; border-radius: 3px; }
        .cat-row:hover { background: rgba(255,255,255,0.08); }
        .cat-row .ck { width: 14px; text-align: center; flex-shrink: 0; }
        .cat-row.done { color: #7d7; }
        .cat-row.todo { color: #999; }

        /* plan */
        #plan-alt { width: 100%; height: 120px; background: #0a0a12; border-radius: 4px; display: block; }
    </style>
</head>
<body>
    <div id="menubar">
        <span class="menubar-title">🔭 Sky Map</span>
        <button class="menu-btn" data-panel="controls">Controls</button>
        <button class="menu-btn" data-panel="folders">Folders</button>
        <button class="menu-btn" data-panel="stats">Statistics</button>
        <button class="menu-btn" data-panel="catalog">Catalog</button>
        <button class="menu-btn" id="plan-mode-btn" title="When on, clicking empty sky opens a planning popup for that point">Plan mode</button>
        <button class="menu-btn" id="poster-btn" title="Download the current map view as a PNG">Export poster</button>
    </div>

    <div id="controls">
        <div id="controls-header" title="Drag to move · click to collapse/expand">
            <h3>Astrophotography Archive Map</h3>
            <button id="controls-collapse">▾</button>
        </div>
        <div id="controls-body">
        <p id="target-count"></p>

        <div class="control-group">
            <label>Search Target (e.g. 6044, abell 39, M42):</label>
            <div class="search-container">
                <input type="text" id="search-box" class="search-input" placeholder="Go to target...">
                <button id="search-submit" class="search-btn">Go</button>
            </div>
        </div>

        <div class="control-group">
            <label>Background Survey:</label>
            <select id="survey-selector" class="survey-select">
                <option value="P/DSS2/color">DSS2 Color (Global Optical)</option>
                <option value="P/DSS2/red">DSS2 Red (Best H-Alpha Nebula Contrast)</option>
                <option value="CDS/P/DESI-Legacy-Surveys/DR10/color">DESI Legacy DR10 (Ultra High-Res Galaxy/Deep Sky)</option>
                <option value="CDS/P/PanSTARRS/DR1/color-i-r-g">PanSTARRS DR1 Color (High-Res Northern Sky)</option>
                <option value="P/Mellinger/color">Mellinger All-Sky Optical (Wide-Field/Milky Way)</option>
                <option value="CDS/P/Finkbeiner">Finkbeiner H-Alpha (Global 100% Full Sky Nebula)</option>
                <option value="P/2MASS/color">2MASS Infrared (Star Crowding Reducer)</option>
                <option value="P/allWISE/color">AllWISE Infrared (Dust & Gas Profiles)</option>
            </select>
        </div>

        <div class="control-group">
            <label class="checkbox-container">
                <input type="checkbox" id="grid-toggle"> RA/DEC grid lines
            </label>
            <label class="checkbox-container">
                <input type="checkbox" id="constellation-toggle"> Constellation lines
            </label>
            <label class="checkbox-container">
                <input type="checkbox" id="mosaic-toggle" checked> Show mosaic frames
            </label>
            <label class="checkbox-container" title="Draws an outline around the part of the sky that rises above the minimum altitude during darkness tonight. Site location comes from your image headers (or the Python script).">
                <input type="checkbox" id="tonight-toggle"> Outline sky available tonight
            </label>
            <label class="checkbox-container" title="Mirror the map east–west">
                <input type="checkbox" id="flip-h-toggle"> Flip horizontal
            </label>
            <label class="checkbox-container" title="Mirror the map north–south">
                <input type="checkbox" id="flip-v-toggle"> Flip vertical
            </label>
        </div>

        <div id="legend"></div>

        <button id="export-csv" class="search-btn" style="width:100%;margin-top:10px;">Export target list (CSV)</button>
        </div>
    </div>

    <div id="folders">
        <div id="folders-header" title="Drag to move · click to collapse/expand">
            <h3>Image Folders</h3>
            <button id="folders-collapse">▾</button>
        </div>
        <div id="folders-body">
            <div id="folders-tools">
                <span title="Outline color for highlighted folders">Highlight:</span>
                <input type="color" id="highlight-color">
                <button id="folders-clear" class="search-btn" style="padding:3px 8px;">Clear</button>
                <button id="tree-collapse" class="search-btn" style="padding:3px 8px;" title="Collapse all folders">Collapse</button>
            </div>
            <div id="rescan-bar" title="Re-scan disk for added/deleted files, then rebuild the map"></div>
            <div id="folder-tree"></div>
        </div>
    </div>

    <div id="header-panel">
        <div id="header-panel-bar">
            <span class="title" id="header-panel-title"></span>
            <button id="header-panel-close" title="Close">✕</button>
        </div>
        <pre id="header-panel-pre"></pre>
    </div>

    <div id="stats" class="panel">
        <div class="panel-head" data-drag="stats"><h3>Statistics</h3><button class="panel-x" data-panel="stats">✕</button></div>
        <div class="panel-body" id="stats-body"></div>
    </div>

    <div id="catalog" class="panel">
        <div class="panel-head" data-drag="catalog"><h3>Catalog completion</h3><button class="panel-x" data-panel="catalog">✕</button></div>
        <div class="panel-body" id="catalog-body"></div>
    </div>

    <div id="plan" class="panel">
        <div class="panel-head" data-drag="plan"><h3>Plan point</h3><button class="panel-x" data-panel="plan">✕</button></div>
        <div class="panel-body" id="plan-body"></div>
    </div>

    <div id="details-panel"></div>
    <div id="toast"></div>
    <div id="aladin-lite-div"></div>

    <script type="text/javascript">
        const DATA = __DATA_JSON__;
        const FOLDERS = __FOLDERS_JSON__;
        const STATS = __STATS_JSON__;
        const SETTINGS = __SETTINGS_JSON__;

        // --- USER-DEFINED SEARCH ALIASES ---
        const ALIASES = {
            'abell 39': 'PN A66 39',
            'abell 33': 'PN A66 33'
        };

        const TONIGHT_COLOR = '#ffaa00';
        const LS_KEY = 'sky_mapper_ui_v2';

        A.init.then(() => {
            console.log("✅ Aladin Lite initialized");

            const aladin = A.aladin('#aladin-lite-div', {
                survey: "P/DSS2/color",
                fov: SETTINGS.initialFov,
                showReticle: true,
                showZoomControl: true,
                showFullscreenControl: true,
                showCooGrid: false
            });
            aladin.gotoRaDec(SETTINGS.initialRa, SETTINGS.initialDec);
            window.aladin = aladin;   // console access for debugging

            const els = {
                grid:           document.getElementById('grid-toggle'),
                constellations: document.getElementById('constellation-toggle'),
                mosaics:        document.getElementById('mosaic-toggle'),
                tonight:        document.getElementById('tonight-toggle'),
                flipH:          document.getElementById('flip-h-toggle'),
                flipV:          document.getElementById('flip-v-toggle'),
                survey:         document.getElementById('survey-selector'),
                count:          document.getElementById('target-count'),
                panel:          document.getElementById('details-panel'),
            };

            // ---------- small helpers ----------
            const esc = s => String(s).replace(/[&<>"]/g,
                ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[ch]));

            let toastTimer = null;
            function toastMsg(msg) {
                const t = document.getElementById('toast');
                t.innerText = msg;
                t.style.opacity = '1';
                clearTimeout(toastTimer);
                toastTimer = setTimeout(() => { t.style.opacity = '0'; }, 2200);
            }

            function fmtExposure(sec) {
                if (!sec) return '—';
                const h = Math.floor(sec / 3600);
                const m = Math.round((sec % 3600) / 60);
                return h ? h + 'h ' + String(m).padStart(2, '0') + 'm' : m + 'm';
            }
            function raToHms(ra) {
                const totalH = ra / 15;
                const h = Math.floor(totalH);
                const m = Math.floor((totalH - h) * 60);
                const s = ((totalH - h) * 60 - m) * 60;
                return h + 'h ' + String(m).padStart(2, '0') + 'm ' + String(Math.round(s)).padStart(2, '0') + 's';
            }
            function decToDms(dec) {
                const sign = dec < 0 ? '-' : '+';
                const a = Math.abs(dec);
                const d = Math.floor(a);
                const m = Math.floor((a - d) * 60);
                const s = ((a - d) * 60 - m) * 60;
                return sign + d + '° ' + String(m).padStart(2, '0') + "' " + String(Math.round(s)).padStart(2, '0') + '"';
            }

            // ---------- per-rig layers (plus one gray layer for dimmed targets) ----------
            const layers = {};
            for (const color of [...new Set(DATA.clusters.map(c => c.color))]) {
                const overlay = A.graphicOverlay({ color: color, lineWidth: 2 });
                aladin.addOverlay(overlay);
                const catalog = A.catalog({
                    name: 'Labels ' + color, sourceSize: 0, color: color,
                    displayLabel: true, labelColumn: 'name',
                    labelColor: color, labelFont: 'bold 15px sans-serif'
                });
                aladin.addCatalog(catalog);
                layers[color] = { overlay: overlay, catalog: catalog };
            }
            function clearLayer(l) {
                for (const part of [l.overlay, l.catalog]) {
                    if (part && typeof part.removeAll === 'function') part.removeAll();
                }
            }

            // ---------- rig legend with show/hide checkboxes ----------
            const rigEnabled = {}, rigCheckboxes = {};
            const legendDiv = document.getElementById('legend');
            for (const rig of DATA.rigs) {
                rigEnabled[rig.color] = true;
                const row = document.createElement('label');
                row.className = 'legend-row';
                row.innerHTML = '<input type="checkbox" checked> <span class="swatch" style="background:'
                    + rig.color + '"></span> <span>' + esc(rig.label) + '</span>';
                const cb = row.querySelector('input');
                cb.addEventListener('change', () => {
                    rigEnabled[rig.color] = cb.checked;
                    saveState();
                    rebuild();
                });
                rigCheckboxes[rig.color] = cb;
                legendDiv.appendChild(row);
            }

            // ---------- "imageable tonight" altitude math ----------
            const D2R = Math.PI / 180, R2D = 180 / Math.PI;
            function julianDate(d) { return d.getTime() / 86400000 + 2440587.5; }
            function gmstDeg(jd) {
                const d = jd - 2451545.0;
                return (((280.46061837 + 360.98564736629 * d) % 360) + 360) % 360;
            }
            function sunRaDec(jd) {
                const n = jd - 2451545.0;
                const L = (280.460 + 0.9856474 * n) % 360;
                const g = ((357.528 + 0.9856003 * n) % 360) * D2R;
                const lam = (L + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * D2R;
                const eps = (23.439 - 0.0000004 * n) * D2R;
                return {
                    ra: ((Math.atan2(Math.cos(eps) * Math.sin(lam), Math.cos(lam)) * R2D) + 360) % 360,
                    dec: Math.asin(Math.sin(eps) * Math.sin(lam)) * R2D
                };
            }
            function altDeg(raDeg, decDeg, jd) {
                const ha = (gmstDeg(jd) + SETTINGS.lon - raDeg) * D2R;
                const lat = SETTINGS.lat * D2R, dec = decDeg * D2R;
                return Math.asin(Math.sin(lat) * Math.sin(dec)
                     + Math.cos(lat) * Math.cos(dec) * Math.cos(ha)) * R2D;
            }
            function darkWindowJDs() {
                // Tonight = local noon to local noon. Sample every 10 min, keep
                // times when the sun is below -18° (fallback -12°).
                const start = new Date();
                start.setHours(12, 0, 0, 0);
                if (Date.now() < start.getTime()) start.setDate(start.getDate() - 1);
                const astro = [], nautical = [];
                for (let m = 0; m < 24 * 60; m += 10) {
                    const jd = julianDate(new Date(start.getTime() + m * 60000));
                    const s = sunRaDec(jd);
                    const alt = altDeg(s.ra, s.dec, jd);
                    if (alt < -18) astro.push(jd);
                    if (alt < -12) nautical.push(jd);
                }
                if (!astro.length && nautical.length) toastMsg('No astronomical dark tonight — using nautical twilight');
                return astro.length ? astro : nautical;
            }

            // The sky available tonight is the union of "altitude ≥ minAlt" caps
            // (radius 90−minAlt around the zenith) as the zenith sweeps along
            // dec = site latitude from dusk LST to dawn LST. Its boundary is a
            // single closed contour when a celestial pole lies inside the region
            // (true whenever |lat| ≥ 90 − capRadius, e.g. lat ≥ 30° for a 30°
            // altitude limit) — that is the case drawn here.
            let tonightOverlay = null;
            function ensureTonightOutline() {
                if (tonightOverlay) { tonightOverlay.show(); return; }
                const win = darkWindowJDs();
                if (!win.length) {
                    toastMsg('No darkness tonight at your site');
                    els.tonight.checked = false;
                    return;
                }
                const lst0 = ((gmstDeg(win[0]) + SETTINGS.lon) % 360 + 360) % 360;
                const lst1 = ((gmstDeg(win[win.length - 1]) + SETTINGS.lon) % 360 + 360) % 360;
                const sweep = ((lst1 - lst0) % 360 + 360) % 360;
                const capR = 90 - SETTINGS.minAlt;
                const poleSign = SETTINGS.lat >= 0 ? 1 : -1;
                const phi = SETTINGS.lat * D2R;
                const cosR = Math.cos(capR * D2R);

                // Southernmost (north sites) / northernmost (south sites) dec where
                // the cap centred (deltaRA away from the zenith track end) crosses
                // this meridian: solve sinφ·sinδ + cosφ·cosΔα·cosδ = cos(capR).
                function capEdgeDec(deltaRaDeg) {
                    const Ac = Math.sin(phi);
                    const Bc = Math.cos(phi) * Math.cos(deltaRaDeg * D2R);
                    const Rm = Math.sqrt(Ac * Ac + Bc * Bc);
                    if (Rm < 1e-9 || Math.abs(cosR / Rm) > 1) return null;
                    const psi = Math.atan2(Bc, Ac);
                    const sols = [Math.asin(cosR / Rm) - psi, Math.PI - Math.asin(cosR / Rm) - psi]
                        .filter(d => d >= -Math.PI / 2 && d <= Math.PI / 2)
                        .map(d => d * R2D);
                    if (!sols.length) return null;
                    return poleSign > 0 ? Math.min(...sols) : Math.max(...sols);
                }

                const pts = [];
                for (let a = 0; a <= 360; a += 2) {
                    const off = ((a - lst0) % 360 + 360) % 360;
                    let dec;
                    if (off <= sweep) {
                        dec = SETTINGS.lat - poleSign * capR;
                    } else {
                        const offEnd = ((a - lst1) % 360 + 360) % 360;
                        const cand = [
                            capEdgeDec(Math.min(off, 360 - off)),
                            capEdgeDec(Math.min(offEnd, 360 - offEnd)),
                        ].filter(v => v !== null);
                        dec = cand.length
                            ? (poleSign > 0 ? Math.min(...cand) : Math.max(...cand))
                            : poleSign * 89.9;
                    }
                    pts.push([a, Math.max(-89.9, Math.min(89.9, dec))]);
                }
                pts[pts.length - 1] = [360, pts[0][1]];

                tonightOverlay = A.graphicOverlay({ color: TONIGHT_COLOR, lineWidth: 3 });
                aladin.addOverlay(tonightOverlay);
                tonightOverlay.add(A.polyline(pts, { color: TONIGHT_COLOR, lineWidth: 3 }));
                toastMsg('Orange outline = sky above ' + SETTINGS.minAlt + '° altitude sometime tonight');
            }

            // ---------- (re)draw all footprints + labels from current UI state ----------
            function clusterVisible(c) {
                if (rigEnabled[c.color] === false) return false;
                if (c.mosaic && !els.mosaics.checked) return false;
                return true;
            }
            function rebuild() {
                for (const l of Object.values(layers)) clearLayer(l);
                let shown = 0;
                for (const c of DATA.clusters) {
                    if (!clusterVisible(c)) continue;
                    const layer = layers[c.color];
                    layer.overlay.addFootprints([A.polygon(c.corners)]);
                    layer.catalog.addSources([A.source(c.ra, c.label_dec, { name: c.label, cid: c.id })]);
                    shown++;
                }
                els.count.innerText = 'Showing ' + shown + ' of ' + DATA.clusters.length + ' targets';
            }

            // ---------- label declutter: hide labels when zoomed way out ----------
            function updateLabelVisibility() {
                let fov = SETTINGS.initialFov;
                try { fov = aladin.getFov()[0]; } catch (e) {}
                const hide = fov > SETTINGS.labelHideFov;
                for (const l of Object.values(layers)) {
                    try { hide ? l.catalog.hide() : l.catalog.show(); } catch (e) {}
                }
            }
            try { aladin.on('zoomChanged', updateLabelVisibility); } catch (e) {}

            // ---------- details panel ----------
            function findClusterAt(ra, dec) {
                let best = null, bestArea = Infinity;
                for (const c of DATA.clusters) {
                    if (!clusterVisible(c)) continue;
                    const r = c.rect;
                    if (dec < r[2] || dec > r[3]) continue;
                    const width = r[1] - r[0];
                    const delta = (((ra - r[0]) % 360) + 360) % 360;   // RA-wrap safe
                    if (delta > width) continue;
                    const area = width * (r[3] - r[2]);
                    if (area < bestArea) { best = c; bestArea = area; }
                }
                return best;
            }
            function hideDetails() { els.panel.style.display = 'none'; }
            function showDetails(c) {
                const rows = Object.entries(c.filters)
                    .sort((a, b) => b[1].seconds - a[1].seconds)
                    .map(([f, v]) => '<tr><td>' + esc(f) + '</td><td>' + v.frames + '</td><td>'
                                   + fmtExposure(v.seconds) + '</td></tr>')
                    .join('');
                const dates = c.date_min
                    ? (c.date_min === c.date_max ? c.date_min : c.date_min + ' → ' + c.date_max)
                    : 'unknown';
                els.panel.innerHTML =
                    '<span class="details-close" title="Close">✕</span>' +
                    '<h3 style="color:' + c.color + '">' + esc(c.name) + '</h3>' +
                    (c.mosaic ? '<span class="badge">MOSAIC</span>' : '') +
                    (c.alt_names.length ? '<div class="muted">also: ' + esc(c.alt_names.join(', ')) + '</div>' : '') +
                    '<div class="stat"><b>' + c.count + '</b> frames · <b>' + fmtExposure(c.total_seconds) + '</b> total</div>' +
                    '<table class="filters"><tr><th>Filter</th><th>Frames</th><th>Time</th></tr>' + rows + '</table>' +
                    '<div class="muted">Dates: ' + esc(dates) + '</div>' +
                    (c.rigs.length ? '<div class="muted">Rig: ' + esc(c.rigs.join(', ')) + '</div>' : '') +
                    '<div class="muted">' + raToHms(c.ra) + '   ' + decToDms(c.dec) + '</div>' +
                    '<button class="search-btn goto-btn" style="margin-top:8px;width:100%;">Center on target</button>' +
                    '<canvas id="det-alt" style="width:100%;height:90px;margin-top:8px;display:none;"></canvas>' +
                    '<div id="det-alt-info" class="muted"></div>' +
                    '<div id="ref-data" class="ref-data muted">Loading reference data…</div>';
                els.panel.dataset.cid = c.id;
                els.panel.querySelector('.details-close').addEventListener('click', hideDetails);
                els.panel.querySelector('.goto-btn').addEventListener('click', () => {
                    aladin.gotoRaDec(c.ra, c.dec);
                    aladin.setFov(4);
                });
                els.panel.style.display = 'block';
                if (SETTINGS.lat !== null) {
                    const cv = document.getElementById('det-alt');
                    cv.style.display = 'block'; cv.width = cv.clientWidth || 240; cv.height = 90;
                    const r = altitudeCurve(cv, c.ra, c.dec);
                    if (r) document.getElementById('det-alt-info').innerText =
                        'Max altitude tonight ' + r.maxAlt + '° at ' + hhmm(r.maxTime);
                }
                loadReferenceData(c);
            }

            // ---------- object reference data (SIMBAD via /api/objectinfo) ----------
            function refRow(label, value) {
                return value ? '<tr><td class="rk">' + esc(label) + '</td><td>' + value + '</td></tr>' : '';
            }
            function renderRefData(c, info) {
                const box = els.panel.querySelector('#ref-data');
                if (!box || els.panel.dataset.cid != c.id) return;   // panel moved on
                if (!info || !info.ok) { box.textContent = 'Reference data unavailable.'; return; }

                let rows = '';
                if (info.resolved) {
                    let typ = esc(info.otype || info.otype_code || '');
                    if (info.morph) typ += ' <span style="color:#888">(' + esc(info.morph) + ')</span>';
                    rows += refRow('Type', typ);
                    rows += refRow('Constellation', esc(info.constellation || ''));
                    if (info.mag) {
                        const m = ['V', 'B'].filter(b => info.mag[b]).map(b => b + ' ' + esc(info.mag[b])).join(', ');
                        rows += refRow('Magnitude', m);
                    }
                    rows += refRow('Apparent size', esc(info.size || ''));
                    rows += refRow('Distance', esc(info.distance || ''));
                    rows += refRow('Redshift', info.z ? esc(info.z) : '');
                    rows += refRow('Radial vel.', info.rv ? esc(info.rv) + ' km/s' : '');
                    if (info.common && info.common.length)
                        rows += refRow('Also known as', esc(info.common.join(', ')));
                    rows += refRow('SIMBAD id', esc(info.main_id || ''));
                } else {
                    rows += refRow('Constellation', esc(info.constellation || ''));
                    rows += '<tr><td colspan="2" style="color:#888">Not resolved by SIMBAD' +
                            ' — try the links below.</td></tr>';
                }
                const links = [];
                if (info.simbad_url) links.push('<a href="' + info.simbad_url + '" target="_blank" rel="noopener">SIMBAD</a>');
                if (info.ned_url) links.push('<a href="' + info.ned_url + '" target="_blank" rel="noopener">NED</a>');

                box.innerHTML =
                    '<div class="ref-head">Reference data</div>' +
                    '<table class="ref-table">' + rows + '</table>' +
                    (links.length ? '<div class="ref-links">More: ' + links.join(' · ') + '</div>' : '');
            }
            async function loadReferenceData(c) {
                // reuse the search-box aliases (e.g. Abell 39 → the planetary nebula PN A66 39,
                // not the galaxy cluster SIMBAD returns for "Abell 39")
                const qName = ALIASES[(c.name || '').trim().toLowerCase()] || c.name;
                try {
                    const resp = await fetch('/api/objectinfo?name=' + encodeURIComponent(qName) +
                                             '&ra=' + c.ra + '&dec=' + c.dec);
                    renderRefData(c, await resp.json());
                } catch (e) {
                    const box = els.panel.querySelector('#ref-data');
                    if (box && els.panel.dataset.cid == c.id) box.textContent = 'Reference data unavailable.';
                }
            }

            // ---------- map click: open details, or plan popup in plan mode ----------
            aladin.on('click', (object) => {
                if (!object || object.ra === undefined || object.dec === undefined) return;
                const c = findClusterAt(object.ra, object.dec);
                if (c) { showDetails(c); }
                else { hideDetails(); if (planMode) openPlanPopup(object.ra, object.dec); }
                revealFoldersAt(object.ra, object.dec);
            });

            // ---------- constellation lines (lazy-loaded) ----------
            let constellationOverlay = null;
            async function ensureConstellations() {
                if (constellationOverlay) { constellationOverlay.show(); return; }
                let data = null;
                for (const url of [SETTINGS.constellationLocal, SETTINGS.constellationUrl]) {
                    if (!url) continue;
                    try {
                        const resp = await fetch(url);
                        if (resp.ok) { data = await resp.json(); break; }
                    } catch (e) { /* try next source */ }
                }
                if (!data || !data.features) {
                    toastMsg('Could not load constellation data');
                    els.constellations.checked = false;
                    return;
                }
                constellationOverlay = A.graphicOverlay({ color: '#557788', lineWidth: 1 });
                aladin.addOverlay(constellationOverlay);
                for (const feat of data.features) {
                    const geom = feat.geometry;
                    if (!geom) continue;
                    const lines = geom.type === 'MultiLineString' ? geom.coordinates : [geom.coordinates];
                    for (const line of lines) {
                        let seg = [], prevRa = null;
                        for (const pt of line) {
                            const ra = ((pt[0] % 360) + 360) % 360;
                            // break polylines that jump across the RA 0/360 seam
                            if (prevRa !== null && Math.abs(ra - prevRa) > 180) {
                                if (seg.length > 1) constellationOverlay.add(A.polyline(seg, { color: '#557788', lineWidth: 1 }));
                                seg = [];
                            }
                            seg.push([ra, pt[1]]);
                            prevRa = ra;
                        }
                        if (seg.length > 1) constellationOverlay.add(A.polyline(seg, { color: '#557788', lineWidth: 1 }));
                    }
                }
                if (!els.constellations.checked) constellationOverlay.hide();
                console.log('✅ Constellation lines loaded');
            }

            // ---------- CSV export ----------
            document.getElementById('export-csv').addEventListener('click', () => {
                const rows = [['name', 'ra_deg', 'dec_deg', 'ra_hms', 'dec_dms', 'frames',
                               'integration_hours', 'filters', 'mosaic', 'rigs', 'first_date', 'last_date']];
                for (const c of DATA.clusters) {
                    const filters = Object.entries(c.filters)
                        .map(([f, v]) => f + ': ' + fmtExposure(v.seconds)).join('; ');
                    rows.push([c.name, c.ra.toFixed(5), c.dec.toFixed(5), raToHms(c.ra), decToDms(c.dec),
                               c.count, (c.total_seconds / 3600).toFixed(2), filters,
                               c.mosaic ? 'yes' : 'no', c.rigs.join('; '),
                               c.date_min || '', c.date_max || '']);
                }
                const csv = rows.map(r => r.map(v => '"' + String(v).replace(/"/g, '""') + '"').join(',')).join('\n');
                const blob = new Blob([csv], { type: 'text/csv' });
                const a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = 'sky_targets.csv';
                a.click();
                URL.revokeObjectURL(a.href);
                toastMsg('Exported ' + DATA.clusters.length + ' targets');
            });

            // ---------- search with auto-prefix and aliases ----------
            const executeSearch = () => {
                let query = document.getElementById('search-box').value.trim().toLowerCase();
                if (!query) return;
                if (/^\d+$/.test(query)) {
                    query = "ngc " + query;
                }
                query = ALIASES[query] || query;
                aladin.gotoObject(query, {
                    error: () => console.log("⚠ Could not resolve target: " + query)
                });
            };
            document.getElementById('search-submit').addEventListener('click', executeSearch);
            document.getElementById('search-box').addEventListener('keypress', (e) => {
                if (e.key === 'Enter') executeSearch();
            });

            // ---------- collapsible, draggable control panel ----------
            const controls = document.getElementById('controls');
            const controlsBody = document.getElementById('controls-body');
            const controlsHeader = document.getElementById('controls-header');
            const collapseBtn = document.getElementById('controls-collapse');
            function setCollapsed(collapsed) {
                controlsBody.style.display = collapsed ? 'none' : 'block';
                collapseBtn.innerText = collapsed ? '▸' : '▾';
            }
            function placePanel(x, y) {
                x = Math.max(0, Math.min(x, window.innerWidth - 80));
                y = Math.max(0, Math.min(y, window.innerHeight - 40));
                controls.style.left = x + 'px';
                controls.style.top = y + 'px';
                controls.style.maxHeight = (window.innerHeight - y - 20) + 'px';
            }
            placePanel(SETTINGS.panelX, SETTINGS.panelY);

            // Drag by the title bar; a press that barely moves counts as a
            // collapse/expand click instead.
            let drag = null;
            controlsHeader.addEventListener('pointerdown', (e) => {
                drag = { x0: e.clientX, y0: e.clientY,
                         left: controls.offsetLeft, top: controls.offsetTop, moved: false };
                try { controlsHeader.setPointerCapture(e.pointerId); } catch (err) {}
                e.preventDefault();
            });
            controlsHeader.addEventListener('pointermove', (e) => {
                if (!drag) return;
                const dx = e.clientX - drag.x0, dy = e.clientY - drag.y0;
                if (Math.abs(dx) + Math.abs(dy) > 4) drag.moved = true;
                if (drag.moved) placePanel(drag.left + dx, drag.top + dy);
            });
            controlsHeader.addEventListener('pointerup', () => {
                if (!drag) return;
                const moved = drag.moved;
                drag = null;
                if (!moved) setCollapsed(controlsBody.style.display !== 'none');
                saveState();
            });

            // ---------- control wiring ----------
            els.grid.addEventListener('change', (e) => {
                if (e.target.checked) aladin.showCooGrid(); else aladin.hideCooGrid();
                saveState();
            });
            els.survey.addEventListener('change', (e) => {
                aladin.setBaseImageLayer(aladin.createImageSurvey(e.target.value));
                saveState();
            });
            els.constellations.addEventListener('change', (e) => {
                if (e.target.checked) ensureConstellations();
                else if (constellationOverlay) constellationOverlay.hide();
                saveState();
            });
            els.mosaics.addEventListener('change', () => { saveState(); rebuild(); });
            // A vertical mirror = horizontal mirror + 180° rotation, so the two
            // checkboxes combine into a longitude-reverse plus optional rotation.
            // setRotation(0) is rejected as falsy by Aladin, hence 360 for upright.
            function applyFlips() {
                aladin.reverseLongitude(els.flipH.checked !== els.flipV.checked);
                aladin.setRotation(els.flipV.checked ? 180 : 360);
            }
            els.flipH.addEventListener('change', () => { applyFlips(); saveState(); });
            els.flipV.addEventListener('change', () => { applyFlips(); saveState(); });
            els.tonight.addEventListener('change', () => {
                if (els.tonight.checked && SETTINGS.lat === null) {
                    toastMsg('No site location — set OBSERVER_LAT/LON in the script');
                    els.tonight.checked = false;
                    return;
                }
                if (els.tonight.checked) ensureTonightOutline();
                else if (tonightOverlay) tonightOverlay.hide();
                saveState();
            });

            // ---------- folder browser: click a folder to outline its frames ----------
            const treeDiv = document.getElementById('folder-tree');
            const colorInput = document.getElementById('highlight-color');
            const foldersPanel = document.getElementById('folders');
            const foldersBody = document.getElementById('folders-body');
            const foldersHeader = document.getElementById('folders-header');
            const foldersCollapseBtn = document.getElementById('folders-collapse');
            colorInput.value = SETTINGS.highlightColor;

            const tree = FOLDERS.tree;
            const kids = {};
            tree.forEach((n, id) => { (kids[n.parent] = kids[n.parent] || []).push(id); });
            for (const k in kids) kids[k].sort((a, b) => tree[a].label.localeCompare(tree[b].label));
            const roots = (kids['null'] || []).slice();

            // Auto-expand the rig roots only when there's a single rig; with
            // several, start collapsed so every rig is visible at the top.
            const expanded = new Set(roots.length === 1 ? roots : []);
            const selected = new Set();           // selected folder node ids
            const selectedFiles = new Map();      // file path -> footprint index

            const highlightOverlay = A.graphicOverlay({ color: SETTINGS.highlightColor, lineWidth: 3 });
            aladin.addOverlay(highlightOverlay);

            function subtreeFps(id, out) {
                for (const f of tree[id].fp) out.add(f);
                for (const c of (kids[id] || [])) subtreeFps(c, out);
            }
            function forceRedraw() {
                // removeAll() alone doesn't repaint an emptied overlay, so the map
                // would keep showing cleared outlines — nudge Aladin to redraw.
                try { aladin.view.requestRedraw(); }
                catch (e) { try { aladin.view.forceRedraw(); } catch (e2) {} }
            }
            function refreshHighlights() {
                highlightOverlay.removeAll();
                const col = colorInput.value;
                const out = new Set();
                for (const id of selected) subtreeFps(id, out);
                for (const fp of selectedFiles.values()) out.add(fp);
                for (const idx of out) {
                    const c = FOLDERS.footprints[idx];
                    highlightOverlay.add(A.polyline([...c, c[0]], { color: col, lineWidth: 3 }));
                }
                forceRedraw();
            }
            function centerOnFolder(id) {
                const out = new Set();
                subtreeFps(id, out);
                if (!out.size) return;
                let x = 0, y = 0, z = 0, n = 0;
                for (const idx of out) {
                    const c = FOLDERS.footprints[idx];
                    const ra = (c[0][0] + c[1][0] + c[2][0] + c[3][0]) / 4 * Math.PI / 180;
                    const dec = (c[0][1] + c[1][1] + c[2][1] + c[3][1]) / 4 * Math.PI / 180;
                    x += Math.cos(dec) * Math.cos(ra); y += Math.cos(dec) * Math.sin(ra); z += Math.sin(dec); n++;
                }
                const ra = ((Math.atan2(y, x) * 180 / Math.PI) % 360 + 360) % 360;
                const dec = Math.asin(z / Math.sqrt(x * x + y * y + z * z)) * 180 / Math.PI;
                aladin.gotoRaDec(ra, dec);
                try { if (aladin.getFov()[0] > 20) aladin.setFov(8); } catch (e) {}
            }

            function renderNode(id, depth, frag) {
                const n = tree[id];
                const row = document.createElement('div');
                row.className = 'folder-row' + (selected.has(id) ? ' selected' : '');
                row.dataset.nodeId = id;
                row.style.paddingLeft = (depth * 14) + 'px';
                const hasKids = !!(kids[id] && kids[id].length);
                const files = n.files || [];
                const expandable = hasKids || files.length;

                const twirl = document.createElement('span');
                twirl.className = 'folder-twirl';
                twirl.textContent = expandable ? (expanded.has(id) ? '▾' : '▸') : '';
                if (expandable) twirl.addEventListener('click', (e) => {
                    e.stopPropagation();
                    expanded.has(id) ? expanded.delete(id) : expanded.add(id);
                    render();
                });
                row.appendChild(twirl);

                const name = document.createElement('span');
                name.className = 'folder-name';
                name.textContent = n.label;
                row.appendChild(name);

                const count = document.createElement('span');
                count.className = 'folder-count';
                count.textContent = n.total;
                row.appendChild(count);

                row.addEventListener('click', () => {
                    if (selected.has(id)) { selected.delete(id); }
                    else { selected.add(id); centerOnFolder(id); }
                    refreshHighlights();
                    render();
                });
                frag.appendChild(row);

                if (expandable && expanded.has(id)) {
                    for (const c of kids[id] || []) renderNode(c, depth + 1, frag);
                    for (const [fname, fpIdx] of files) renderFile(n.adir, fname, fpIdx, depth + 1, frag);
                }
            }

            function renderFile(adir, fname, fpIdx, depth, frag) {
                const path = adir + '/' + fname;
                const row = document.createElement('div');
                row.className = 'folder-row file-row' + (selectedFiles.has(path) ? ' selected' : '');
                row.style.paddingLeft = (depth * 14) + 'px';

                const icon = document.createElement('span');
                icon.className = 'file-icon';
                icon.textContent = '·';
                row.appendChild(icon);

                const name = document.createElement('span');
                name.className = 'folder-name';
                name.textContent = fname;
                row.appendChild(name);

                row.addEventListener('click', () => {
                    if (selectedFiles.has(path)) selectedFiles.delete(path);
                    else { selectedFiles.set(path, fpIdx); centerOnFp(fpIdx); }
                    refreshHighlights();
                    showHeader(path, fname);
                    render();
                });
                frag.appendChild(row);
            }

            function centerOnFp(fpIdx) {
                const c = FOLDERS.footprints[fpIdx];
                const ra = (c[0][0] + c[1][0] + c[2][0] + c[3][0]) / 4;
                const dec = (c[0][1] + c[1][1] + c[2][1] + c[3][1]) / 4;
                aladin.gotoRaDec(ra, dec);
                try { if (aladin.getFov()[0] > 20) aladin.setFov(8); } catch (e) {}
            }
            function render() {
                const frag = document.createDocumentFragment();
                for (const r of roots) renderNode(r, 0, frag);
                treeDiv.innerHTML = '';
                treeDiv.appendChild(frag);
            }

            // ---------- map click → jump the browser to that target's folder(s) ----------
            // Reverse index: which leaf folders contain a file at each footprint.
            const fpToNodes = new Map();
            tree.forEach((n, id) => {
                for (const [, fpIdx] of (n.files || [])) {
                    if (!fpToNodes.has(fpIdx)) fpToNodes.set(fpIdx, new Set());
                    fpToNodes.get(fpIdx).add(id);
                }
            });
            function pointInQuad(ra, dec, quad) {
                let inside = false;
                for (let i = 0, j = 3; i < 4; j = i++) {
                    const xi = quad[i][0], yi = quad[i][1], xj = quad[j][0], yj = quad[j][1];
                    if (((yi > dec) !== (yj > dec)) &&
                        (ra < (xj - xi) * (dec - yi) / (yj - yi) + xi)) inside = !inside;
                }
                return inside;
            }
            function revealFoldersAt(ra, dec) {
                const nodeSet = new Set();
                for (let idx = 0; idx < FOLDERS.footprints.length; idx++) {
                    if (pointInQuad(ra, dec, FOLDERS.footprints[idx])) {
                        const ns = fpToNodes.get(idx);
                        if (ns) for (const id of ns) nodeSet.add(id);
                    }
                }
                if (!nodeSet.size) return;
                setFoldersCollapsed(false);
                for (const id of nodeSet) {           // expand each match's ancestors
                    let p = tree[id].parent;
                    while (p !== null && p !== undefined) { expanded.add(p); p = tree[p].parent; }
                }
                render();
                let first = null;
                for (const id of nodeSet) {
                    const row = treeDiv.querySelector('.folder-row[data-node-id="' + id + '"]');
                    if (!row) continue;
                    row.classList.add('flash');
                    setTimeout(() => row.classList.remove('flash'), 2000);
                    if (!first) first = row;
                }
                if (first) {
                    const r = first.getBoundingClientRect(), t = treeDiv.getBoundingClientRect();
                    treeDiv.scrollTop += (r.top - t.top) - t.height / 2;
                }
                toastMsg('Imaged in ' + nodeSet.size + ' folder' + (nodeSet.size > 1 ? 's' : ''));
            }

            // ---------- rescan buttons (one per rig + All) ----------
            const rescanBar = document.getElementById('rescan-bar');
            function buildRescanBar() {
                rescanBar.innerHTML = '<span>Rescan:</span>';
                const mk = (label, rig) => {
                    const b = document.createElement('button');
                    b.className = 'rescan-btn';
                    b.textContent = label;
                    b.addEventListener('click', () => doRescan(rig, b));
                    rescanBar.appendChild(b);
                };
                for (const r of DATA.rigs) mk(r.label, r.label);
                if (DATA.rigs.length > 1) mk('All', 'all');
            }
            async function doRescan(rig, btn) {
                const btns = rescanBar.querySelectorAll('.rescan-btn');
                btns.forEach(b => b.disabled = true);
                const label = btn.textContent;
                btn.textContent = 'Rescanning…';
                toastMsg('Rescanning ' + rig + ' — checking disk, this can take a moment…');
                try {
                    const resp = await fetch('/api/rescan?rig=' + encodeURIComponent(rig));
                    const r = await resp.json();
                    if (!r.ok) throw new Error(r.error || 'rescan failed');
                    toastMsg('Rescanned ' + r.rigs + ': ' + r.processed + ' new/changed, '
                             + r.deleted + ' removed → reloading');
                    setTimeout(() => window.location.reload(), 900);
                } catch (e) {
                    btns.forEach(b => b.disabled = false);
                    btn.textContent = label;
                    toastMsg('Rescan failed: ' + e.message);
                }
            }
            buildRescanBar();

            colorInput.addEventListener('input', () => { refreshHighlights(); saveState(); });
            document.getElementById('folders-clear').addEventListener('click', () => {
                selected.clear(); selectedFiles.clear(); refreshHighlights(); render();
            });
            document.getElementById('tree-collapse').addEventListener('click', () => {
                expanded.clear(); render();           // collapse the whole tree to rig roots
                treeDiv.scrollTop = 0;
            });

            // ---------- header popup: read a file's full FITS/XISF header live ----------
            const headerPanel = document.getElementById('header-panel');
            const headerTitle = document.getElementById('header-panel-title');
            const headerPre = document.getElementById('header-panel-pre');
            let headerReqId = 0;
            async function showHeader(path, fname) {
                const reqId = ++headerReqId;
                headerTitle.textContent = fname;
                headerPre.textContent = 'Loading header…';
                headerPanel.style.display = 'flex';
                try {
                    const resp = await fetch('/api/header?path=' + encodeURIComponent(path));
                    const text = await resp.text();
                    if (reqId !== headerReqId) return;   // a newer click superseded this
                    headerPre.textContent = resp.ok ? text : ('Could not read header (' + resp.status + ').\n' +
                        'Header reading needs the script-launched server on port ' + location.port + '.');
                } catch (e) {
                    if (reqId !== headerReqId) return;
                    headerPre.textContent = 'Could not reach the header service.\n' +
                        'Make sure the map was opened via the server the script started (python sky_mapper12.py).';
                }
            }
            document.getElementById('header-panel-close').addEventListener('click', () => {
                headerPanel.style.display = 'none';
            });
            // drag the header popup by its title bar
            const headerBar = document.getElementById('header-panel-bar');
            let hdrag = null;
            headerBar.addEventListener('pointerdown', (e) => {
                if (e.target.id === 'header-panel-close') return;
                hdrag = { x0: e.clientX, y0: e.clientY, left: headerPanel.offsetLeft, top: headerPanel.offsetTop };
                headerPanel.style.transform = 'none';
                headerPanel.style.left = hdrag.left + 'px';
                headerPanel.style.top = hdrag.top + 'px';
                try { headerBar.setPointerCapture(e.pointerId); } catch (err) {}
                e.preventDefault();
            });
            headerBar.addEventListener('pointermove', (e) => {
                if (!hdrag) return;
                headerPanel.style.left = (hdrag.left + e.clientX - hdrag.x0) + 'px';
                headerPanel.style.top = Math.max(0, hdrag.top + e.clientY - hdrag.y0) + 'px';
            });
            headerBar.addEventListener('pointerup', () => { hdrag = null; });

            // folder panel: draggable + collapsible (same behavior as the main panel)
            function setFoldersCollapsed(c) {
                foldersBody.style.display = c ? 'none' : 'flex';
                foldersCollapseBtn.innerText = c ? '▸' : '▾';
            }
            function placeFolders(x, y) {
                x = Math.max(0, Math.min(x, window.innerWidth - 80));
                y = Math.max(0, Math.min(y, window.innerHeight - 40));
                foldersPanel.style.left = x + 'px';
                foldersPanel.style.top = y + 'px';
                foldersPanel.style.maxHeight = (window.innerHeight - y - 20) + 'px';
            }
            placeFolders(SETTINGS.folderPanelX, SETTINGS.folderPanelY);
            let fdrag = null;
            foldersHeader.addEventListener('pointerdown', (e) => {
                fdrag = { x0: e.clientX, y0: e.clientY,
                          left: foldersPanel.offsetLeft, top: foldersPanel.offsetTop, moved: false };
                try { foldersHeader.setPointerCapture(e.pointerId); } catch (err) {}
                e.preventDefault();
            });
            foldersHeader.addEventListener('pointermove', (e) => {
                if (!fdrag) return;
                const dx = e.clientX - fdrag.x0, dy = e.clientY - fdrag.y0;
                if (Math.abs(dx) + Math.abs(dy) > 4) fdrag.moved = true;
                if (fdrag.moved) placeFolders(fdrag.left + dx, fdrag.top + dy);
            });
            foldersHeader.addEventListener('pointerup', () => {
                if (!fdrag) return;
                const moved = fdrag.moved;
                fdrag = null;
                if (!moved) setFoldersCollapsed(foldersBody.style.display !== 'none');
                saveState();
            });
            render();

            // console hook for debugging the folder browser
            window.folderBrowser = { tree, kids, selected, selectedFiles, expanded,
                                     highlightOverlay, refreshHighlights, showHeader,
                                     revealFoldersAt };

            // ================= menu + panels =================
            const PANELS = {
                controls: { el: document.getElementById('controls'), show: 'block' },
                folders:  { el: document.getElementById('folders'),  show: 'flex' },
                stats:    { el: document.getElementById('stats'),     show: 'flex' },
                catalog:  { el: document.getElementById('catalog'),   show: 'flex' },
                plan:     { el: document.getElementById('plan'),      show: 'flex' },
            };
            const openPanels = new Set(['controls']);
            let planMode = false;
            const panelPos = {};

            function applyPanels() {
                for (const [name, p] of Object.entries(PANELS)) {
                    p.el.style.display = openPanels.has(name) ? p.show : 'none';
                    const btn = document.querySelector('.menu-btn[data-panel="' + name + '"]');
                    if (btn) btn.classList.toggle('active', openPanels.has(name));
                }
            }
            const DEFAULT_POS = {
                stats:   () => [Math.max(360, window.innerWidth - 340), 46],
                catalog: () => [Math.max(360, window.innerWidth - 340), 46],
                plan:    () => [Math.round(window.innerWidth / 2 - 160), 80],
            };
            function placeGeneric(name) {
                const el = PANELS[name].el;
                const p = panelPos[name] || (DEFAULT_POS[name] ? DEFAULT_POS[name]() : [380, 60]);
                el.style.left = Math.max(0, Math.min(p[0], window.innerWidth - 80)) + 'px';
                el.style.top = Math.max(38, Math.min(p[1], window.innerHeight - 40)) + 'px';
            }
            function openPanel(name, on) {
                if (on === undefined) on = !openPanels.has(name);
                if (on) { openPanels.add(name); placeGeneric(name); } else openPanels.delete(name);
                if (on && name === 'stats') renderStats();
                if (on && name === 'catalog') ensureCatalog();
                applyPanels(); saveState();
            }
            document.querySelectorAll('.menu-btn[data-panel]').forEach(btn =>
                btn.addEventListener('click', () => openPanel(btn.dataset.panel)));
            document.querySelectorAll('.panel-x[data-panel]').forEach(btn =>
                btn.addEventListener('click', () => openPanel(btn.dataset.panel, false)));

            // generic drag for the stats/catalog/plan panels
            document.querySelectorAll('.panel-head[data-drag]').forEach(head => {
                const name = head.dataset.drag, el = PANELS[name].el;
                let dg = null;
                head.addEventListener('pointerdown', (e) => {
                    if (e.target.classList.contains('panel-x')) return;
                    dg = { x0: e.clientX, y0: e.clientY, left: el.offsetLeft, top: el.offsetTop };
                    try { head.setPointerCapture(e.pointerId); } catch (err) {}
                    e.preventDefault();
                });
                head.addEventListener('pointermove', (e) => {
                    if (!dg) return;
                    el.style.left = Math.max(0, dg.left + e.clientX - dg.x0) + 'px';
                    el.style.top = Math.max(38, dg.top + e.clientY - dg.y0) + 'px';
                });
                head.addEventListener('pointerup', () => {
                    if (!dg) return;
                    panelPos[name] = [el.offsetLeft, el.offsetTop];
                    dg = null; saveState();
                });
            });

            const planBtn = document.getElementById('plan-mode-btn');
            planBtn.addEventListener('click', () => {
                planMode = !planMode;
                planBtn.classList.toggle('active', planMode);
                toastMsg(planMode ? 'Plan mode ON — click anywhere on the sky' : 'Plan mode off');
                saveState();
            });
            document.getElementById('poster-btn').addEventListener('click', exportPoster);

            // ================= statistics =================
            function renderStats() {
                const S = STATS, h = (s) => fmtExposure(s);
                const filterRows = Object.entries(S.by_filter).map(([f, s]) =>
                    '<tr><td>' + esc(f) + '</td><td class="r">' + h(s) + '</td></tr>').join('');
                const rigRows = Object.entries(S.by_rig).map(([r, s]) =>
                    '<tr><td>' + esc(r) + '</td><td class="r">' + h(s) + '</td></tr>').join('');
                const months = Object.entries(S.by_month);
                const mMax = Math.max(1, ...months.map(m => m[1]));
                const monthRows = months.map(([m, s]) =>
                    '<tr><td>' + m + '</td><td><div class="barwrap"><span class="bar" style="width:' +
                    Math.round(s / mMax * 120) + 'px"></span><span style="color:#9bd">' + h(s) + '</span></div></td></tr>').join('');
                const topRows = S.top_targets.map(t =>
                    '<tr><td>' + esc(t.name) + '</td><td class="r">' + h(t.seconds) + ' · ' + t.frames + 'f</td></tr>').join('');
                const cond = S.cond_count
                    ? '<div class="stat-sec">Capture conditions (median)</div><table class="kv">' +
                      (S.median_alt != null ? '<tr><td>Altitude</td><td class="r">' + S.median_alt + '°</td></tr>' : '') +
                      (S.median_airmass != null ? '<tr><td>Airmass</td><td class="r">' + S.median_airmass + '</td></tr>' : '') +
                      (S.median_foctemp != null ? '<tr><td>Focuser temp</td><td class="r">' + S.median_foctemp + '°C</td></tr>' : '') +
                      '</table>' : '';
                document.getElementById('stats-body').innerHTML =
                    '<div class="stat-big">' +
                    '<div><div class="num">' + h(S.total_seconds) + '</div><div class="lbl">Integration</div></div>' +
                    '<div><div class="num">' + S.total_frames.toLocaleString() + '</div><div class="lbl">Light frames</div></div>' +
                    '<div><div class="num">' + S.targets + '</div><div class="lbl">Targets</div></div>' +
                    '<div><div class="num">' + S.nights + '</div><div class="lbl">Nights</div></div>' +
                    '</div>' +
                    '<div class="muted">' + (S.date_min || '?') + ' → ' + (S.date_max || '?') + '</div>' +
                    '<div class="stat-sec">Integration by filter</div><table class="kv">' + filterRows + '</table>' +
                    '<div class="stat-sec">By rig</div><table class="kv">' + rigRows + '</table>' +
                    '<div class="stat-sec">By month</div><table class="kv">' + monthRows + '</table>' +
                    '<div class="stat-sec">Most-imaged targets</div><table class="kv">' + topRows + '</table>' +
                    cond;
            }

            // ================= Messier / Caldwell catalog =================
            let catalogData = null, catLayers = null;
            function clusterCoversPoint(ra, dec) {
                for (const c of DATA.clusters) {
                    const r = c.rect;
                    if (dec < r[2] || dec > r[3]) continue;
                    if (pointInQuad(ra, dec, c.corners)) return true;
                }
                return false;
            }
            async function ensureCatalog() {
                if (catalogData) { applyCatalogVis(); return; }
                const body = document.getElementById('catalog-body');
                if (!SETTINGS.catalogPresent) { body.textContent = 'Catalog file missing — run generate_catalog.py.'; return; }
                body.textContent = 'Loading catalog…';
                try {
                    const r = await fetch(SETTINGS.catalogFile);
                    catalogData = await r.json();
                } catch (e) { body.textContent = 'Could not load catalog.'; return; }
                for (const o of catalogData) o.done = clusterCoversPoint(o.ra, o.dec);
                buildCatalogOverlays();
                renderCatalogPanel();
                applyCatalogVis();
            }
            function buildCatalogOverlays() {
                const mk = (color) => { const c = A.catalog({ name: 'cat', sourceSize: 9, shape: 'circle',
                    color: color, displayLabel: true, labelColumn: 'name', labelColor: color, labelFont: '10px sans-serif' });
                    aladin.addCatalog(c); return c; };
                catLayers = { M: { done: mk('#5d5'), todo: mk('#888') }, C: { done: mk('#5d5'), todo: mk('#888') } };
                for (const o of catalogData)
                    catLayers[o.cat][o.done ? 'done' : 'todo'].addSources([A.source(o.ra, o.dec, { name: o.id })]);
            }
            function applyCatalogVis() {
                if (!catLayers) return;
                const showM = document.getElementById('cat-m').checked;
                const showC = document.getElementById('cat-c').checked;
                const todoOnly = document.getElementById('cat-todo').checked;
                const set = (layer, on) => { try { on ? layer.show() : layer.hide(); } catch (e) {} };
                set(catLayers.M.done, showM && !todoOnly);
                set(catLayers.M.todo, showM);
                set(catLayers.C.done, showC && !todoOnly);
                set(catLayers.C.todo, showC);
            }
            function renderCatalogPanel() {
                const mDone = catalogData.filter(o => o.cat === 'M' && o.done).length;
                const mAll = catalogData.filter(o => o.cat === 'M').length;
                const cDone = catalogData.filter(o => o.cat === 'C' && o.done).length;
                const cAll = catalogData.filter(o => o.cat === 'C').length;
                document.getElementById('catalog-body').innerHTML =
                    '<div class="cat-toggles">' +
                    '<label class="checkbox-container"><input type="checkbox" id="cat-m" checked> Messier</label>' +
                    '<label class="checkbox-container"><input type="checkbox" id="cat-c" checked> Caldwell</label>' +
                    '<label class="checkbox-container"><input type="checkbox" id="cat-todo"> Only to-do</label>' +
                    '</div>' +
                    '<div class="cat-prog">Messier <span class="pct">' + mDone + '/' + mAll + '</span> &nbsp; ' +
                    'Caldwell <span class="pct">' + cDone + '/' + cAll + '</span></div>' +
                    '<div id="catalog-list"></div>';
                ['cat-m', 'cat-c', 'cat-todo'].forEach(id => document.getElementById(id)
                    .addEventListener('change', () => { applyCatalogVis(); renderCatalogList(); }));
                renderCatalogList();
            }
            function renderCatalogList() {
                const showM = document.getElementById('cat-m').checked;
                const showC = document.getElementById('cat-c').checked;
                const todoOnly = document.getElementById('cat-todo').checked;
                const list = document.getElementById('catalog-list');
                const frag = document.createDocumentFragment();
                for (const o of catalogData) {
                    if (o.cat === 'M' && !showM) continue;
                    if (o.cat === 'C' && !showC) continue;
                    if (todoOnly && o.done) continue;
                    const row = document.createElement('div');
                    row.className = 'cat-row ' + (o.done ? 'done' : 'todo');
                    row.innerHTML = '<span class="ck">' + (o.done ? '✓' : '○') + '</span>' +
                        '<span>' + esc(o.id) + (o.common ? ' · ' + esc(o.common) : ' · ' + esc(o.desig)) +
                        ' <span style="color:#777">(' + esc(o.otype) + ')</span></span>';
                    row.addEventListener('click', () => { aladin.gotoRaDec(o.ra, o.dec); aladin.setFov(2.5); });
                    frag.appendChild(row);
                }
                list.innerHTML = '';
                list.appendChild(frag);
            }

            // ================= altitude curve (shared by plan + details) =================
            function altitudeCurve(canvas, ra, dec) {
                if (SETTINGS.lat === null) return null;
                const ctx = canvas.getContext('2d');
                const W = canvas.width, H = canvas.height;
                ctx.clearRect(0, 0, W, H);
                const start = new Date(); start.setHours(16, 0, 0, 0);
                const span = 18 * 3600 * 1000, N = 96;
                const X = i => 4 + i * (W - 8) / N;
                const Y = a => H - 14 - (Math.max(-10, Math.min(90, a)) + 10) / 100 * (H - 18);
                let maxAlt = -90, maxI = 0; const alts = [], dark = [];
                for (let i = 0; i <= N; i++) {
                    const jd = julianDate(new Date(start.getTime() + i * span / N));
                    const a = altDeg(ra, dec, jd); alts.push(a);
                    const s = sunRaDec(jd); dark.push(altDeg(s.ra, s.dec, jd) < -18);
                    if (a > maxAlt) { maxAlt = a; maxI = i; }
                }
                ctx.fillStyle = 'rgba(60,90,150,0.30)';
                for (let i = 0; i < N; i++) if (dark[i]) ctx.fillRect(X(i), 0, (W - 8) / N + 1, H - 12);
                ctx.strokeStyle = '#333'; ctx.beginPath(); ctx.moveTo(4, Y(0)); ctx.lineTo(W - 4, Y(0)); ctx.stroke();
                ctx.strokeStyle = '#664'; ctx.setLineDash([3, 3]); ctx.beginPath();
                ctx.moveTo(4, Y(SETTINGS.minAlt)); ctx.lineTo(W - 4, Y(SETTINGS.minAlt)); ctx.stroke(); ctx.setLineDash([]);
                ctx.strokeStyle = '#6cf'; ctx.lineWidth = 2; ctx.beginPath();
                alts.forEach((a, i) => { const px = X(i), py = Y(a); i ? ctx.lineTo(px, py) : ctx.moveTo(px, py); }); ctx.stroke();
                ctx.fillStyle = '#fd6'; ctx.beginPath(); ctx.arc(X(maxI), Y(maxAlt), 3, 0, 7); ctx.fill();
                ctx.fillStyle = '#888'; ctx.font = '9px sans-serif';
                ctx.fillText('90°', 2, 9); ctx.fillText(SETTINGS.minAlt + '°', 2, Y(SETTINGS.minAlt) - 2);
                [['18h', 2 / 18], ['00h', 8 / 18], ['06h', 14 / 18]].forEach(([t, f]) =>
                    ctx.fillText(t, X(N * f) - 6, H - 2));
                return { maxAlt: Math.round(maxAlt), maxTime: new Date(start.getTime() + maxI * span / N) };
            }
            const hhmm = d => String(d.getHours()).padStart(2, '0') + ':' + String(d.getMinutes()).padStart(2, '0');

            // ================= target file export =================
            function downloadText(filename, text, mime) {
                const blob = new Blob([text], { type: mime || 'text/plain' });
                const a = document.createElement('a');
                a.href = URL.createObjectURL(blob); a.download = filename; a.click();
                URL.revokeObjectURL(a.href);
            }
            function ninaTargetCSV(name, ra, dec) {
                // Telescopius-style target list — importable via NINA's Telescopius CSV import
                const raS = raToHms(ra).replace(/\s+/g, '');
                const decS = decToDms(dec).replace(/\s+/g, '');
                const esccsv = v => '"' + String(v).replace(/"/g, '""') + '"';
                return 'Familiar Name,Catalogue Entry,Right Ascension,Declination\n' +
                    [name, name, raS, decS].map(esccsv).join(',') + '\n';
            }

            // ================= plan popup (click empty sky in plan mode) =================
            function openPlanPopup(ra, dec) {
                openPanels.add('plan'); placeGeneric('plan'); applyPanels();
                const body = document.getElementById('plan-body');
                const safe = (raToHms(ra) + ' ' + decToDms(dec));
                body.innerHTML =
                    '<div class="muted">' + esc(safe) + '</div>' +
                    '<div class="muted">' + ra.toFixed(5) + ', ' + dec.toFixed(5) + '</div>' +
                    '<canvas id="plan-alt"></canvas>' +
                    '<div id="plan-alt-info" class="muted"></div>' +
                    '<div id="plan-where" class="muted">…</div>' +
                    '<button class="search-btn" id="plan-copy" style="width:100%;margin-top:8px;">Copy coordinates</button>' +
                    '<button class="search-btn" id="plan-nina" style="width:100%;margin-top:6px;">Download NINA target (CSV)</button>';
                const cv = document.getElementById('plan-alt');
                cv.width = cv.clientWidth || 296; cv.height = 120;
                const r = altitudeCurve(cv, ra, dec);
                document.getElementById('plan-alt-info').innerText = r
                    ? 'Tonight: max altitude ' + r.maxAlt + '° at ' + hhmm(r.maxTime)
                    : 'Set OBSERVER_LAT/LON for altitude.';
                document.getElementById('plan-copy').addEventListener('click', () => {
                    navigator.clipboard.writeText(ra.toFixed(6) + ', ' + dec.toFixed(6))
                        .then(() => toastMsg('Copied coordinates')).catch(() => {});
                });
                document.getElementById('plan-nina').addEventListener('click', () => {
                    const nm = 'SkyMap ' + raToHms(ra).replace(/\s+/g, '') + decToDms(dec).replace(/\s+/g, '');
                    downloadText('nina_target.csv', ninaTargetCSV(nm, ra, dec), 'text/csv');
                    toastMsg('NINA target CSV downloaded');
                });
                fetch('/api/objectinfo?ra=' + ra + '&dec=' + dec).then(r => r.json()).then(info => {
                    const w = document.getElementById('plan-where');
                    if (!w) return;
                    let s = info.constellation ? 'Constellation: ' + esc(info.constellation) : '';
                    if (info.simbad_url) s += (s ? ' · ' : '') + '<a href="' + info.simbad_url +
                        '" target="_blank" rel="noopener" style="color:#6cf">SIMBAD here</a>';
                    if (info.ned_url) s += ' · <a href="' + info.ned_url +
                        '" target="_blank" rel="noopener" style="color:#6cf">NED here</a>';
                    w.innerHTML = s;
                }).catch(() => {});
            }

            // ================= poster export =================
            function exportPoster() {
                try {
                    const r = aladin.getViewDataURL ? aladin.getViewDataURL({ format: 'image/png' }) : null;
                    if (r) {
                        Promise.resolve(r).then(url => { downloadHref(url, 'sky_coverage.png'); toastMsg('Poster exported'); });
                        return;
                    }
                } catch (e) {}
                const cv = document.querySelector('#aladin-lite-div canvas');
                if (cv) { downloadHref(cv.toDataURL('image/png'), 'sky_coverage.png'); toastMsg('Poster exported'); }
                else toastMsg('Export not supported by this Aladin build');
            }
            function downloadHref(href, name) {
                const a = document.createElement('a'); a.href = href; a.download = name; a.click();
            }

            // ---------- remember UI state between sessions ----------
            function saveState() {
                try {
                    localStorage.setItem(LS_KEY, JSON.stringify({
                        survey: els.survey.value,
                        grid: els.grid.checked,
                        constellations: els.constellations.checked,
                        mosaics: els.mosaics.checked,
                        tonight: els.tonight.checked,
                        flipH: els.flipH.checked,
                        flipV: els.flipV.checked,
                        rigs: rigEnabled,
                        collapsed: controlsBody.style.display === 'none',
                        panelX: controls.offsetLeft,
                        panelY: controls.offsetTop,
                        highlightColor: colorInput.value,
                        folderCollapsed: foldersBody.style.display === 'none',
                        folderPanelX: foldersPanel.offsetLeft,
                        folderPanelY: foldersPanel.offsetTop,
                        openPanels: Array.from(openPanels).filter(n => n !== 'plan'),
                        planMode: planMode,
                        panelPos: panelPos
                    }));
                } catch (e) {}
            }
            function restoreState() {
                let s = null;
                try { s = JSON.parse(localStorage.getItem(LS_KEY) || 'null'); } catch (e) {}
                if (!s) return;
                if (s.survey && [...els.survey.options].some(o => o.value === s.survey)) {
                    els.survey.value = s.survey;
                    aladin.setBaseImageLayer(aladin.createImageSurvey(s.survey));
                }
                for (const key of ['grid', 'constellations', 'mosaics', 'tonight', 'flipH', 'flipV']) {
                    if (typeof s[key] === 'boolean') els[key].checked = s[key];
                }
                if (els.flipH.checked || els.flipV.checked) applyFlips();
                if (s.rigs) {
                    for (const color in s.rigs) {
                        if (rigCheckboxes[color]) {
                            rigCheckboxes[color].checked = !!s.rigs[color];
                            rigEnabled[color] = !!s.rigs[color];
                        }
                    }
                }
                if (typeof s.panelX === 'number' && typeof s.panelY === 'number') {
                    placePanel(s.panelX, s.panelY);
                }
                if (s.collapsed) setCollapsed(true);
                if (typeof s.highlightColor === 'string') {
                    colorInput.value = s.highlightColor;
                    refreshHighlights();
                }
                if (typeof s.folderPanelX === 'number' && typeof s.folderPanelY === 'number') {
                    placeFolders(s.folderPanelX, s.folderPanelY);
                }
                if (s.folderCollapsed) setFoldersCollapsed(true);
                if (els.grid.checked) aladin.showCooGrid();
                if (els.constellations.checked) ensureConstellations();
                if (els.tonight.checked && SETTINGS.lat === null) els.tonight.checked = false;
                if (els.tonight.checked) ensureTonightOutline();
                if (s.panelPos) Object.assign(panelPos, s.panelPos);
                if (Array.isArray(s.openPanels)) {
                    openPanels.clear();
                    s.openPanels.forEach(n => { if (PANELS[n] && n !== 'plan') openPanels.add(n); });
                }
                if (typeof s.planMode === 'boolean') {
                    planMode = s.planMode;
                    planBtn.classList.toggle('active', planMode);
                }
            }

            restoreState();
            applyPanels();
            if (openPanels.has('stats')) renderStats();
            if (openPanels.has('catalog')) ensureCatalog();
            for (const n of openPanels) if (DEFAULT_POS[n]) placeGeneric(n);
            rebuild();
            updateLabelVisibility();
            console.log("✅ Sky map ready: " + DATA.clusters.length + " clusters");
        }).catch(err => console.error("Aladin error:", err));
    </script>
</body>
</html>
"""


def resolve_observer_site(entries: list) -> tuple[float | None, float | None]:
    """Configured OBSERVER_LAT/LON wins; otherwise the most common header site."""
    if OBSERVER_LAT is not None and OBSERVER_LON is not None:
        return OBSERVER_LAT, OBSERVER_LON
    sites = Counter(tuple(e['site']) for e in entries if e.get('site'))
    if sites:
        lat, lon = sites.most_common(1)[0][0]
        print(f"\n📍 Observing site from image headers: lat {lat}, lon {lon}")
        return lat, lon
    return None, None


def generate_aladin_html(payload: list, folder_data: dict, stats: dict, site_lat, site_lon) -> None:
    initial_ra, initial_dec, initial_fov = 0.0, 0.0, 180
    if payload:
        vsum = [0.0, 0.0, 0.0]
        for c in payload:
            v = _radec_to_vector(c['ra'], c['dec'])
            vsum[0] += v[0]; vsum[1] += v[1]; vsum[2] += v[2]
        initial_ra, initial_dec = _vector_to_radec(vsum)
        initial_fov = 180 if len(payload) > 1 else 60

    data = {
        'clusters': payload,
        'rigs': [{'color': color, 'label': label} for _, color, label in SEARCH_DIRS],
    }
    settings = {
        'initialRa': round(initial_ra, 6),
        'initialDec': round(initial_dec, 6),
        'initialFov': initial_fov,
        'lat': site_lat,
        'lon': site_lon,
        'minAlt': MIN_ALTITUDE_DEG,
        'labelHideFov': LABEL_HIDE_FOV_DEG,
        'panelX': PANEL_X,
        'panelY': PANEL_Y,
        'folderPanelX': FOLDER_PANEL_X,
        'folderPanelY': FOLDER_PANEL_Y,
        'highlightColor': HIGHLIGHT_COLOR,
        'constellationLocal': os.path.basename(CONSTELLATION_FILE),
        'constellationUrl': CONSTELLATION_URL,
        'catalogFile': os.path.basename(CATALOG_FILE),
        'catalogPresent': os.path.exists(CATALOG_FILE),
    }

    html = (HTML_TEMPLATE
            .replace('__DATA_JSON__', json.dumps(data).replace('</', '<\\/'))
            .replace('__FOLDERS_JSON__', json.dumps(folder_data).replace('</', '<\\/'))
            .replace('__STATS_JSON__', json.dumps(stats).replace('</', '<\\/'))
            .replace('__SETTINGS_JSON__', json.dumps(settings)))

    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"\n✅ Map generated → {OUTPUT_HTML}")


def ensure_constellation_data() -> None:
    if os.path.exists(CONSTELLATION_FILE):
        return
    try:
        print("⬇ Downloading constellation line data (one-time)…")
        with urllib.request.urlopen(CONSTELLATION_URL, timeout=15) as resp, \
                open(CONSTELLATION_FILE, 'wb') as f:
            f.write(resp.read())
        print(f"   Saved → {CONSTELLATION_FILE}")
    except Exception as e:
        print(f"  ⚠ Could not download constellation data ({e}).")
        print("    The map will fetch it from the CDN instead when the checkbox is enabled.")


# =============================================================================
# OBJECT REFERENCE DATA (SIMBAD)
# =============================================================================

# SIMBAD object-type codes → friendly labels for the common amateur targets.
OTYPE_NAMES = {
    'G': 'Galaxy', 'GiG': 'Galaxy in Group', 'GiC': 'Galaxy in Cluster', 'GiP': 'Galaxy in Pair',
    'IG': 'Interacting Galaxies', 'PaG': 'Pair of Galaxies', 'GrG': 'Group of Galaxies',
    'ClG': 'Cluster of Galaxies', 'SCG': 'Supercluster of Galaxies',
    'AGN': 'Active Galaxy Nucleus', 'Sy1': 'Seyfert 1 Galaxy', 'Sy2': 'Seyfert 2 Galaxy',
    'SyG': 'Seyfert Galaxy', 'QSO': 'Quasar', 'rG': 'Radio Galaxy', 'LIN': 'LINER Galaxy',
    'SBG': 'Starburst Galaxy', 'H2G': 'HII Galaxy', 'EmG': 'Emission-line Galaxy', 'LSB': 'Low Surface Brightness Galaxy',
    'GlC': 'Globular Cluster', 'OpC': 'Open Cluster', 'Cl*': 'Star Cluster', 'As*': 'Stellar Association',
    'PN': 'Planetary Nebula', 'SNR': 'Supernova Remnant', 'HII': 'HII / Emission Nebula',
    'RNe': 'Reflection Nebula', 'DNe': 'Dark Nebula', 'GNe': 'Galactic Nebula', 'Neb': 'Nebula',
    'EmO': 'Emission Object', 'Cld': 'Cloud', 'MoC': 'Molecular Cloud', 'ISM': 'Interstellar Matter',
    'Y*O': 'Young Stellar Object', 'out': 'Outflow', 'HH': 'Herbig-Haro Object',
    '*': 'Star', '**': 'Double / Multiple Star', 'V*': 'Variable Star', 'C*': 'Carbon Star',
    'Em*': 'Emission-line Star', 'WR*': 'Wolf-Rayet Star', 'Be*': 'Be Star', 'RG*': 'Red Giant',
    'No*': 'Nova', 'SN*': 'Supernova', 'Ce*': 'Cepheid', 'RR*': 'RR Lyrae', 'Mi*': 'Mira Variable',
    'Sy*': 'Symbiotic Star', 'XB': 'X-ray Binary', 'CV*': 'Cataclysmic Variable',
}

_objinfo_cache = None
_OBJINFO_LOCK = threading.Lock()


def query_simbad(name: str):
    """Fetch and parse SIMBAD's ASCII record for `name`. Returns a dict of fields,
    or None if SIMBAD doesn't resolve it. (No astroquery dependency.)"""
    url = ("https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident="
           + quote(name))
    try:
        with urllib.request.urlopen(url, timeout=8) as r:
            txt = r.read().decode('utf-8', 'replace')
    except Exception:
        return None

    d = {}
    m = re.search(r'^Object\s+(.+?)\s+---\s+(\S+)', txt, re.M)
    if not m:
        return None                      # not resolved
    d['main_id'] = m.group(1).strip()
    d['otype_code'] = m.group(2).strip()
    d['otype'] = OTYPE_NAMES.get(d['otype_code'], d['otype_code'])

    m = re.search(r'Morphological type:\s*(\S+)', txt)
    if m and m.group(1) != '~':
        d['morph'] = m.group(1)
    m = re.search(r'Angular size:\s*([\d.]+)\s+([\d.]+)', txt)
    if m:
        d['size'] = f"{float(m.group(1)):.1f} × {float(m.group(2)):.1f}′"
    for band in ('V', 'B', 'J', 'H', 'K'):
        m = re.search(r'Flux ' + band + r'\s*:\s*([-\d.]+)', txt)
        if m:
            d.setdefault('mag', {})[band] = m.group(1)
    m = re.search(r'Redshift:\s*([-\d.]+)', txt)
    if m:
        d['z'] = m.group(1)
    m = re.search(r'Radial Velocity:\s*([-\d.]+)', txt)
    if m:
        d['rv'] = m.group(1)
    m = re.search(r'Distance:\s*([\d.]+)\s+(\w+)', txt)
    if m:
        d['distance'] = f"{m.group(1)} {m.group(2)}"
    names = re.findall(r'NAME\s+([A-Za-z][\w\'\.\- ]*?)(?:\s{2,}|\n)', txt)
    common = sorted({n.strip() for n in names if len(n.strip()) > 2})
    if common:
        d['common'] = common[:6]
    return d


def simbad_lookup(name: str):
    """SIMBAD lookup with a persistent on-disk cache (caches misses too)."""
    global _objinfo_cache
    if _objinfo_cache is None:
        try:
            with open(OBJECT_INFO_CACHE) as f:
                _objinfo_cache = json.load(f)
        except Exception:
            _objinfo_cache = {}
    key = name.strip().lower()
    if key in _objinfo_cache:
        return _objinfo_cache[key]
    data = query_simbad(name)
    _objinfo_cache[key] = data
    try:
        with open(OBJECT_INFO_CACHE, 'w') as f:
            json.dump(_objinfo_cache, f)
    except Exception:
        pass
    return data


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
            result['simbad_url'] = (f'https://simbad.cds.unistra.fr/simbad/sim-coo?Coord='
                                    f'{ra}+{dec}&Radius=5&Radius.unit=arcmin')
            result['ned_url'] = (f'https://ned.ipac.caltech.edu/conesearch?in_csys=Equatorial'
                                 f'&in_equinox=J2000&ra={ra}&dec={dec}&radius=5')
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


def ensure_web_server() -> None:
    url = f"http://localhost:{WEB_SERVER_PORT}/{os.path.basename(OUTPUT_HTML)}"

    port_in_use = True
    try:
        socket.create_connection(('127.0.0.1', WEB_SERVER_PORT), timeout=0.3).close()
    except OSError:
        port_in_use = False

    if port_in_use:
        if _server_is_ours(WEB_SERVER_PORT):
            print(f"\n🌐 Web server already running → {url}  (just refresh your browser)")
            return
        # Something else holds our port (e.g. an old plain http.server without the
        # header endpoint). Replace it so /api/header works.
        print(f"\n♻  Replacing existing server on port {WEB_SERVER_PORT} …")
        os.system(f"kill $(lsof -ti :{WEB_SERVER_PORT}) 2>/dev/null")
        time.sleep(0.6)

    subprocess.Popen(
        [sys.executable, os.path.abspath(__file__), '--serve', str(WEB_SERVER_PORT)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    print(f"\n🌐 Started web server → {url}")
    print(f"   (Keeps running in the background; stop it with:  kill $(lsof -ti :{WEB_SERVER_PORT}) )")
    webbrowser.open(url)


# =============================================================================
# MAIN
# =============================================================================

def scan_files(cache: dict, only_rigs=None) -> tuple[int, int, Counter]:
    """Scan rig folders, processing new/changed files and pruning cache entries
    for files deleted on disk. `only_rigs` limits the scan to those rig labels."""
    updated = 0
    deleted = 0
    skipped = Counter()

    for search_dir, color, rig_label in SEARCH_DIRS:
        if only_rigs is not None and rig_label not in only_rigs:
            continue
        print(f"\n📂 Scanning [{rig_label}] {search_dir} …")
        all_files = []
        for ext in ('*.fits', '*.fit', '*.xisf'):
            all_files.extend(glob.glob(os.path.join(search_dir, '**', ext), recursive=True))
        print(f"   Found {len(all_files)} file(s).")

        # Prune cache entries for files under this rig dir that no longer exist on
        # disk (only when the rig dir itself is present — don't prune a missing
        # drive). Match by path prefix so skipped entries are pruned too.
        if os.path.isdir(search_dir):
            existing = set(all_files)
            prefix = os.path.join(search_dir, '')   # ensure trailing separator
            gone = [fp for fp in cache if fp.startswith(prefix) and fp not in existing]
            for fp in gone:
                del cache[fp]
            deleted += len(gone)
            if gone:
                print(f"   Pruned {len(gone)} deleted file(s) from cache.")

        for filepath in all_files:
            try:
                mtime = os.path.getmtime(filepath)
            except OSError:
                continue
            cached = cache.get(filepath)

            if cached and cached.get('mtime') == mtime and cached.get('ver') == CACHE_SCHEMA:
                if 'skip' not in cached and cached.get('color') != color:
                    cached['color'] = color
                    cached['search_dir'] = search_dir
                    cached['rig'] = rig_label
                    updated += 1
                continue

            if filepath.lower().endswith(('.fits', '.fit')):
                data = get_fits_data(filepath)
            elif filepath.lower().endswith('.xisf'):
                data = get_xisf_data(filepath)
            else:
                continue

            entry = {'ver': CACHE_SCHEMA, 'mtime': mtime,
                     'search_dir': search_dir, 'rig': rig_label}
            entry.update(data)
            if 'skip' in data:
                skipped[data['skip']] += 1
            else:
                entry['color'] = color
                entry['mosaic'] = is_mosaic(filepath, data['target'])
                print(f"   Processed: {os.path.basename(filepath)} → {entry['target']}")
            cache[filepath] = entry
            updated += 1

    return updated, deleted, skipped


def collect_entries(cache: dict) -> tuple[list, Counter]:
    color_map = {sd: c for sd, c, _ in SEARCH_DIRS}
    rig_map = {sd: l for sd, _, l in SEARCH_DIRS}
    default_color = SEARCH_DIRS[0][1] if SEARCH_DIRS else '#ffffff'

    entries = []
    skip_counts = Counter()
    for filepath, entry in cache.items():
        if entry.get('skip'):
            skip_counts[entry['skip']] += 1
            continue
        if 'ra' not in entry:
            continue
        if should_exclude_target(entry.get('target', '')):
            skip_counts['excluded_keyword'] += 1
            continue
        sd = entry.get('search_dir', '')
        # Drop entries for files deleted since the last scan (but keep entries
        # whose whole search dir is missing — e.g. a disconnected drive).
        if sd and os.path.isdir(sd) and not os.path.exists(filepath):
            skip_counts['deleted'] += 1
            continue

        e = dict(entry)
        e['color'] = entry.get('color') or color_map.get(sd, default_color)
        e['rig'] = entry.get('rig') or rig_map.get(sd, '')
        entries.append(e)

    return entries, skip_counts


RESCAN_LOCK = threading.Lock()


def rebuild_catalog(only_rigs=None) -> dict:
    """Scan (optionally just some rigs), prune deleted files, recluster, and
    regenerate the map HTML. Returns a summary dict. Used by both the CLI and
    the /api/rescan endpoint, so it's serialized by RESCAN_LOCK there."""
    cache = load_cache()

    updated, deleted, skipped_run = scan_files(cache, only_rigs=only_rigs)
    if updated or deleted:
        save_cache(cache)
        print(f"\n💾 Cache updated: {updated} new/changed, {deleted} removed.")
        if skipped_run:
            print("   Newly skipped: " + ", ".join(f"{n} {reason}" for reason, n in skipped_run.items()))
    else:
        print("\nNo changes found. Using existing cache.")

    entries, skip_counts = collect_entries(cache)
    clusters = cluster_by_distance(entries, GROUPING_TOLERANCE_DEG)
    payload = build_cluster_payload(clusters)

    print(f"\n🔭 {len(entries)} light frames → {len(payload)} clustered targets.")
    if skip_counts:
        print("   Not shown: " + ", ".join(f"{n} {reason}" for reason, n in sorted(skip_counts.items())))
    mosaic_count = sum(1 for c in payload if c['mosaic'])
    if mosaic_count:
        print(f"   {mosaic_count} cluster(s) flagged as mosaic panels (toggle visibility in the map).")

    folder_data = build_folder_data(cache)
    print(f"   Folder browser: {len(folder_data['tree'])} folders, "
          f"{len(folder_data['footprints'])} unique footprints.")

    stats = build_stats(cache)

    site_lat, site_lon = resolve_observer_site(entries)
    ensure_constellation_data()
    generate_aladin_html(payload, folder_data, stats, site_lat, site_lon)

    return {
        'rigs': 'all rigs' if only_rigs is None else ', '.join(sorted(only_rigs)),
        'processed': updated, 'deleted': deleted,
        'targets': len(payload), 'lights': len(entries),
        'folders': len(folder_data['tree']),
    }


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
