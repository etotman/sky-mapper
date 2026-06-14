"""Generate the Aladin map HTML from template.html, plus site/constellation helpers."""
import os, json
import urllib.request
from collections import Counter
from config import *
from skymap_catalog import _radec_to_vector, _vector_to_radec

_TEMPLATE_FILE = os.path.join(SCRIPT_DIR, "template.html")

def _read_template():
    with open(_TEMPLATE_FILE, encoding="utf-8") as f:
        return f.read()

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

    html = (_read_template()
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

