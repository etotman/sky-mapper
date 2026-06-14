"""Configuration constants. Override machine-specific values in local_config.py
(git-ignored). See local_config.example.py."""
import os

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

# Folder names to skip during scanning (case-insensitive, exact directory-name
# match). Any file whose path contains one of these as a path component is ignored.
EXCLUDE_FOLDERS: set = set()

# Target-name substrings to exclude from the map (case-insensitive). Keep these
# distinctive so they don't accidentally match real deep-sky target names.
EXCLUDE_KEYWORDS = [
    "flatwizard",
    "snapshot",
]

# Solar-system bodies move against the fixed star field, so it makes no sense to
# plot frames whose target IS one of these. Matched as the EXACT target name
# (case-insensitive) — so deep-sky objects named after planets, e.g. the Saturn
# Nebula (NGC 7009) or the Ghost of Jupiter (NGC 3242), are NOT excluded.
EXCLUDE_TARGETS = {
    "moon", "sun", "mercury", "venus", "mars", "jupiter",
    "saturn", "uranus", "neptune", "pluto",
}

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
SERVER_VERSION  = 18           # bump when the server's API code changes, to force a
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
# Local overrides: your real rig folders / site / preferences live in
# local_config.py (git-ignored) and may redefine any constant above.
# ---------------------------------------------------------------------------
try:
    from local_config import *  # noqa: F401,F403
except ImportError:
    pass
