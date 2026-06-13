# Copy this file to local_config.py and edit it with your own settings.
# local_config.py is git-ignored, so your real paths stay out of the repo.
# Anything you set here overrides the matching constant in sky_mapper12.py.

# One entry per imaging rig / image folder: (absolute_path, hex_color, label)
SEARCH_DIRS = [
    (r"/path/to/rig_a/images/", "#00ff00", "Rig A"),   # green
    (r"/path/to/rig_b/images/", "#ff4444", "Rig B"),   # red
    # (r"/Volumes/NAS/rig_c/", "#44aaff", "Rig C"),
]

# Optional — observing site for the altitude/visibility features.
# Leave commented out to auto-detect from your FITS/XISF headers (SITELAT/SITELONG).
# OBSERVER_LAT = 42.0    # degrees north
# OBSERVER_LON = -71.0   # degrees east (US longitudes are negative)
# MIN_ALTITUDE_DEG = 30.0

# Optional — default outline color for folder-browser highlights.
# HIGHLIGHT_COLOR = "#ff00ff"

# Optional — web server port (the map is served at http://localhost:<port>/).
# WEB_SERVER_PORT = 8001

# Optional — target name keywords to exclude from the map.
# EXCLUDE_KEYWORDS = ["flatwizard", "snapshot", "tsuchinshan"]
