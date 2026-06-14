"""SIMBAD object reference-data lookup (cached on disk)."""
import re, json, threading
import urllib.request
from urllib.parse import quote
from config import *

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

