"""One-off: build the Messier, Caldwell, and GI750 catalogs and write
messier_caldwell.json (ra/dec/type/mag/size/common name/score) for the sky map
overlay. Messier and Caldwell are resolved via SIMBAD; GI750 is a curated, scored
target list shipped as gi750.csv. Run once; the JSON is committed so the map needs
no live lookups."""
import os, re, csv, json, sys, time, urllib.request, urllib.parse
from astropy.coordinates import SkyCoord
import astropy.units as u

# Keep redirected/piped output from crashing on non-ASCII object names on Windows,
# where stdout otherwise defaults to the legacy locale codepage (e.g. cp1252).
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'messier_caldwell.json')

# Drop objects too far south to be observable from the northern hemisphere.
MIN_DEC = -40.0

# Caldwell number -> designation SIMBAD understands (Caldwell isn't a SIMBAD catalog).
CALDWELL = {
 1:'NGC 188',2:'NGC 40',3:'NGC 4236',4:'NGC 7023',5:'IC 342',6:'NGC 6543',7:'NGC 2403',
 8:'NGC 559',9:'Sh2-155',10:'NGC 663',11:'NGC 7635',12:'NGC 6946',13:'NGC 457',
 14:'NGC 869',15:'NGC 6826',16:'NGC 7243',17:'NGC 147',18:'NGC 185',19:'IC 5146',
 20:'NGC 7000',21:'NGC 4449',22:'NGC 7662',23:'NGC 891',24:'NGC 1275',25:'NGC 2419',
 26:'NGC 4244',27:'NGC 6888',28:'NGC 752',29:'NGC 5005',30:'NGC 7331',31:'IC 405',
 32:'NGC 4631',33:'NGC 6992',34:'NGC 6960',35:'NGC 4889',36:'NGC 4559',37:'NGC 6885',
 38:'NGC 4565',39:'NGC 2392',40:'NGC 3626',41:'Melotte 25',42:'NGC 7006',43:'NGC 7814',
 44:'NGC 7479',45:'NGC 5248',46:'NGC 2261',47:'NGC 6934',48:'NGC 2775',49:'NGC 2237',
 50:'NGC 2244',51:'IC 1613',52:'NGC 4697',53:'NGC 3115',54:'NGC 2506',55:'NGC 7009',
 56:'NGC 246',57:'NGC 6822',58:'NGC 2360',59:'NGC 3242',60:'NGC 4038',61:'NGC 4039',
 62:'NGC 247',63:'NGC 7293',64:'NGC 2362',65:'NGC 253',66:'NGC 5694',67:'NGC 1097',
 68:'NGC 6729',69:'NGC 6302',70:'NGC 300',71:'NGC 2477',72:'NGC 55',73:'NGC 1851',
 74:'NGC 3132',75:'NGC 6124',76:'NGC 6231',77:'NGC 5128',78:'NGC 6541',79:'NGC 3201',
 80:'NGC 5139',81:'NGC 6352',82:'NGC 6193',83:'NGC 4945',84:'NGC 5286',85:'IC 2391',
 86:'NGC 6397',87:'NGC 1261',88:'NGC 5823',89:'NGC 6087',90:'NGC 2867',91:'NGC 3532',
 92:'NGC 3372',93:'NGC 6752',94:'NGC 4755',95:'NGC 6025',96:'NGC 2516',97:'NGC 3766',
 98:'NGC 4609',99:'Coalsack',100:'IC 2944',101:'NGC 6744',102:'IC 2602',103:'NGC 2070',
 104:'NGC 362',105:'NGC 4833',106:'NGC 104',107:'NGC 6101',108:'NGC 4372',109:'NGC 3195',
}

OTYPE = {  # condensed SIMBAD type -> short label
 'G':'Galaxy','GiG':'Galaxy','GiC':'Galaxy','IG':'Galaxies','AGN':'Active Galaxy','Sy1':'Galaxy',
 'Sy2':'Galaxy','SyG':'Galaxy','rG':'Galaxy','LIN':'Galaxy','SBG':'Galaxy','EmG':'Galaxy','H2G':'Galaxy',
 'GlC':'Globular Cluster','OpC':'Open Cluster','Cl*':'Star Cluster','As*':'Association',
 'PN':'Planetary Nebula','SNR':'Supernova Remnant','HII':'Emission Nebula','RNe':'Reflection Nebula',
 'DNe':'Dark Nebula','GNe':'Nebula','Neb':'Nebula','EmO':'Nebula','*':'Star','**':'Double Star','Cld':'Nebula',
}

def resolve(ident):
    url = "https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident=" + urllib.parse.quote(ident)
    try:
        txt = urllib.request.urlopen(url, timeout=12).read().decode('utf-8','replace')
    except Exception as e:
        return None
    m = re.search(r'^Object\s+(.+?)\s+---\s+(\S+)', txt, re.M)
    if not m:
        return None
    main_id, otype_code = m.group(1).strip(), m.group(2).strip()
    cm = re.search(r'Coordinates\(ICRS[^)]*\):\s*([\d.]+ [\d.]+ [\d.]+)\s+([+\-][\d.]+ [\d.]+ [\d.]+)', txt)
    if not cm:
        return None
    c = SkyCoord(cm.group(1) + ' ' + cm.group(2), unit=(u.hourangle, u.deg))
    out = {'desig': main_id, 'ra': round(c.ra.deg,5), 'dec': round(c.dec.deg,5),
           'otype': OTYPE.get(otype_code, otype_code)}
    sm = re.search(r'Angular size:\s*([\d.]+)', txt)
    if sm: out['size'] = round(float(sm.group(1)),1)
    for band in ('V','B'):
        fm = re.search(r'Flux '+band+r'\s*:\s*([-\d.]+)', txt)
        if fm: out['mag'] = fm.group(1); break
    nm = re.findall(r'NAME\s+([A-Za-z][\w\'\.\- ]*?)(?:\s{2,}|\n)', txt)
    nm = [n.strip() for n in nm if len(n.strip())>2]
    if nm: out['common'] = nm[0]
    return out

catalog = []
for n in range(1, 111):
    info = resolve(f'M {n}')
    if info:
        info['id'] = f'M{n}'; info['cat'] = 'M'; catalog.append(info)
        print(f"M{n}: {info['desig']} {info.get('common','')}")
    else:
        print(f"M{n}: FAILED")
    time.sleep(0.15)

for n, desig in CALDWELL.items():
    info = resolve(desig)
    if info:
        info['id'] = f'C{n}'; info['cat'] = 'C'; catalog.append(info)
        print(f"C{n} ({desig}): {info['desig']} {info.get('common','')}")
    else:
        print(f"C{n} ({desig}): FAILED")
    time.sleep(0.15)

# --- GI750 catalog: a curated list of 750 deep-sky targets, each scored 1–5 ---
# The list (id, name, score, RA, Dec, size, magnitude) ships as gi750.csv so this
# script needs no spreadsheet. Coordinates/size/magnitude/score come from the CSV;
# the object type and a common name are enriched from SIMBAD (best-effort).
GI_TYPE = {
    'G':'Galaxy','GiG':'Galaxy','GiC':'Galaxy','IG':'Galaxies','AGN':'Active Galaxy','Sy1':'Galaxy',
    'Sy2':'Galaxy','SyG':'Galaxy','rG':'Galaxy','LIN':'Galaxy','SBG':'Galaxy','EmG':'Galaxy','H2G':'Galaxy',
    'bCG':'Galaxy','GiP':'Galaxy Pair','PaG':'Galaxy Pair','ClG':'Galaxy Cluster','GrG':'Galaxy Group',
    'CGG':'Compact Galaxy Group',
    'GlC':'Globular Cluster','OpC':'Open Cluster','Cl*':'Star Cluster','As*':'Association',
    'MGr':'Moving Group','Cl?':'Star Cluster','C?*':'Star Cluster',
    'PN':'Planetary Nebula','PN?':'Planetary Nebula','SNR':'Supernova Remnant','SR?':'Supernova Remnant',
    'HII':'Emission Nebula','RNe':'Reflection Nebula','DNe':'Dark Nebula','GNe':'Nebula','Neb':'Nebula',
    'EmO':'Nebula','Cld':'Nebula','ISM':'Nebula','MoC':'Molecular Cloud','glb':'Dark Nebula',
    'SFR':'Star Forming Region','HH':'Herbig-Haro Object','cor':'Dark Nebula','bub':'Nebula',
    '*':'Star','**':'Double Star','V*':'Variable Star','WR*':'Wolf-Rayet Star',
}

def simbad_type(ident):
    """Best-effort (otype_label, common_name) from SIMBAD; (None, None) on failure."""
    url = "https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident=" + urllib.parse.quote(ident)
    try:
        txt = urllib.request.urlopen(url, timeout=12).read().decode('utf-8', 'replace')
    except Exception:
        return None, None
    m = re.search(r'^Object\s+(.+?)\s+---\s+(\S+)', txt, re.M)
    if not m:
        return None, None
    label = GI_TYPE.get(m.group(2).strip(), m.group(2).strip())
    nm = [n.strip() for n in re.findall(r'NAME\s+([A-Za-z][\w\'\.\- ]*?)(?:\s{2,}|\n)', txt) if len(n.strip()) > 2]
    return label, (nm[0] if nm else None)

def _sexa_to_deg(s, is_ra):
    sign = -1.0 if s.strip().startswith('-') else 1.0
    parts = re.findall(r'[\d.]+', s)
    a, m, sec = float(parts[0]), float(parts[1]), float(parts[2]) if len(parts) > 2 else 0.0
    deg = a + m / 60 + sec / 3600
    return round((deg * 15) if is_ra else (sign * deg), 5)

GI750_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gi750.csv')
with open(GI750_CSV, newline='') as f:
    for r in csv.DictReader(f):
        name = r['name'].strip()
        otype, common = simbad_type(name)
        entry = {'id': 'GI' + str(int(r['id'])), 'cat': 'G', 'desig': name,
                 'ra': _sexa_to_deg(r['ra'], True), 'dec': _sexa_to_deg(r['dec'], False),
                 'otype': otype or 'Deep Sky', 'score': int(r['score'])}
        try:
            if r.get('size'):
                entry['size'] = round(float(r['size']), 1)
        except ValueError:
            pass
        if r.get('mag'):
            entry['mag'] = r['mag']
        if common and common.lower() != name.lower():
            entry['common'] = common
        catalog.append(entry)
        print(f"{entry['id']}: {name} [{entry['otype']}] score {entry['score']}")
        time.sleep(0.1)

# Northern-hemisphere only — but GI750 is a hand-picked target list, so keep it whole.
catalog = [o for o in catalog if o['cat'] == 'G' or o['dec'] >= MIN_DEC]
with open(OUTPUT, 'w') as f:
    json.dump(catalog, f)
print(f"\nWrote {len(catalog)} objects "
      f"({sum(1 for o in catalog if o['cat']=='M')} Messier, "
      f"{sum(1 for o in catalog if o['cat']=='C')} Caldwell, "
      f"{sum(1 for o in catalog if o['cat']=='G')} GI750)")
