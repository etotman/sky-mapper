"""File readers (FITS/XISF), header parsing, target/calibration helpers, cache."""
import os, re, math, json, warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from config import *

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
    target_lower = target_name.strip().lower()
    if target_lower in EXCLUDE_TARGETS:            # exact match (solar-system bodies)
        return True
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

            # Prefer the catalog target coordinates (OBJCTRA/OBJCTDEC) over the
            # mount-reported RA/DEC: the latter carries the mount's pointing error
            # (seen up to ~0.5° on un-plate-solved frames). A solved WCS, when
            # present, still wins over both (handled below).
            ra_val  = _first(header, 'OBJCTRA', 'RA')
            dec_val = _first(header, 'OBJCTDEC', 'DEC')
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

