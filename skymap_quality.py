"""Depth, filter and image quality per target, measured from the pixels.

For every (rig, night, target) group of light frames a few sample frames are
measured (first / middle / last by filename):

* sky level per Bayer channel (bias-subtracted, per second, normalised to gain 100,
  an f/2 beam and a 3.76 um pixel, so different scopes compare directly)
  -> sky brightness and, for one-shot-colour cameras with no FILTER keyword, a
  filter class: a dual-band (Ha+OIII) filter makes the green sky ~6x darker, so
  frames whose green sky is below QUALITY_DUALBAND_SKY_G_MAX are classed dual-band;
* star FWHM (arcsec), elongation and count from a 2x2 superpixel luminance image.

Measurements are cached per file (path + mtime) in QUALITY_CACHE, so only new
frames are read on later runs. Nothing here writes to the image folders.

Per target the frames are then combined into an *effective* integration time
for each filter class: each frame counts t * (sky_ref / sky) * (fwhm_ref / fwhm)^2,
i.e. the time it would have needed under the rig's reference sky and seeing to
reach the same point-source signal-to-noise. Depth is reported relative to one
hour at the reference conditions: depth_mag = 1.25 * log10(effective_hours).
"""
import os, json, math, zlib, re
from collections import defaultdict
import numpy as np
from astropy.io import fits
from config import *
from skymap_io import read_xisf_header, _xisf_fitskw, should_exclude_target

QUALITY_SCHEMA = 3
QUALITY_CACHE = globals().get('QUALITY_CACHE', os.path.join(SCRIPT_DIR, 'quality_cache.json'))
SAMPLES_PER_GROUP = globals().get('QUALITY_SAMPLES_PER_GROUP', 3)
# Bias level in ADU per unit of the camera's OFFSET setting (ZWO ASI533: 2792 ADU at offset 70).
BIAS_PER_OFFSET = globals().get('QUALITY_BIAS_PER_OFFSET', {'ZWO ASI533MC Pro': 2792 / 70})
# Green-channel sky (ADU/s at gain 100) below which an OSC frame with no FILTER keyword
# is classed dual-band. Measured on the RASA 8 + ASI533 archive: dual-band nights
# 4-16, broadband nights 40-300.
DUALBAND_SKY_G_MAX = globals().get('QUALITY_DUALBAND_SKY_G_MAX', 30.0)
# Focal ratio for rigs whose headers lack FOCRATIO/APTDIA (ASIAIR writes neither).
DEFAULT_FRATIO = globals().get('QUALITY_DEFAULT_FRATIO', {'Rig A': 2.0})
# ZWO-style gain is in 0.1 dB: signal in ADU scales as 10**(gain/200).
GAIN_DB_TENTHS = globals().get('QUALITY_GAIN_DB_TENTHS', {'ZWO ASI533MC Pro'})

try:
    import lz4.block as _lz4
except ImportError:          # optional: only needed for lz4-compressed XISF (NINA's default)
    _lz4 = None


# ----------------------------------------------------------------------------- reading

def _read_xisf_pixels(path):
    hdr = read_xisf_header(path)
    m = re.search(r'<Image\b[^>]*>', hdr)
    if not m:
        return None, hdr
    tag = m.group(0)
    attr = dict(re.findall(r'(\w+)="([^"]*)"', tag))
    w, h, ch = (int(v) for v in attr['geometry'].split(':'))
    if ch != 1 or attr.get('sampleFormat') != 'UInt16':
        return None, hdr
    _, pos, size = attr['location'].split(':')
    with open(path, 'rb') as f:
        f.seek(int(pos)); raw = f.read(int(size))
    comp = attr.get('compression')
    if comp:
        codec, usize, *item = comp.split(':')
        usize = int(usize)
        if codec.startswith('zlib'):
            raw = zlib.decompress(raw)
        elif codec.startswith('lz4'):
            if _lz4 is None:
                return None, hdr
            raw = _lz4.decompress(raw, uncompressed_size=usize)
        else:
            return None, hdr
        if codec.endswith('+sh'):
            isz = int(item[0]) if item else 2
            raw = np.frombuffer(raw, np.uint8).reshape(isz, -1).T.tobytes()
    return np.frombuffer(raw, '<u2').reshape(h, w).astype(np.float32), hdr


def _load(path):
    """Return (pixels float32, keyword getter) or (None, None)."""
    if path.lower().endswith(('.fit', '.fits')):
        with fits.open(path, ignore_missing_end=True, memmap=False) as hd:
            if hd[0].data is None or hd[0].data.ndim != 2:
                return None, None
            h = hd[0].header
            return hd[0].data.astype(np.float32), (lambda k: h.get(k))
    data, hdr = _read_xisf_pixels(path)
    if data is None:
        return None, None
    def get(k):
        v = _xisf_fitskw(hdr, k)
        if v is None:
            return None
        v = v.strip().strip("'").strip()
        try:
            return float(v)
        except ValueError:
            return v
    return data, get


def _num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ----------------------------------------------------------------------------- measuring

def _block_bg(img, b=64):
    """Coarse background: median of b x b blocks, bilinearly upsampled."""
    H, W = img.shape
    nh, nw = H // b, W // b
    bm = np.median(img[:nh * b, :nw * b].reshape(nh, b, nw, b).transpose(0, 2, 1, 3).reshape(nh, nw, -1), axis=2)
    yi = np.clip((np.arange(H) + .5) / b - .5, 0, nh - 1.001); xi = np.clip((np.arange(W) + .5) / b - .5, 0, nw - 1.001)
    y0 = yi.astype(int); x0 = xi.astype(int); fy = (yi - y0)[:, None]; fx = (xi - x0)[None]
    return ((bm[y0][:, x0] * (1 - fx) + bm[y0][:, x0 + 1] * fx) * (1 - fy)
            + (bm[y0 + 1][:, x0] * (1 - fx) + bm[y0 + 1][:, x0 + 1] * fx) * fy), bm


def _stars(lum, bg, scale):
    """Star FWHM (arcsec), elongation and count from second moments of isolated peaks."""
    r = lum - bg
    sig = 1.4826 * np.median(np.abs(r[::4, ::4] - np.median(r[::4, ::4])))
    if not sig > 0:
        return None
    m = r > 10 * sig
    for dy in (-2, -1, 0, 1, 2):
        for dx in (-2, -1, 0, 1, 2):
            if dy or dx:
                m &= r >= np.roll(np.roll(r, dy, 0), dx, 1)
    m[:8] = m[-8:] = False; m[:, :8] = m[:, -8:] = False
    ys, xs = np.nonzero(m)
    n_det = len(ys)
    if n_det < 5:
        return {'nstars': n_det}
    sat = np.percentile(lum, 99.999) * 0.9
    peaks = r[ys, xs]
    order = np.argsort(-peaks)
    yy, xx = np.mgrid[-4:5, -4:5]
    fw, el = [], []
    for i in order:
        y, x = ys[i], xs[i]
        if lum[y, x] >= sat:
            continue
        c = r[y - 4:y + 5, x - 4:x + 5]
        if c.max() > c[4, 4] or (c[(np.abs(yy) == 4) | (np.abs(xx) == 4)] > 0.3 * c[4, 4]).any():
            continue                                   # neighbour in the box: not isolated
        c = np.clip(c, 0, None); s = c.sum()
        if s <= 0:
            continue
        my = (c * yy).sum() / s; mx = (c * xx).sum() / s
        vyy = (c * (yy - my) ** 2).sum() / s - 1 / 12; vxx = (c * (xx - mx) ** 2).sum() / s - 1 / 12
        vxy = (c * (yy - my) * (xx - mx)).sum() / s
        tr = vyy + vxx; det = math.sqrt(max(((vyy - vxx) / 2) ** 2 + vxy ** 2, 0))
        l1, l2 = tr / 2 + det, tr / 2 - det
        if l2 <= 0:
            continue
        fw.append(2.3548 * math.sqrt(tr / 2) * scale); el.append(1 - math.sqrt(l2 / l1))
        if len(fw) >= 150:
            break
    if len(fw) < 5:
        return {'nstars': n_det}
    return {'nstars': n_det, 'fwhm': round(float(np.median(fw)), 2), 'elong': round(float(np.median(el)), 3)}


def measure_frame(path, rig=''):
    data, get = _load(path)
    if data is None:
        return {'err': 'unreadable'}
    exptime = _num(get('EXPTIME')) or _num(get('EXPOSURE'))
    xpix = _num(get('XPIXSZ')); flen = _num(get('FOCALLEN')); xbin = _num(get('XBINNING'), 1.0) or 1.0
    if not exptime or not xpix or not flen:
        return {'err': 'missing exptime/pixel size/focal length'}
    scale = 206.265 * xpix * xbin / flen            # arcsec per pixel
    inst = str(get('INSTRUME') or '').strip()
    offset = _num(get('OFFSET'))
    bias = BIAS_PER_OFFSET.get(inst, None)
    bias = bias * offset if (bias is not None and offset is not None) else None
    bayer = str(get('BAYERPAT') or '').strip().upper()
    gain = _num(get('GAIN'))
    gnorm = 10 ** ((100 - gain) / 200) if (inst in GAIN_DB_TENTHS and gain is not None) else 1.0
    fratio = _num(get('FOCRATIO'))
    if not fratio and _num(get('APTDIA')):
        fratio = flen / _num(get('APTDIA'))
    if not fratio:
        fratio = DEFAULT_FRATIO.get(rig)
    if fratio:
        gnorm *= (fratio / 2.0) ** 2 * (3.76 / (xpix * xbin)) ** 2
    else:
        bias = None                      # unknown beam: don't report a sky level
    out = {'exptime': exptime, 'gain': gain, 'flen': flen, 'fratio': fratio, 'filter': (str(get('FILTER') or '').strip() or None)}
    if bayer in ('RGGB', 'BGGR', 'GRBG', 'GBRG'):
        H, W = data.shape; data = data[:H // 2 * 2, :W // 2 * 2]
        ch = {(0, 0): data[0::2, 0::2], (0, 1): data[0::2, 1::2], (1, 0): data[1::2, 0::2], (1, 1): data[1::2, 1::2]}
        pos = {bayer[0]: (0, 0), bayer[1]: (0, 1), bayer[3]: (1, 1)}
        R, B = ch[pos['R']], ch[pos['B']]
        greens = [ch[k] for k in ch if k not in (pos['R'], pos['B'])]
        G = (greens[0] + greens[1]) / 2
        lum = R + greens[0] + greens[1] + B
        scale *= 2
        if bias is not None:
            sky = {k: float(np.percentile(_block_bg(v - bias, 32)[1], 20)) / exptime * gnorm
                   for k, v in (('r', R), ('g', G), ('b', B))}
            out.update(sky_r=round(sky['r'], 4), sky_g=round(sky['g'], 4), sky_b=round(sky['b'], 4))
            if sky['g'] > 0:
                out['rg'] = round(sky['r'] / sky['g'], 3)
    else:
        lum = data
        if bias is not None:
            out['sky'] = round(float(np.percentile(_block_bg(data - bias, 64)[1], 20)) / exptime * gnorm, 4)
    bg, _ = _block_bg(lum, 32)
    st = _stars(lum, bg, scale)
    if st:
        out.update(st)
    return out


# ----------------------------------------------------------------------------- grouping

def _group_key(path, entry):
    return (entry.get('rig') or '', os.path.dirname(path), entry.get('date_obs') or '', entry.get('target') or '')


def update_quality_cache(cache, log=print):
    """Measure sample frames of every (rig, folder, night, target) group. Returns the quality cache."""
    try:
        with open(QUALITY_CACHE) as f:
            qc = json.load(f)
    except (OSError, ValueError):
        qc = {}
    groups = defaultdict(list)
    for p, e in cache.items():
        if e.get('skip') or 'ra' not in e or not os.path.exists(p) or should_exclude_target(e.get('target', '')):
            continue
        groups[_group_key(p, e)].append(p)
    todo = []
    for paths in groups.values():
        paths.sort()
        n = len(paths)
        idx = sorted({0, n // 2, n - 1})[:SAMPLES_PER_GROUP] if n > 1 else [0]
        for i in idx:
            p = paths[i]
            mt = cache[p].get('mtime')
            q = qc.get(p)
            if not (q and q.get('ver') == QUALITY_SCHEMA and q.get('mtime') == mt):
                todo.append(p)
    if todo:
        log(f"\n📏 Measuring depth/quality on {len(todo)} sample frame(s) …")
    for k, p in enumerate(todo, 1):
        try:
            r = measure_frame(p, cache[p].get('rig') or '')
        except Exception as ex:        # never let one bad file stop the map
            r = {'err': f'{type(ex).__name__}: {ex}'}
        r.update(ver=QUALITY_SCHEMA, mtime=cache[p].get('mtime'))
        qc[p] = r
        if k % 25 == 0 or k == len(todo):
            log(f"   {k}/{len(todo)}")
            with open(QUALITY_CACHE, 'w') as f:
                json.dump(qc, f)
    # forget measurements of files that left the main cache
    for p in [p for p in qc if p not in cache]:
        del qc[p]
    with open(QUALITY_CACHE, 'w') as f:
        json.dump(qc, f)
    return qc


def _filter_class(sky_g, header_filter):
    if header_filter:
        return header_filter
    if sky_g is None:
        return 'unknown'
    return 'Dual-band (est.)' if sky_g < DUALBAND_SKY_G_MAX else 'Broadband (est.)'


def group_metrics(cache, qc):
    """Per-group medians of the sampled frames: {group_key: {...}} plus per-file group lookup."""
    groups = defaultdict(list)
    for p, e in cache.items():
        if e.get('skip') or 'ra' not in e or should_exclude_target(e.get('target', '')):
            continue
        groups[_group_key(p, e)].append(p)
    gm, file_group = {}, {}
    for key, paths in groups.items():
        qs = [qc[p] for p in paths if p in qc and 'err' not in qc[p]]
        for p in paths:
            file_group[p] = key
        if not qs:
            continue
        def med(k):
            v = [q[k] for q in qs if q.get(k) is not None]
            return float(np.median(v)) if v else None
        hdr_f = cache[paths[0]].get('filter')
        cls = _filter_class(med('sky_g'), hdr_f)
        flen = med('flen')
        gm[key] = {'optics': (key[0], int(round(flen, -1))) if flen else (key[0], 0), 'cls': cls, 'fwhm': med('fwhm'), 'elong': med('elong'), 'nstars': med('nstars'),
                   'sky': med('sky_g') if med('sky_g') is not None else med('sky'), 'rg': med('rg')}
    return gm, file_group


def references(gm):
    """Reference sky (10th percentile) per (rig, filter class) and seeing (median FWHM)
    per optics (rig, focal length): a 3910 mm f/11 frame is never graded against f/2."""
    sky, fw = defaultdict(list), defaultdict(list)
    for (rig, *_), g in gm.items():
        if g['sky'] and g['sky'] > 0:
            sky[(rig, g['cls'])].append(g['sky'])
        if g['fwhm']:
            fw[g['optics']].append(g['fwhm'])
    return ({k: float(np.percentile(v, 10)) for k, v in sky.items()},
            {k: float(np.median(v)) for k, v in fw.items()})


def cluster_quality(members, gm, file_group, refs):
    """Depth / filter / quality summary for one map target (cluster)."""
    sky_ref, fw_ref = refs
    classes = {}
    fws, rels, els, skyrel, nights, measured = [], [], [], [], set(), 0
    for m in members:
        key = file_group.get(m.get('path'))
        g = gm.get(key) if key else None
        t = _num(m.get('exptime'), 0.0)
        rig = m.get('rig') or ''
        # Frames that couldn't be measured (mostly processed stacks / masters, whose
        # exposure time already repeats their subs) are listed but never add depth.
        cls = g['cls'] if g else 'not measured'
        slot = classes.setdefault(cls, {'frames': 0, 'seconds': 0.0, 'eff_seconds': 0.0})
        slot['frames'] += 1; slot['seconds'] += t
        if g:
            measured += 1
            w = 1.0
            sr = sky_ref.get((rig, cls))
            if sr and g['sky'] and g['sky'] > 0:
                w *= sr / g['sky']; skyrel.append(g['sky'] / sr)
            fr = fw_ref.get(g['optics'])
            if fr and g['fwhm']:
                rels.append(g['fwhm'] / fr)
                w *= (fr / g['fwhm']) ** 2; fws.append(g['fwhm'])
            if g['elong'] is not None:
                els.append(g['elong'])
            slot['eff_seconds'] += t * min(w, 4.0)
            nights.add(key[2])
    for k, v in classes.items():
        v['seconds'] = round(v['seconds']); v['eff_seconds'] = round(v['eff_seconds'])
        eh = v['eff_seconds'] / 3600
        v['depth_mag'] = round(1.25 * math.log10(eh), 2) if (eh > 0 and k != 'not measured') else None
    if not measured:
        return {'classes': classes, 'measured': 0}
    fwhm = float(np.median(fws)) if fws else None
    elong = float(np.median(els)) if els else None
    rel = float(np.median(rels)) if rels else None
    grade = None
    if rel:
        if rel <= 0.92 and (elong is None or elong < 0.2):
            grade = 'A'
        elif rel > 1.2 or (elong is not None and elong > 0.3):
            grade = 'C'
        else:
            grade = 'B'
    scored = [kv for kv in classes.items() if kv[1]['depth_mag'] is not None]
    best = max(scored, key=lambda kv: kv[1]['eff_seconds']) if scored else (None, {'depth_mag': None})
    return {'classes': classes, 'measured': measured,
            'fwhm': round(fwhm, 2) if fwhm else None,
            'fwhm_rel': round(rel, 2) if rel else None,
            'elong': round(elong, 3) if elong is not None else None,
            'sky_rel': round(float(np.median(skyrel)), 2) if skyrel else None,
            'grade': grade, 'depth_mag': best[1]['depth_mag'], 'depth_class': best[0]}
