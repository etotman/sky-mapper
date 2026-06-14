"""Clustering, footprints, folder tree, and archive statistics."""
import os, math
from collections import Counter
from config import *
from skymap_io import should_exclude_target

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
            # Forward slashes so the browser can join "adir/file" consistently on
            # Windows too; the header endpoint normalizes via os.path.realpath.
            node['adir'] = os.path.dirname(filepath).replace(os.sep, '/')

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


