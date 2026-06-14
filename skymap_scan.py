"""Scan rig folders, prune deleted files, recluster, and rebuild the map."""
import os, glob, threading
from collections import Counter
from config import *
from skymap_io import (load_cache, save_cache, get_fits_data, get_xisf_data,
                       is_mosaic, should_exclude_target)
from skymap_catalog import (cluster_by_distance, build_cluster_payload,
                            build_folder_data, build_stats)
from skymap_web import (resolve_observer_site, ensure_constellation_data,
                        generate_aladin_html)

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
        excluded_lower = {f.lower() for f in EXCLUDE_FOLDERS}
        for ext in ('*.fits', '*.fit', '*.xisf'):
            for path in glob.glob(os.path.join(search_dir, '**', ext), recursive=True):
                parts = os.path.normpath(path).split(os.sep)
                if not any(p.lower() in excluded_lower for p in parts):
                    all_files.append(path)
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
