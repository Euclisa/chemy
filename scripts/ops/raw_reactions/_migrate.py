"""
One-time migration from the legacy flat raw_reactions/ layout to per-preset subdirs.

Legacy layout:
    data/raw_reactions/raw_reactions.jsonl              -> default/raw.jsonl
    data/raw_reactions/wiki_raw_reactions.jsonl         -> wiki_uncommon_rp/raw.jsonl
    data/raw_reactions/top_rare_raw_reactions.jsonl     -> top_rare_rp/raw.jsonl
    data/raw_reactions/wiki_products_raw_reactions.jsonl    -> wiki_uncommon_p/raw.jsonl
    data/raw_reactions/annotated_products_raw_reactions.jsonl -> annotated_uncommon_p/raw.jsonl
    data/raw_reactions/raw_reactions_verdict.jsonl      -> partitioned per preset

Run once:
    python -m scripts.ops.raw_reactions._migrate --data-dir data
"""

import json
import os
import shutil


_LEGACY_RAW_FILES = {
    "raw_reactions.jsonl": "default",
    "wiki_raw_reactions.jsonl": "wiki_uncommon_rp",
    "top_rare_raw_reactions.jsonl": "top_rare_rp",
    "wiki_products_raw_reactions.jsonl": "wiki_uncommon_p",
    "annotated_products_raw_reactions.jsonl": "annotated_uncommon_p",
}

_LEGACY_VERDICT_FILE = "raw_reactions_verdict.jsonl"


def migrate(data_dir: str, dry_run: bool = False):
    raw_dir = os.path.join(data_dir, "raw_reactions")
    if not os.path.isdir(raw_dir):
        print(f"Nothing to migrate: {raw_dir} does not exist.")
        return

    # --- migrate raw files ---
    for legacy_name, preset_name in _LEGACY_RAW_FILES.items():
        src = os.path.join(raw_dir, legacy_name)
        if not os.path.exists(src):
            print(f"  skip (not found): {legacy_name}")
            continue
        dest_dir = os.path.join(raw_dir, preset_name)
        dest = os.path.join(dest_dir, "raw.jsonl")
        print(f"  {legacy_name} -> {preset_name}/raw.jsonl")
        if not dry_run:
            os.makedirs(dest_dir, exist_ok=True)
            shutil.move(src, dest)

    # --- partition verdict file ---
    legacy_verdict = os.path.join(raw_dir, _LEGACY_VERDICT_FILE)
    if not os.path.exists(legacy_verdict):
        print(f"  skip (not found): {_LEGACY_VERDICT_FILE}")
        return

    # Build a map: cid -> which preset's raw file contains that cid.
    # Read from legacy locations (before files are moved).
    cid_to_presets = {}
    for legacy_name, preset_name in _LEGACY_RAW_FILES.items():
        raw_fn = os.path.join(raw_dir, legacy_name)
        if not os.path.exists(raw_fn):
            # Already moved (real run, second pass) — fall back to new location.
            raw_fn = os.path.join(raw_dir, preset_name, "raw.jsonl")
        if not os.path.exists(raw_fn):
            continue
        with open(raw_fn) as f:
            for line in f:
                entry = json.loads(line)
                cid = entry.get('cid')
                if cid is not None:
                    cid_to_presets.setdefault(cid, set()).add(preset_name)

    verdict_entries = {}  # preset_name -> list of entries
    unrouted = []

    with open(legacy_verdict) as f:
        for line in f:
            entry = json.loads(line.strip())
            cid = entry.get('cid')
            destinations = cid_to_presets.get(cid, set())
            if not destinations:
                unrouted.append(entry)
            for preset_name in destinations:
                verdict_entries.setdefault(preset_name, []).append(entry)

    for preset_name, entries in verdict_entries.items():
        dest_dir = os.path.join(raw_dir, preset_name)
        dest = os.path.join(dest_dir, "verdict.jsonl")
        print(f"  verdict -> {preset_name}/verdict.jsonl ({len(entries)} entries)")
        if not dry_run:
            os.makedirs(dest_dir, exist_ok=True)
            with open(dest, 'w') as f:
                for entry in entries:
                    f.write(json.dumps(entry) + '\n')

    if unrouted:
        print(f"  WARNING: {len(unrouted)} verdict entries could not be routed to any preset (no matching raw cid).")

    if not dry_run:
        backup = legacy_verdict + ".migrated_backup"
        shutil.move(legacy_verdict, backup)
        print(f"  Legacy verdict file backed up to: {backup}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data", help="Path to the data directory")
    parser.add_argument("--dry-run", action="store_true", help="Print what would happen without moving files")
    args = parser.parse_args()

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Migrating {args.data_dir}/raw_reactions/ ...")
    migrate(args.data_dir, dry_run=args.dry_run)
    print("Done.")
