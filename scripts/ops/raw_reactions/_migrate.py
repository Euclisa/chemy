"""
Migration to the run-based / global-verdict layout.

Previous layout:
    data/raw_reactions/<preset>/raw.jsonl
    data/raw_reactions/<preset>/verdict.jsonl   (one per preset)
    data/raw_reactions/repairs.jsonl

New layout:
    data/raw_reactions/<preset>/raw_1.jsonl     (run-based; old raw.jsonl becomes run 1)
    data/raw_reactions/verdict/verdict_1.jsonl  (global vault; all preset verdicts merged)
    data/raw_reactions/repairs.jsonl            (unchanged)

Run once:
    python -m scripts.ops.raw_reactions._migrate --data-dir data
"""

import json
import os
import shutil


def _iter_preset_dirs(raw_dir):
    """Yield names of preset subdirectories (those that are plain directories)."""
    for name in sorted(os.listdir(raw_dir)):
        path = os.path.join(raw_dir, name)
        if os.path.isdir(path) and name != 'verdict':
            yield name, path


def migrate(data_dir: str, dry_run: bool = False):
    raw_dir = os.path.join(data_dir, "raw_reactions")
    if not os.path.isdir(raw_dir):
        print(f"Nothing to migrate: {raw_dir} does not exist.")
        return

    # --- migrate raw files: raw.jsonl -> raw_1.jsonl ---
    for preset_name, preset_dir in _iter_preset_dirs(raw_dir):
        src = os.path.join(preset_dir, "raw.jsonl")
        dst = os.path.join(preset_dir, "raw_1.jsonl")
        if not os.path.exists(src):
            continue
        if os.path.exists(dst):
            print(f"  skip (raw_1.jsonl already exists): {preset_name}/raw_1.jsonl")
            continue
        print(f"  {preset_name}/raw.jsonl -> {preset_name}/raw_1.jsonl")
        if not dry_run:
            shutil.move(src, dst)

    # --- merge all preset-local verdict.jsonl into global vault ---
    verdict_dir = os.path.join(raw_dir, "verdict")
    verdict_shard = os.path.join(verdict_dir, "verdict_1.jsonl")

    all_entries = []
    migrated_presets = []

    for preset_name, preset_dir in _iter_preset_dirs(raw_dir):
        src = os.path.join(preset_dir, "verdict.jsonl")
        if not os.path.exists(src):
            continue
        with open(src) as f:
            entries = [json.loads(line) for line in f if line.strip()]
        print(f"  collecting {len(entries)} entries from {preset_name}/verdict.jsonl")
        all_entries.extend(entries)
        migrated_presets.append((preset_name, src))

    if not all_entries:
        print("  no preset-local verdict files found; nothing to merge.")
        return

    if os.path.exists(verdict_shard):
        print(f"  skip (global verdict shard already exists): verdict/verdict_1.jsonl")
    else:
        print(f"  writing {len(all_entries)} entries -> verdict/verdict_1.jsonl")
        if not dry_run:
            os.makedirs(verdict_dir, exist_ok=True)
            with open(verdict_shard, 'w') as f:
                for entry in all_entries:
                    f.write(json.dumps(entry) + '\n')

    # --- remove preset-local verdict files after successful migration ---
    for preset_name, src in migrated_presets:
        backup = src + ".migrated_backup"
        print(f"  archiving {preset_name}/verdict.jsonl -> {preset_name}/verdict.jsonl.migrated_backup")
        if not dry_run:
            shutil.move(src, backup)

    print(f"Migration complete. {len(all_entries)} verdict entries in global vault.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data", help="Path to the data directory")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would happen without moving files")
    args = parser.parse_args()

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Migrating {args.data_dir}/raw_reactions/ ...")
    migrate(args.data_dir, dry_run=args.dry_run)
    print("Done.")
