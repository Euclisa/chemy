"""
Entry point for raw reactions extraction and/or validation.

Usage examples:
    # fetch + validate in one step (default)
    python -m scripts.run.fetch_reactions wiki_uncommon_rp

    # fetch only
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --fetch-only

    # validate only (when raw file already exists)
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --validate-only

    # fetch + validate for every preset
    python -m scripts.run.fetch_reactions --all

    # parallel workers
    python -m scripts.run.fetch_reactions top_rare_rp --workers 4
"""
import sys

from ._base import base_parser, make_chemy
from scripts.ops.raw_reactions.presets import PRESETS


def main():
    parser = base_parser(description="Fetch and/or validate raw LLM reactions.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("preset", nargs="?", choices=list(PRESETS),
                       help="Preset to run")
    group.add_argument("--all", action="store_true",
                       help="Run all presets in sequence")

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--fetch-only", action="store_true",
                      help="Only fetch raw reactions, skip validation")
    mode.add_argument("--validate-only", action="store_true",
                      help="Only validate an existing raw file")

    args = parser.parse_args()

    chemy = make_chemy(args)
    pipeline = chemy.raw_reactions

    preset_names = list(PRESETS) if args.all else [args.preset]

    for name in preset_names:
        if args.fetch_only:
            pipeline.fetch(name, args.workers)
        elif args.validate_only:
            pipeline.validate(name, args.workers)
        else:
            pipeline.run(name, args.workers)


if __name__ == "__main__":
    sys.exit(main())
