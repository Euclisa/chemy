"""
Entry point for raw reactions extraction and/or validation.

Usage examples:
    # fetch + validate in one step (default, continues latest run)
    python -m scripts.run.fetch_reactions wiki_uncommon_rp

    # start a fresh run from scratch
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --run new

    # use a specific existing run
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --run 2

    # fetch only
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --fetch-only

    # validate only (aggregates all runs for the preset)
    python -m scripts.run.fetch_reactions wiki_uncommon_rp --validate-only

    # repair borderline rejected verdicts (global, no preset needed)
    python -m scripts.run.fetch_reactions --repair-only

    # fetch + validate for every preset
    python -m scripts.run.fetch_reactions --all

    # parallel workers
    python -m scripts.run.fetch_reactions top_rare_rp --workers 4
"""
import sys

from ._base import base_parser, make_chemy
from scripts.ops.raw_reactions.validator import ACCEPT_THR
from scripts.ops.raw_reactions.presets import DEFAULT_PRESET_NAMES, get_preset


def main():
    parser = base_parser(description="Fetch and/or validate raw LLM reactions.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("preset", nargs="?",
                       help="Preset to run, in the form <base>_<rp|r|p>")
    group.add_argument("--all", action="store_true",
                       help="Run the default preset suite in sequence")
    group.add_argument("--repair-only", action="store_true",
                       help="Repair borderline rejected verdicts (global, no preset needed)")

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--fetch-only", action="store_true",
                      help="Only fetch raw reactions, skip validation")
    mode.add_argument("--validate-only", action="store_true",
                      help="Only validate all raw runs for the preset")

    parser.add_argument("--run", metavar="{new|INT}",
                        help="Run to use: 'new' creates a fresh run, INT resumes a specific run, "
                             "omit to continue the latest run (only valid when fetching)")
    parser.add_argument("--with-repair", action="store_true",
                        help="Run global repair after fetch/validation")
    parser.add_argument("--repair-lower-bound", type=float, default=0.1,
                        help="Minimum confidence for repair candidates (default: 0.1)")
    parser.add_argument("--acceptance-threshold", type=float, default=ACCEPT_THR,
                        help=f"Confidence threshold for acceptance (default: {ACCEPT_THR})")

    args = parser.parse_args()

    if args.repair_only:
        if args.run:
            parser.error("--run is not meaningful with --repair-only")
        if args.fetch_only or args.validate_only:
            parser.error("--fetch-only/--validate-only are not applicable with --repair-only")

    if args.run and args.validate_only:
        parser.error("--run is not meaningful with --validate-only")

    if args.with_repair and args.fetch_only:
        parser.error("--with-repair requires validation output; use it without --fetch-only")

    chemy = make_chemy(args)
    pipeline = chemy.raw_reactions

    if args.repair_only:
        pipeline.repair(
            args.workers,
            lower_bound=args.repair_lower_bound,
            acceptance_threshold=args.acceptance_threshold,
        )
        return

    preset_names = DEFAULT_PRESET_NAMES if args.all else [args.preset]

    for name in preset_names:
        try:
            get_preset(name)
        except ValueError as e:
            parser.error(str(e))

    for name in preset_names:
        if args.fetch_only:
            pipeline.fetch(name, args.workers, run=args.run)
        elif args.validate_only:
            pipeline.validate(
                name,
                args.workers,
                acceptance_threshold=args.acceptance_threshold,
            )
        else:
            pipeline.run(
                name,
                args.workers,
                acceptance_threshold=args.acceptance_threshold,
                run=args.run,
            )

    if args.with_repair:
        pipeline.repair(
            args.workers,
            lower_bound=args.repair_lower_bound,
            acceptance_threshold=args.acceptance_threshold,
        )


if __name__ == "__main__":
    sys.exit(main())
