"""
Thermochemistry pipeline entry point.

Usage examples:
    python -m scripts.run.thermo collect-llm --workers 4
    python -m scripts.run.thermo extract-burcat --burcat-file BURCAT.THR.txt
    python -m scripts.run.thermo solve --mode all
    python -m scripts.run.thermo validate --mode hybrid
    python -m scripts.run.thermo experiment --policy core
"""
import sys

from ._base import base_parser, make_chemy


def main():
    parser = base_parser(description="Run thermochemistry extraction and solves.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser(
        "collect-llm",
        help="Collect LLM reaction thermochemistry estimates.",
    )

    extract = subparsers.add_parser(
        "extract-burcat",
        help="Extract BURCAT formation references into thermo evidence.",
    )
    extract.add_argument(
        "--burcat-file",
        default="BURCAT.THR.txt",
        metavar="PATH",
        help="BURCAT thermo file path (default: BURCAT.THR.txt)",
    )

    solve = subparsers.add_parser(
        "solve",
        help="Solve direct formation thermochemistry from evidence.",
    )
    solve.add_argument(
        "--mode",
        choices=["llm", "hybrid", "all"],
        default="hybrid",
        help="Solve mode (default: hybrid)",
    )

    validate = subparsers.add_parser(
        "validate",
        help="Compare solved formation values with BURCAT references.",
    )
    validate.add_argument(
        "--mode",
        choices=["llm", "hybrid"],
        default="hybrid",
        help="Solve output to validate (default: hybrid)",
    )

    experiment = subparsers.add_parser(
        "experiment",
        help="Run thermo reference-policy experiments.",
    )
    experiment.add_argument(
        "--policy",
        choices=["elements-only", "burcat-all", "burcat-split", "core"],
        default="core",
        help="Reference policy to run (default: core)",
    )
    experiment.add_argument(
        "--train-fraction",
        type=float,
        default=0.8,
        help="Train fraction for burcat-split (default: 0.8)",
    )
    experiment.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Deterministic split seed for burcat-split (default: 42)",
    )
    experiment.add_argument(
        "--name",
        default=None,
        help="Optional experiment directory name for non-core policies.",
    )

    args = parser.parse_args()
    chemy = make_chemy(args)

    if args.command == "collect-llm":
        chemy.thermo_llm.get_reactions_thermo(args.workers)
    elif args.command == "extract-burcat":
        chemy.thermo_burcat.extract_formation_references(args.burcat_file)
    elif args.command == "solve":
        modes = ["llm", "hybrid"] if args.mode == "all" else [args.mode]
        reaction_estimates = chemy.thermo_llm.normalize_reaction_estimates()
        for mode in modes:
            chemy.thermo_experiments.run_legacy_solve(mode, reaction_estimates)
    elif args.command == "validate":
        chemy.thermo_experiments.validate_legacy_solve(mode=args.mode)
    elif args.command == "experiment":
        reaction_estimates = chemy.thermo_llm.normalize_reaction_estimates()
        chemy.thermo_experiments.run_policy(
            args.policy,
            train_fraction=args.train_fraction,
            seed=args.seed,
            name=args.name,
            reaction_estimates=reaction_estimates,
        )


if __name__ == "__main__":
    sys.exit(main())
