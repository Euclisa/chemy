"""Utilities for synonym filtering audits and regeneration.

Usage examples:
    python -m scripts.run.synonyms --audit
    python -m scripts.run.synonyms --regenerate
"""
import json
import sys

from ._base import base_parser, make_chemy


def main():
    parser = base_parser(description="Audit or regenerate compound synonym filtering.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--audit", action="store_true", help="Print the current synonym audit report")
    mode.add_argument(
        "--regenerate",
        action="store_true",
        help="Migrate legacy decisions and rewrite mapped compounds through the synonym filter",
    )
    args = parser.parse_args()

    chemy = make_chemy(args)
    if args.regenerate:
        audit = chemy.compounds.regenerate_mapped_synonyms()
    else:
        audit = chemy.compounds.audit_synonym_filtering()

    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    sys.exit(main())
