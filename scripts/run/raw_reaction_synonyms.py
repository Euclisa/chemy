"""Mine and apply synonym candidates from unmapped raw reaction names.

Usage examples:
    python -m scripts.run.raw_reaction_synonyms mine --preset wiki_crc_rp
    python -m scripts.run.raw_reaction_synonyms mine --limit 100 --workers 4
    python -m scripts.run.raw_reaction_synonyms apply
    python -m scripts.run.raw_reaction_synonyms apply --audit-jsonl data/raw_reactions/synonym_enrichment/candidates.jsonl
"""
import json
import sys

from ._base import base_parser, make_chemy
from scripts.ops.raw_reactions.synonym_enrichment import (
    DEFAULT_BATCH_SIZE,
    DEFAULT_MIN_VOTES,
    DEFAULT_ROUNDS,
)


def main():
    parser = base_parser(
        description="Mine raw reaction unmapped names and enrich compound synonyms."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    mine = subparsers.add_parser(
        "mine",
        help="Generate synonym enrichment audit files without mutating compounds",
    )
    mine.add_argument("--preset", help="Raw reaction preset to scan; omit to scan all presets")
    mine.add_argument("--model", help="LLM model to use; defaults to the raw reaction default")
    mine.add_argument("--rounds", type=int, default=DEFAULT_ROUNDS)
    mine.add_argument("--min-votes", type=int, default=DEFAULT_MIN_VOTES)
    mine.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    mine.add_argument("--limit", type=int, help="Only mine the top N unmapped names")
    mine.add_argument("--audit-jsonl", help="Output JSONL audit path")
    mine.add_argument("--audit-md", help="Output Markdown audit path")

    apply = subparsers.add_parser(
        "apply",
        help="Apply approved synonym additions from an audit JSONL file",
    )
    apply.add_argument("--audit-jsonl", help="Input JSONL audit path")

    args = parser.parse_args()
    chemy = make_chemy(args)
    enricher = chemy.raw_reaction_synonyms

    if args.command == "mine":
        summary = enricher.mine(
            args.preset,
            model=args.model,
            rounds=args.rounds,
            min_votes=args.min_votes,
            batch_size=args.batch_size,
            max_workers=args.workers,
            limit=args.limit,
            audit_jsonl_fn=args.audit_jsonl,
            audit_md_fn=args.audit_md,
        )
    else:
        summary = enricher.apply(args.audit_jsonl)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    sys.exit(main())
