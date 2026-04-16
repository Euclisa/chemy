import json
from collections import defaultdict
from pathlib import Path

from .helpers import build_arg_parser, load_config, read_jsonl


def _mean(values):
    return sum(values) / len(values) if values else 0.0


def _rate(numer, denom):
    return numer / denom if denom else 0.0


def _init_model_summary():
    return {
        "generation_runs": 0,
        "candidate_n": 0,
        "schema_success_rate": 0.0,
        "parsed_rate": 0.0,
        "target_compound_rate": 0.0,
        "target_position_rate": 0.0,
        "balance_success_rate": 0.0,
        "unique_rid_n": 0,
        "unique_output_id_n": 0,
        "duplicate_rid_rate": 0.0,
        "gold_audited_n": 0,
        "gold_precision": 0.0,
        "gold_accepted_n": 0,
        "gold_accepted_unique_rid_n": 0,
        "input_tokens": 0,
        "completion_tokens": 0,
        "accepted_unique_rids_per_1k_tokens": 0.0,
    }


def build_summary(config, *, data_dir):
    data_dir = Path(data_dir)
    generations = read_jsonl(data_dir / "generations.jsonl")
    candidates = read_jsonl(data_dir / "candidates.jsonl")
    gold = read_jsonl(data_dir / "gold_audits.jsonl")
    gold_by_candidate = {entry["candidate_id"]: entry for entry in gold}

    by_model = defaultdict(_init_model_summary)
    parse_denominators = defaultdict(int)
    parse_accepts = defaultdict(int)
    candidates_by_model = defaultdict(list)

    for generation in generations:
        model = generation["model"]
        by_model[model]["generation_runs"] += 1
        usage = generation.get("usage") or {}
        by_model[model]["input_tokens"] += usage.get("input_tokens", 0)
        by_model[model]["completion_tokens"] += usage.get("completion_tokens", 0)
        stats = generation.get("parse_stats") or {}
        reviewed_stats = generation.get("reviewed_parse_stats") or {}
        source_stats = reviewed_stats if reviewed_stats.get("accepted") else stats
        accepted = source_stats.get("accepted", 0)
        parse_accepts[model] += accepted
        parse_denominators[model] += accepted
        parse_denominators[model] += source_stats.get("invalid_schema", 0)
        parse_denominators[model] += source_stats.get("malformed_json", 0)

    for candidate in candidates:
        model = candidate["model"]
        by_model[model]["candidate_n"] += 1
        candidates_by_model[model].append(candidate)

    for model, rows in candidates_by_model.items():
        parsed = [row for row in rows if row.get("parsed")]
        target_present = [row for row in rows if row.get("target_compound_present")]
        target_correct = [row for row in rows if row.get("target_position_correct")]
        balanced = [row for row in rows if row.get("balance_success")]
        rids = [row.get("rid") for row in rows if row.get("rid")]
        unique_rids = set(rids)
        output_ids = {row.get("output_id") for row in rows if row.get("output_id")}

        by_model[model]["schema_success_rate"] = _rate(
            parse_accepts[model],
            parse_denominators[model],
        )
        by_model[model]["parsed_rate"] = _rate(len(parsed), len(rows))
        by_model[model]["target_compound_rate"] = _rate(len(target_present), len(rows))
        by_model[model]["target_position_rate"] = _rate(len(target_correct), len(rows))
        by_model[model]["balance_success_rate"] = _rate(len(balanced), len(rows))
        by_model[model]["unique_rid_n"] = len(unique_rids)
        by_model[model]["unique_output_id_n"] = len(output_ids)
        by_model[model]["duplicate_rid_rate"] = 1 - _rate(len(unique_rids), len(rids))

        audited = []
        accepted = []
        for row in rows:
            audit = gold_by_candidate.get(row["candidate_id"])
            if audit is None:
                continue
            audited.append(audit)
            if audit.get("gold_confidence", 0.0) >= config["acceptance_threshold"]:
                accepted.append((row, audit))
        by_model[model]["gold_audited_n"] = len(audited)
        by_model[model]["gold_accepted_n"] = len(accepted)
        by_model[model]["gold_precision"] = _rate(len(accepted), len(audited))
        by_model[model]["gold_accepted_unique_rid_n"] = len(
            {row.get("rid") for row, _ in accepted if row.get("rid")}
        )
        total_tokens = (
            by_model[model]["input_tokens"] + by_model[model]["completion_tokens"]
        )
        by_model[model]["accepted_unique_rids_per_1k_tokens"] = (
            by_model[model]["gold_accepted_unique_rid_n"] / total_tokens * 1000
            if total_tokens else 0.0
        )

    by_band = defaultdict(lambda: defaultdict(list))
    for candidate in candidates:
        key = f"{candidate.get('preset')}:{candidate.get('complexity_band')}"
        by_band[key][candidate["model"]].append(candidate)

    bands = {}
    for band_key, model_rows in by_band.items():
        bands[band_key] = {}
        for model, rows in model_rows.items():
            bands[band_key][model] = {
                "candidate_n": len(rows),
                "parsed_rate": _mean([1.0 if row.get("parsed") else 0.0 for row in rows]),
                "target_position_rate": _mean(
                    [1.0 if row.get("target_position_correct") else 0.0 for row in rows]
                ),
                "balance_success_rate": _mean(
                    [1.0 if row.get("balance_success") else 0.0 for row in rows]
                ),
            }

    return {
        "models": dict(by_model),
        "bands": bands,
        "inputs": {
            "generations": len(generations),
            "candidates": len(candidates),
            "gold_audits": len(gold),
        },
    }


def main(argv=None):
    parser = build_arg_parser("Score fresh-generation model bakeoff outputs.")
    parser.add_argument(
        "--out",
        default=None,
        help="Optional JSON path for summary output.",
    )
    args = parser.parse_args(argv)
    config = load_config(args.config)
    summary = build_summary(config, data_dir=args.data_dir)
    payload = json.dumps(summary, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            f.write(payload)
            f.write("\n")
    print(payload)


if __name__ == "__main__":
    main()
