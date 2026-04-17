import json
from collections import defaultdict
from pathlib import Path

from .helpers import build_arg_parser, load_config, read_jsonl


_PARSE_FAILURE_KEYS = ("null_or_empty", "malformed_json", "invalid_schema")
_DEFAULT_VALIDATION_THRESHOLDS = [round(i / 10, 1) for i in range(11)]


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
        "gold_coverage": 0.0,
        "mean_gold_confidence": 0.0,
        "macro_mean_gold_confidence": 0.0,
        "expected_unique_valid_rid_n": 0.0,
        "input_tokens": 0,
        "completion_tokens": 0,
        "expected_unique_valid_rids_per_1k_tokens": 0.0,
    }


def _gold_metrics(rows, gold_by_candidate):
    """Soft gold metrics for a list of candidate rows.

    Returns gold_audited_n, mean_gold_confidence, and expected_unique_valid_rid_n
    (sum of max confidence per unique rid — expected number of valid unique reactions).
    """
    confidences = []
    max_conf_by_rid = {}

    for row in rows:
        audit = gold_by_candidate.get(row.get("candidate_id"))
        if audit is None:
            continue
        confidence = audit.get("gold_confidence", 0.0)
        confidences.append(confidence)
        rid = row.get("rid")
        if rid is not None:
            if confidence > max_conf_by_rid.get(rid, -1.0):
                max_conf_by_rid[rid] = confidence

    return {
        "gold_audited_n": len(confidences),
        "mean_gold_confidence": _mean(confidences),
        "expected_unique_valid_rid_n": sum(max_conf_by_rid.values()),
    }


def _threshold_key(threshold):
    return f"{threshold:.1f}"


def _validation_confidence(row):
    return row.get("validation_confidence", row.get("confidence", 0.0)) or 0.0


def _validation_key(row):
    candidate_id = row.get("candidate_id") or row.get("sample_id")
    validation_model = row.get("validation_model") or row.get("source")
    return candidate_id, validation_model


def _latest_self_validations(rows):
    latest = {}
    for row in rows:
        key = _validation_key(row)
        if all(key):
            latest[key] = row
    return latest


def _gold_confidences(rows, gold_by_candidate):
    values = []
    for row in rows:
        audit = gold_by_candidate.get(row.get("candidate_id"))
        if audit is not None:
            values.append(audit.get("gold_confidence", 0.0))
    return values


def _self_validation_summary(
    candidates_by_model,
    validation_rows,
    gold_by_candidate,
    model_summaries,
    thresholds,
):
    latest_validations = _latest_self_validations(validation_rows)
    models = {}
    threshold_summary = {
        _threshold_key(threshold): {}
        for threshold in thresholds
    }

    for model, rows in candidates_by_model.items():
        validated = []
        for row in rows:
            validation = latest_validations.get((row.get("candidate_id"), model))
            if validation is not None:
                validated.append((row, validation))

        validation_confidences = [
            _validation_confidence(validation)
            for _, validation in validated
        ]
        gold_for_validated = _gold_confidences(
            [row for row, _ in validated],
            gold_by_candidate,
        )
        models[model] = {
            "candidate_n": len(rows),
            "validated_n": len(validated),
            "validation_coverage": _rate(len(validated), len(rows)),
            "mean_validation_confidence": _mean(validation_confidences),
            "gold_audited_validated_n": len(gold_for_validated),
            "gold_mean_for_validated": _mean(gold_for_validated),
        }

        total_gold_mass = sum(gold_for_validated)
        total_tokens = (
            model_summaries[model].get("input_tokens", 0)
            + model_summaries[model].get("completion_tokens", 0)
        )
        for threshold in thresholds:
            accepted = [
                row for row, validation in validated
                if _validation_confidence(validation) >= threshold
            ]
            rejected = [
                row for row, validation in validated
                if _validation_confidence(validation) < threshold
            ]
            accepted_gold = _gold_confidences(accepted, gold_by_candidate)
            rejected_gold = _gold_confidences(rejected, gold_by_candidate)
            gm = _gold_metrics(accepted, gold_by_candidate)
            accepted_rids = {
                row.get("rid")
                for row in accepted
                if row.get("rid") is not None
            }
            accepted_expected_valid = sum(accepted_gold)
            accepted_expected_invalid = sum(1.0 - value for value in accepted_gold)
            rejected_expected_valid = sum(rejected_gold)

            threshold_summary[_threshold_key(threshold)][model] = {
                "candidate_n": len(rows),
                "validated_n": len(validated),
                "accepted_candidate_n": len(accepted),
                "accepted_rate": _rate(len(accepted), len(validated)),
                "accepted_candidate_coverage": _rate(len(accepted), len(rows)),
                "accepted_unique_rid_n": len(accepted_rids),
                "accepted_parsed_rate": _mean(
                    [1.0 if row.get("parsed") else 0.0 for row in accepted]
                ),
                "accepted_balance_success_rate": _mean(
                    [1.0 if row.get("balance_success") else 0.0 for row in accepted]
                ),
                "accepted_gold_audited_n": len(accepted_gold),
                "accepted_mean_gold_confidence": _mean(accepted_gold),
                "accepted_expected_valid_candidate_n": accepted_expected_valid,
                "accepted_expected_unique_valid_rid_n": (
                    gm["expected_unique_valid_rid_n"]
                ),
                "accepted_expected_unique_valid_rids_per_1k_tokens": (
                    gm["expected_unique_valid_rid_n"] / total_tokens * 1000
                    if total_tokens else 0.0
                ),
                "precision_estimate": _mean(accepted_gold),
                "recall_estimate": _rate(accepted_expected_valid, total_gold_mass),
                "false_accept_gold_confidence_sum": accepted_expected_invalid,
                "rejected_valid_gold_confidence_sum": rejected_expected_valid,
            }

    return {
        "models": models,
        "thresholds": threshold_summary,
    }


def build_summary(*, data_dir, validation_thresholds=None):
    data_dir = Path(data_dir)
    generations = read_jsonl(data_dir / "generations.jsonl")
    candidates = read_jsonl(data_dir / "candidates.jsonl")
    gold_raw = read_jsonl(data_dir / "gold_audits.jsonl")
    validation_raw = read_jsonl(data_dir / "self_validation_runs.jsonl")
    thresholds = (
        validation_thresholds
        if validation_thresholds is not None
        else _DEFAULT_VALIDATION_THRESHOLDS
    )
    # gold_audits.jsonl is append-only; last entry per candidate_id wins.
    gold_by_candidate = {}
    for entry in gold_raw:
        if entry.get("candidate_id") is not None:
            gold_by_candidate[entry["candidate_id"]] = entry

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
        parse_denominators[model] += accepted + sum(
            source_stats.get(k, 0) for k in _PARSE_FAILURE_KEYS
        )

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
            parse_accepts[model], parse_denominators[model]
        )
        by_model[model]["parsed_rate"] = _rate(len(parsed), len(rows))
        by_model[model]["target_compound_rate"] = _rate(len(target_present), len(rows))
        by_model[model]["target_position_rate"] = _rate(len(target_correct), len(rows))
        by_model[model]["balance_success_rate"] = _rate(len(balanced), len(rows))
        by_model[model]["unique_rid_n"] = len(unique_rids)
        by_model[model]["unique_output_id_n"] = len(output_ids)
        by_model[model]["duplicate_rid_rate"] = 1 - _rate(len(unique_rids), len(rids))

        gm = _gold_metrics(rows, gold_by_candidate)
        by_model[model]["gold_audited_n"] = gm["gold_audited_n"]
        by_model[model]["gold_coverage"] = _rate(gm["gold_audited_n"], len(rows))
        by_model[model]["mean_gold_confidence"] = gm["mean_gold_confidence"]
        by_model[model]["expected_unique_valid_rid_n"] = gm["expected_unique_valid_rid_n"]
        total_tokens = by_model[model]["input_tokens"] + by_model[model]["completion_tokens"]
        by_model[model]["expected_unique_valid_rids_per_1k_tokens"] = (
            gm["expected_unique_valid_rid_n"] / total_tokens * 1000
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
            gm = _gold_metrics(rows, gold_by_candidate)
            bands[band_key][model] = {
                "candidate_n": len(rows),
                "parsed_rate": _mean([1.0 if row.get("parsed") else 0.0 for row in rows]),
                "target_position_rate": _mean(
                    [1.0 if row.get("target_position_correct") else 0.0 for row in rows]
                ),
                "balance_success_rate": _mean(
                    [1.0 if row.get("balance_success") else 0.0 for row in rows]
                ),
                **gm,
                "gold_coverage": _rate(gm["gold_audited_n"], len(rows)),
            }

    for model in by_model:
        audited_band_confidences = [
            bands[band_key][model]["mean_gold_confidence"]
            for band_key in bands
            if model in bands[band_key] and bands[band_key][model]["gold_audited_n"] > 0
        ]
        by_model[model]["macro_mean_gold_confidence"] = _mean(audited_band_confidences)

    return {
        "models": dict(by_model),
        "bands": bands,
        "self_validation": _self_validation_summary(
            candidates_by_model,
            validation_raw,
            gold_by_candidate,
            by_model,
            thresholds,
        ),
        "inputs": {
            "generations": len(generations),
            "candidates": len(candidates),
            "gold_audits": len(gold_raw),
            "self_validations": len(validation_raw),
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
    summary = build_summary(
        data_dir=args.data_dir,
        validation_thresholds=config.get("validation_thresholds"),
    )
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
