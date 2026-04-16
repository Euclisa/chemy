import json
from pathlib import Path

from .helpers import build_arg_parser, load_config, read_jsonl


def _mean(values):
    return sum(values) / len(values) if values else 0.0


def summarize_validation(entries, *, acceptance_threshold):
    if not entries:
        return {
            "n": 0,
            "accepted_n": 0,
            "accepted_rate": 0.0,
            "mean_validation_confidence": 0.0,
            "mean_original_confidence": 0.0,
        }
    accepted = [
        entry
        for entry in entries
        if entry.get("validation_confidence", 0.0) >= acceptance_threshold
    ]
    return {
        "n": len(entries),
        "accepted_n": len(accepted),
        "accepted_rate": len(accepted) / len(entries),
        "mean_validation_confidence": _mean(
            [entry.get("validation_confidence", 0.0) for entry in entries]
        ),
        "mean_original_confidence": _mean(
            [
                entry.get("original_confidence", 0.0)
                for entry in entries
                if entry.get("original_confidence") is not None
            ]
        ),
    }


def summarize_repairs(entries):
    by_status = {}
    for entry in entries:
        status = entry.get("repair_status", "missing_status")
        by_status[status] = by_status.get(status, 0) + 1
    total = len(entries)
    return {
        "n": total,
        "status_counts": by_status,
        "ok_rate": by_status.get("ok", 0) / total if total else 0.0,
    }


def build_summary(config, *, data_dir):
    data_dir = Path(data_dir)
    acceptance_threshold = config["acceptance_threshold"]
    summary = {"audit": {}, "repairs": {}, "validation": {}}

    for path in sorted((data_dir / "audit").glob("*.jsonl")):
        summary["audit"][path.stem] = summarize_validation(
            read_jsonl(path),
            acceptance_threshold=acceptance_threshold,
        )

    for path in sorted((data_dir / "repairs").glob("*.jsonl")):
        summary["repairs"][path.stem] = summarize_repairs(read_jsonl(path))

    for path in sorted((data_dir / "validation").glob("*.jsonl")):
        summary["validation"][path.stem] = summarize_validation(
            read_jsonl(path),
            acceptance_threshold=acceptance_threshold,
        )

    return summary


def main(argv=None):
    parser = build_arg_parser("Summarize model-quality experiment outputs.")
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
