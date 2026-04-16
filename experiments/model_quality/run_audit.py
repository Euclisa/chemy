import json
from pathlib import Path

from scripts.infra.llm_client import LLMClient

from .helpers import (
    append_jsonl,
    build_arg_parser,
    completed_sample_ids,
    fixed_round_validate_entries,
    load_config,
    read_jsonl,
    usage_context,
)


def run_audit_group(
    llm_client,
    *,
    group_name,
    sample_path,
    output_path,
    model,
    rounds,
    batch_size,
):
    completed = completed_sample_ids(output_path)
    sample = [
        entry
        for entry in read_jsonl(sample_path)
        if entry.get("sample_id") not in completed
    ]
    if not sample:
        return {"group": group_name, "skipped": True, "pending": 0}

    with usage_context(llm_client) as usage:
        results = fixed_round_validate_entries(
            llm_client,
            sample,
            model=model,
            rounds=rounds,
            batch_size=batch_size,
        )

    for result in results:
        result["audit_group"] = group_name
        append_jsonl(result, output_path)

    return {
        "group": group_name,
        "skipped": False,
        "sample_path": str(sample_path),
        "output_path": str(output_path),
        "model": model,
        "rounds": rounds,
        "batch_size": batch_size,
        "completed_before": len(completed),
        "ran": len(results),
        "usage": usage,
    }


def run_audit(config, *, data_dir, llm_client=None, groups=None, model=None):
    data_dir = Path(data_dir)
    llm_client = llm_client or LLMClient()
    model = model or config["strong_model"]
    rounds = config["audit_rounds"]
    batch_size = config["batch_size"]
    samples_dir = data_dir / "samples"
    audit_dir = data_dir / "audit"
    manifests_dir = data_dir / "manifests"
    audit_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)

    if groups is None:
        groups = sorted(
            path.name[len("audit_"):-len(".jsonl")]
            for path in samples_dir.glob("audit_*.jsonl")
        )

    manifest = []
    for group_name in groups:
        sample_path = samples_dir / f"audit_{group_name}.jsonl"
        output_path = audit_dir / f"{group_name}_strong.jsonl"
        summary = run_audit_group(
            llm_client,
            group_name=group_name,
            sample_path=sample_path,
            output_path=output_path,
            model=model,
            rounds=rounds,
            batch_size=batch_size,
        )
        manifest.append(summary)

    manifest_path = manifests_dir / "audit_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest


def main(argv=None):
    parser = build_arg_parser("Run fixed-round strong-model validation audit.")
    parser.add_argument(
        "--group",
        action="append",
        dest="groups",
        help="Audit group to run, e.g. borderline_reject. Repeatable.",
    )
    parser.add_argument("--model", help="Override strong validation model.")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    manifest = run_audit(
        config,
        data_dir=args.data_dir,
        groups=args.groups,
        model=args.model,
    )
    print(
        "Audit complete: "
        + ", ".join(
            f"{entry['group']}={'skipped' if entry.get('skipped') else entry['ran']}"
            for entry in manifest
        )
    )


if __name__ == "__main__":
    main()
