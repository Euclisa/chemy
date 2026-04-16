import json
from pathlib import Path

from scripts.infra.llm_client import LLMClient

from .helpers import (
    append_jsonl,
    build_arg_parser,
    build_repair_candidates,
    completed_sample_ids,
    fixed_round_validate_entries,
    generate_repairs,
    load_config,
    read_jsonl,
    usage_context,
    validation_inputs_from_repairs,
)


def run_repair_generation(
    llm_client,
    *,
    sample_path,
    output_path,
    model,
    batch_size,
):
    completed = completed_sample_ids(output_path)
    sample = [
        entry
        for entry in read_jsonl(sample_path)
        if entry.get("sample_id") not in completed
    ]
    if not sample:
        return {"model": model, "skipped": True, "pending": 0}

    candidates = build_repair_candidates(sample)
    with usage_context(llm_client) as usage:
        results = generate_repairs(
            llm_client,
            candidates,
            model=model,
            batch_size=batch_size,
        )

    for result in results:
        append_jsonl(result, output_path)

    return {
        "model": model,
        "skipped": False,
        "sample_path": str(sample_path),
        "output_path": str(output_path),
        "completed_before": len(completed),
        "ran": len(results),
        "batch_size": batch_size,
        "usage": usage,
    }


def _original_validation_inputs(sample):
    return [
        {
            "sample_id": entry["sample_id"],
            "rid": entry.get("rid"),
            "cid": entry.get("cid"),
            "reaction": entry["reaction"],
            "original_confidence": entry.get("original_confidence", entry.get("confidence")),
            "original_rounds": entry.get("original_rounds", entry.get("rounds")),
            "original_positives": entry.get("original_positives", entry.get("positives")),
            "original_source": entry.get("original_source", entry.get("source")),
        }
        for entry in sample
    ]


def run_validation_route(
    llm_client,
    *,
    route,
    inputs,
    output_path,
    model,
    rounds,
    batch_size,
    acceptance_threshold,
):
    completed = completed_sample_ids(output_path)
    pending = [
        entry
        for entry in inputs
        if entry.get("sample_id") not in completed
    ]
    if not pending:
        return {"route": route, "model": model, "skipped": True, "pending": 0}

    with usage_context(llm_client) as usage:
        results = fixed_round_validate_entries(
            llm_client,
            pending,
            model=model,
            rounds=rounds,
            batch_size=batch_size,
        )

    by_sample_id = {entry["sample_id"]: entry for entry in pending}
    for result in results:
        source = by_sample_id[result["sample_id"]]
        result["route"] = route
        result["acceptance_threshold"] = acceptance_threshold
        result["accepted"] = result["validation_confidence"] >= acceptance_threshold
        if source.get("fixer_model") is not None:
            result["fixer_model"] = source["fixer_model"]
            result["fixed_reaction"] = source.get("fixed_reaction", source["reaction"])
            result["original_reaction"] = source.get("original_reaction")
        append_jsonl(result, output_path)

    return {
        "route": route,
        "model": model,
        "skipped": False,
        "output_path": str(output_path),
        "completed_before": len(completed),
        "ran": len(results),
        "rounds": rounds,
        "batch_size": batch_size,
        "usage": usage,
    }


def run_repairs(config, *, data_dir, llm_client=None, phase="all"):
    data_dir = Path(data_dir)
    llm_client = llm_client or LLMClient()
    cheap_model = config["cheap_model"]
    strong_model = config["strong_model"]
    batch_size = config["batch_size"]
    rounds = config["repair_validation_rounds"]
    acceptance_threshold = config["acceptance_threshold"]

    sample_path = data_dir / "samples" / "repair_borderline_reject.jsonl"
    repairs_dir = data_dir / "repairs"
    validation_dir = data_dir / "validation"
    manifests_dir = data_dir / "manifests"
    repairs_dir.mkdir(parents=True, exist_ok=True)
    validation_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)

    manifest = {"generation": [], "validation": []}

    repair_outputs = {
        "cheap_fixer": repairs_dir / "cheap_fixer.jsonl",
        "strong_fixer": repairs_dir / "strong_fixer.jsonl",
    }

    if phase in {"all", "generate"}:
        manifest["generation"].append(
            run_repair_generation(
                llm_client,
                sample_path=sample_path,
                output_path=repair_outputs["cheap_fixer"],
                model=cheap_model,
                batch_size=batch_size,
            )
        )
        manifest["generation"].append(
            run_repair_generation(
                llm_client,
                sample_path=sample_path,
                output_path=repair_outputs["strong_fixer"],
                model=strong_model,
                batch_size=batch_size,
            )
        )

    if phase in {"all", "validate"}:
        original_inputs = _original_validation_inputs(read_jsonl(sample_path))
        cheap_fix_inputs = validation_inputs_from_repairs(
            read_jsonl(repair_outputs["cheap_fixer"])
        )
        strong_fix_inputs = validation_inputs_from_repairs(
            read_jsonl(repair_outputs["strong_fixer"])
        )

        routes = [
            (
                "original_strong_validator",
                original_inputs,
                strong_model,
                validation_dir / "original_strong_validator.jsonl",
            ),
            (
                "cheap_fix_cheap_validator",
                cheap_fix_inputs,
                cheap_model,
                validation_dir / "cheap_fix_cheap_validator.jsonl",
            ),
            (
                "cheap_fix_strong_validator",
                cheap_fix_inputs,
                strong_model,
                validation_dir / "cheap_fix_strong_validator.jsonl",
            ),
            (
                "strong_fix_strong_validator",
                strong_fix_inputs,
                strong_model,
                validation_dir / "strong_fix_strong_validator.jsonl",
            ),
        ]

        for route, inputs, model, output_path in routes:
            manifest["validation"].append(
                run_validation_route(
                    llm_client,
                    route=route,
                    inputs=inputs,
                    output_path=output_path,
                    model=model,
                    rounds=rounds,
                    batch_size=batch_size,
                    acceptance_threshold=acceptance_threshold,
                )
            )

    manifest_path = manifests_dir / "repair_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest


def main(argv=None):
    parser = build_arg_parser("Run model-quality repair routing experiment.")
    parser.add_argument(
        "--phase",
        choices=["all", "generate", "validate"],
        default="all",
        help="Run repair generation, validation, or both.",
    )
    args = parser.parse_args(argv)

    config = load_config(args.config)
    manifest = run_repairs(config, data_dir=args.data_dir, phase=args.phase)
    generation = ", ".join(
        f"{entry['model']}={'skipped' if entry.get('skipped') else entry['ran']}"
        for entry in manifest["generation"]
    )
    validation = ", ".join(
        f"{entry['route']}={'skipped' if entry.get('skipped') else entry['ran']}"
        for entry in manifest["validation"]
    )
    print(f"Repair generation: {generation or 'not run'}")
    print(f"Repair validation: {validation or 'not run'}")


if __name__ == "__main__":
    main()
