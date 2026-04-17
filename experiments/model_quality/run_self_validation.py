import json
from pathlib import Path

from scripts.chemy import Chemy
from scripts.infra.llm_client import LLMClient
from scripts.infra.parallel_runner import iter_parallel
from scripts.ops.raw_reactions.services import RawReactionValidationService

from .helpers import (
    NotebookLogger,
    append_jsonl,
    build_arg_parser,
    load_config,
    read_jsonl,
    resolve_path,
)


def _validation_key(entry):
    candidate_id = entry.get("candidate_id") or entry.get("sample_id")
    validation_model = entry.get("validation_model") or entry.get("source")
    return candidate_id, validation_model


def _existing_validations(path):
    """Return last validation row per (candidate_id, validation_model)."""
    result = {}
    for entry in read_jsonl(path):
        key = _validation_key(entry)
        if all(key):
            result[key] = entry
    return result


def _validation_job(job, llm_client, parser, validation_rounds):
    model = job["model"]
    entries = job["entries"]
    service = RawReactionValidationService(llm_client, parser)
    results = service.validate_entries(
        entries,
        model=model,
        max_rounds=validation_rounds,
        acceptance_threshold=0.0,
        early_stop=False,
        batch_size=len(entries) or 1,
    )
    if results is None:
        raise ValueError(f"Validation model {model!r} returned malformed output")

    enriched = []
    for candidate, result in zip(entries, results):
        result["candidate_id"] = candidate["candidate_id"]
        result["generator_model"] = candidate["model"]
        result["validation_model"] = model
        result["preset"] = candidate.get("preset")
        result["complexity_band"] = candidate.get("complexity_band")
        result["compound_name"] = candidate.get("compound_name")
        result["position"] = candidate.get("position")
        enriched.append(result)
    return enriched


def _batched(entries, batch_size):
    for i in range(0, len(entries), batch_size):
        yield entries[i:i + batch_size]


def run_self_validation(
    config,
    *,
    data_dir,
    root=None,
    llm_client=None,
    models=None,
    run_limit=None,
    max_workers=None,
    logger=None,
):
    data_dir = Path(data_dir)
    chemy_data_dir = resolve_path(config["chemy_data_dir"], root=root)
    chemy = Chemy(str(chemy_data_dir))
    logger = logger or getattr(chemy, "logger", None) or NotebookLogger()
    llm_client = llm_client or LLMClient(logger=logger)
    models = models or config.get("validation_models") or config["candidate_models"]
    models = set(models)
    max_workers = max_workers or config.get("max_workers", 1)
    validation_rounds = config.get("validation_rounds", 3)
    batch_size = config.get("validation_batch_size") or config.get("batch_size", 5)

    candidates_path = data_dir / "candidates.jsonl"
    validations_path = data_dir / "self_validation_runs.jsonl"
    manifests_dir = data_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    candidates = [
        entry for entry in read_jsonl(candidates_path)
        if entry.get("candidate_id") is not None and entry.get("model") in models
    ]
    existing = _existing_validations(validations_path)

    pending = []
    skipped = 0
    for candidate in candidates:
        key = (candidate["candidate_id"], candidate["model"])
        prior = existing.get(key)
        if prior and prior.get("validation_rounds", 0) >= validation_rounds:
            skipped += 1
            continue
        pending.append(candidate)
        if run_limit is not None and len(pending) >= run_limit:
            break

    grouped = {}
    for candidate in pending:
        grouped.setdefault(candidate["model"], []).append(candidate)

    jobs = []
    for model, entries in grouped.items():
        for batch in _batched(entries, batch_size):
            jobs.append({"model": model, "entries": batch})

    manifest = {
        "models": sorted(models),
        "candidate_count": len(candidates),
        "pending_count": len(pending),
        "validation_rounds": validation_rounds,
        "validation_batch_size": batch_size,
        "run_limit": run_limit,
        "max_workers": max_workers,
        "ran": 0,
        "skipped": skipped,
    }

    for result_batch in iter_parallel(
        llm_client,
        jobs,
        _validation_job,
        logger,
        max_workers=max_workers,
        routine_args=[llm_client, chemy.reaction_llm, validation_rounds],
        description="Running self-validation",
    ):
        for result in result_batch:
            append_jsonl(result, validations_path)
            manifest["ran"] += 1

    manifest_path = manifests_dir / "self_validation_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest


def main(argv=None):
    parser = build_arg_parser("Run self-validation for generated candidates.")
    parser.add_argument(
        "--model",
        action="append",
        dest="models",
        help="Candidate/generator model to self-validate. Repeatable; defaults to config.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of pending candidates to validate.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Maximum parallel validation workers.",
    )
    args = parser.parse_args(argv)
    config = load_config(args.config)
    manifest = run_self_validation(
        config,
        data_dir=args.data_dir,
        models=args.models,
        run_limit=args.limit,
        max_workers=args.workers,
    )
    print(
        f"Self-validation complete: ran={manifest['ran']}, "
        f"skipped={manifest['skipped']}, pending={manifest['pending_count']}"
    )


if __name__ == "__main__":
    main()
