import json
from pathlib import Path

from scripts.chemy import Chemy
from scripts.infra.llm_client import LLMClient
from scripts.infra.parallel_runner import iter_parallel
from scripts.ops.raw_reactions.records import ReactionTask
from scripts.ops.raw_reactions.services import (
    RawReactionEvaluator,
    RawReactionGenerationService,
    generation_run_id,
)

from .helpers import (
    NotebookLogger,
    append_jsonl,
    build_arg_parser,
    load_config,
    read_jsonl,
    resolve_path,
)


def _completed_generation_ids(path):
    return {
        entry["generation_id"]
        for entry in read_jsonl(path)
        if entry.get("generation_id") is not None
    }


def _completed_candidate_ids(path):
    return {
        entry["candidate_id"]
        for entry in read_jsonl(path)
        if entry.get("candidate_id") is not None
    }


def _generation_job(job, llm_client, parser, reactions, self_review):
    task = ReactionTask.from_dict(job["task"])
    model = job["model"]
    generation_id = job["generation_id"]
    generator = RawReactionGenerationService(llm_client)
    evaluator = RawReactionEvaluator(parser, reactions)

    generation = generator.generate(
        task,
        model,
        self_review=self_review,
    ).to_dict()
    generation["generation_id"] = generation_id
    generation["task"] = task.to_dict()
    candidates = []
    for candidate in evaluator.evaluate_generation(task, generation):
        candidate["generation_id"] = generation_id
        candidate["preset"] = task.preset
        candidate["complexity_band"] = task.complexity_band
        candidate["compound_name"] = task.compound_name
        candidate["position"] = task.position
        candidates.append(candidate)
    return {"generation": generation, "candidates": candidates}


def run_generation(
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
    llm_client = llm_client or LLMClient()
    models = models or config["candidate_models"]
    max_workers = max_workers or config.get("max_workers", 1)
    logger = logger or getattr(chemy, "logger", None) or NotebookLogger()

    tasks_path = data_dir / "tasks.jsonl"
    if not tasks_path.exists():
        raise FileNotFoundError(
            f"tasks.jsonl not found at {tasks_path}. "
            "Run `python -m experiments.model_quality.sample` first."
        )
    tasks = [ReactionTask.from_dict(entry) for entry in read_jsonl(tasks_path)]
    generations_path = data_dir / "generations.jsonl"
    candidates_path = data_dir / "candidates.jsonl"
    manifests_dir = data_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    existing_generations = {
        entry["generation_id"]: entry
        for entry in read_jsonl(generations_path)
        if entry.get("generation_id") is not None
    }
    completed_generations = set(existing_generations)
    completed_candidates = _completed_candidate_ids(candidates_path)
    evaluator = RawReactionEvaluator(chemy.reaction_llm, chemy.reactions)

    manifest = {
        "models": models,
        "task_n": len(tasks),
        "generation_self_review": config["generation_self_review"],
        "run_limit": run_limit,
        "max_workers": max_workers,
        "ran": 0,
        "skipped": 0,
        "candidate_n": 0,
    }

    jobs = []
    for model in models:
        for task in tasks:
            if run_limit is not None and len(jobs) >= run_limit:
                break
            generation_id = generation_run_id(task.task_id, model)
            if generation_id in completed_generations:
                generation = existing_generations[generation_id]
                for candidate in evaluator.evaluate_generation(task, generation):
                    if candidate["candidate_id"] in completed_candidates:
                        continue
                    candidate["generation_id"] = generation_id
                    candidate["preset"] = task.preset
                    candidate["complexity_band"] = task.complexity_band
                    candidate["compound_name"] = task.compound_name
                    candidate["position"] = task.position
                    append_jsonl(candidate, candidates_path)
                    completed_candidates.add(candidate["candidate_id"])
                    manifest["candidate_n"] += 1
                manifest["skipped"] += 1
                continue
            jobs.append({
                "task": task.to_dict(),
                "model": model,
                "generation_id": generation_id,
            })
        if run_limit is not None and len(jobs) >= run_limit:
            break

    for result in iter_parallel(
        llm_client,
        jobs,
        _generation_job,
        logger,
        max_workers=max_workers,
        routine_args=[
            llm_client,
            chemy.reaction_llm,
            chemy.reactions,
            config["generation_self_review"],
        ],
        description="Generating model-quality reactions",
    ):
        generation = result["generation"]
        append_jsonl(generation, generations_path)
        manifest["ran"] += 1
        for candidate in result["candidates"]:
            if candidate["candidate_id"] in completed_candidates:
                continue
            append_jsonl(candidate, candidates_path)
            completed_candidates.add(candidate["candidate_id"])
            manifest["candidate_n"] += 1

    manifest_path = manifests_dir / "generation_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest


def main(argv=None):
    parser = build_arg_parser("Run fresh-generation model bakeoff.")
    parser.add_argument(
        "--model",
        action="append",
        dest="models",
        help="Candidate model to run. Repeatable; defaults to config candidate_models.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of new generation runs to execute.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Maximum parallel generation workers.",
    )
    args = parser.parse_args(argv)
    config = load_config(args.config)
    manifest = run_generation(
        config,
        data_dir=args.data_dir,
        models=args.models,
        run_limit=args.limit,
        max_workers=args.workers,
    )
    print(
        f"Generation complete: ran={manifest['ran']}, "
        f"skipped={manifest['skipped']}, candidates={manifest['candidate_n']}"
    )


if __name__ == "__main__":
    main()
