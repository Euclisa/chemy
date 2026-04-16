import json
from collections import defaultdict
from pathlib import Path

from scripts.infra.llm_client import LLMClient
from scripts.infra.parallel_runner import iter_parallel
from scripts.ops.raw_reactions.prompts import format_structured_reaction
from scripts.ops.raw_reactions.services import aggregate_gold_votes, batched

from .helpers import append_jsonl, build_arg_parser, load_config, read_jsonl


GOLD_AUDIT_INSTRUCT = (
    "You are auditing generated chemical reaction rows for a chemistry database. "
    "For each indexed candidate, judge whether the chemistry is documented or standard, "
    "whether the named target compound appears in the requested position, whether phases "
    "are plausible, and whether solvent/catalyst fields are correctly used.\n"
    "Return exactly one JSON object per input line, in the same order. Fields: "
    "index, valid_chemistry, documented_or_standard, target_compound_present, "
    "target_position_correct, phase_ok, roles_ok, solvent_catalyst_ok, "
    "error_category, confidence. Use booleans for *_ok and validity fields. "
    "Use error_category=null when the row is acceptable; otherwise use one of "
    "wrong_chemistry, undocumented, wrong_target, wrong_position, bad_phase, "
    "bad_roles, bad_solvent_catalyst, or unclear. Return ONLY JSONL."
)


class _NotebookLogger:
    def log(self, message):
        print(message)

    def track(self, iterable, *args, **kwargs):
        return iterable


def _format_candidate(candidate, index):
    reaction = format_structured_reaction(candidate["reaction"])
    return (
        f"{index}. target={candidate.get('compound_name')} "
        f"position={candidate.get('position')} reaction={reaction}"
    )


def _parse_gold_response(response, expected_count):
    rows = []
    for raw_line in (response or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            return None
        if not isinstance(row, dict):
            return None
        rows.append(row)
    if len(rows) != expected_count:
        return None
    return rows


def _audit_batch_job(job, llm_client, model):
    batch = job["batch"]
    body = "\n".join(
        _format_candidate(candidate, idx + 1)
        for idx, candidate in enumerate(batch)
    )
    response = llm_client.fetch_answer_str(f"{GOLD_AUDIT_INSTRUCT}\n\n{body}", model)
    parsed = _parse_gold_response(response, len(batch))
    if parsed is None:
        raise ValueError(f"Gold model {model!r} returned malformed audit output")
    return {
        "round": job["round"],
        "votes": [
            {"candidate_id": candidate["candidate_id"], "vote": vote}
            for candidate, vote in zip(batch, parsed)
        ],
    }


def _completed_candidate_ids(path):
    return {
        entry["candidate_id"]
        for entry in read_jsonl(path)
        if entry.get("candidate_id") is not None
    }


def run_gold_audit(
    config,
    *,
    data_dir,
    llm_client=None,
    model=None,
    max_workers=None,
    logger=None,
):
    data_dir = Path(data_dir)
    llm_client = llm_client or LLMClient()
    model = model or config["strong_model"]
    max_workers = max_workers or config.get("max_workers", 1)
    logger = logger or _NotebookLogger()
    candidates_path = data_dir / "candidates.jsonl"
    output_path = data_dir / "gold_audits.jsonl"
    manifests_dir = data_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    completed = _completed_candidate_ids(output_path)
    candidates = [
        entry
        for entry in read_jsonl(candidates_path)
        if entry.get("candidate_id") not in completed
    ]
    votes_by_candidate = defaultdict(list)
    jobs = []
    for round_i in range(config["gold_rounds"]):
        for batch in batched(candidates, config["batch_size"]):
            jobs.append({"round": round_i, "batch": batch})

    for result in iter_parallel(
        llm_client,
        jobs,
        _audit_batch_job,
        logger,
        max_workers=max_workers,
        routine_args=[llm_client, model],
        description="Running gold audit batches",
    ):
        for item in result["votes"]:
            votes_by_candidate[item["candidate_id"]].append(item["vote"])

    by_id = {candidate["candidate_id"]: candidate for candidate in candidates}
    for candidate_id, votes in votes_by_candidate.items():
        candidate = by_id[candidate_id]
        aggregate = aggregate_gold_votes(votes)
        output = {
            "candidate_id": candidate_id,
            "task_id": candidate.get("task_id"),
            "model": candidate.get("model"),
            "rid": candidate.get("rid"),
            "output_id": candidate.get("output_id"),
            "gold_model": model,
            "gold_votes": votes,
            **aggregate,
        }
        append_jsonl(output, output_path)

    manifest = {
        "model": model,
        "rounds": config["gold_rounds"],
        "max_workers": max_workers,
        "completed_before": len(completed),
        "audited": len(votes_by_candidate),
    }
    manifest_path = manifests_dir / "gold_audit_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest


def main(argv=None):
    parser = build_arg_parser("Run strong-model gold audit for generated candidates.")
    parser.add_argument("--model", help="Override configured strong gold model.")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Maximum parallel gold-audit workers.",
    )
    args = parser.parse_args(argv)
    config = load_config(args.config)
    manifest = run_gold_audit(
        config,
        data_dir=args.data_dir,
        model=args.model,
        max_workers=args.workers,
    )
    print(
        f"Gold audit complete: audited={manifest['audited']}, "
        f"completed_before={manifest['completed_before']}"
    )


if __name__ == "__main__":
    main()
