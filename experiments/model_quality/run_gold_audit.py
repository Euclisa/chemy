import json
from collections import defaultdict
from pathlib import Path

from scripts.infra.llm_client import LLMClient
from scripts.infra.parallel_runner import iter_parallel
from scripts.ops.raw_reactions.prompts import format_structured_reaction
from scripts.ops.raw_reactions.services import aggregate_gold_votes, batched

from .helpers import NotebookLogger, append_jsonl, build_arg_parser, load_config, read_jsonl


GOLD_AUDIT_INSTRUCT = (
    "You are a chemistry database auditor. "
    "You will receive a numbered list of chemical reaction candidates.\n\n"
    "For each candidate, output exactly one JSON object on its own line:\n"
    '{"index": <int>, "valid": <bool>}\n\n'
    "Field definitions:\n"
    "- index: echo the candidate's input number exactly\n"
    "- valid: true if this is a real, documented chemical reaction; "
    "false if it is fabricated, chemically impossible, or not reliably documented\n\n"
    "Output one JSON line per candidate, in input order. "
    "No markdown, no code fences, no extra text."
)


def _format_candidate(candidate, index):
    reaction = format_structured_reaction(candidate["reaction"])
    return (
        f"{index}. target={candidate.get('compound_name')} "
        f"position={candidate.get('position')} reaction={reaction}"
    )


def _parse_gold_response(response, expected_count):
    """Parse JSONL audit response into a list of {index, valid} dicts.

    Rows are sorted by their index field so a model that reorders output is
    handled correctly. Returns None on any parse or structural error.
    """
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
    try:
        indexed = [(int(row["index"]), bool(row["valid"])) for row in rows]
    except (KeyError, TypeError, ValueError):
        return None
    if sorted(i for i, _ in indexed) != list(range(1, expected_count + 1)):
        return None
    indexed.sort(key=lambda t: t[0])
    return [{"index": i, "valid": v} for i, v in indexed]


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
    # parsed is sorted by index, which matches the original batch order (1-based).
    return {
        "round": job["round"],
        "votes": [
            {"candidate_id": candidate["candidate_id"], "vote": vote}
            for candidate, vote in zip(batch, parsed)
        ],
    }


def _existing_audits(path):
    """Return {candidate_id: entry} from gold_audits.jsonl (last entry wins)."""
    result = {}
    for entry in read_jsonl(path):
        cid = entry.get("candidate_id")
        if cid is not None:
            result[cid] = entry
    return result


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
    # Resolve logger first so it can be wired into LLMClient for error reporting.
    logger = logger or NotebookLogger()
    llm_client = llm_client or LLMClient(logger=logger)
    model = model or config["strong_model"]
    gold_rounds = config["gold_rounds"]
    max_workers = max_workers or config.get("max_workers", 1)
    candidates_path = data_dir / "candidates.jsonl"
    output_path = data_dir / "gold_audits.jsonl"
    manifests_dir = data_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    # Candidates whose existing audit already has gold_rounds votes are done.
    # Candidates with fewer votes (or none) are re-audited from scratch so that
    # increasing gold_rounds triggers fresh complete runs.
    existing = _existing_audits(output_path)
    completed = {
        cid for cid, audit in existing.items()
        if len(audit.get("gold_votes") or []) >= gold_rounds
    }

    all_candidates = read_jsonl(candidates_path)
    candidates = [
        entry for entry in all_candidates
        if entry.get("candidate_id") not in completed
    ]
    by_id = {c["candidate_id"]: c for c in candidates}

    jobs = []
    for round_i in range(gold_rounds):
        for batch in batched(candidates, config["batch_size"]):
            jobs.append({"round": round_i, "batch": batch})

    votes_by_candidate = defaultdict(list)
    written_this_run = set()
    batches_received = 0

    def _flush_candidate(candidate_id):
        candidate = by_id[candidate_id]
        votes = votes_by_candidate[candidate_id]
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
        written_this_run.add(candidate_id)

    for result in iter_parallel(
        llm_client,
        jobs,
        _audit_batch_job,
        logger,
        max_workers=max_workers,
        routine_args=[llm_client, model],
        description="Running gold audit batches",
    ):
        batches_received += 1
        for item in result["votes"]:
            cid = item["candidate_id"]
            votes_by_candidate[cid].append(item["vote"])
            # Write as soon as this candidate has all its votes — crash-safe.
            if cid not in written_this_run and len(votes_by_candidate[cid]) >= gold_rounds:
                _flush_candidate(cid)

    # Flush any candidates that received some votes but never hit gold_rounds
    # (e.g. a partial run that was interrupted before the final batch came in).
    for cid in votes_by_candidate:
        if cid not in written_this_run:
            _flush_candidate(cid)

    batch_failures = len(jobs) - batches_received
    if batch_failures:
        logger.log(
            f"WARNING: {batch_failures}/{len(jobs)} audit batches failed — "
            "check ERROR lines above for details. "
            "Re-run gold audit to retry failed batches."
        )

    manifest = {
        "model": model,
        "rounds": gold_rounds,
        "max_workers": max_workers,
        "candidate_count": len(all_candidates),
        "pending_count": len(candidates),
        "completed_before": len(completed),
        "audited": len(written_this_run),
        "batch_failures": batch_failures,
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
        f"completed_before={manifest['completed_before']}, "
        f"batch_failures={manifest['batch_failures']}"
    )


if __name__ == "__main__":
    main()
