import argparse
import hashlib
import json
import os
import random
from contextlib import nullcontext
from copy import deepcopy
from pathlib import Path

from scripts.infra.store import JsonlStore
from scripts.ops.raw_reactions.prompts import (
    REPAIR_BATCH_INSTRUCT,
    VALIDATE_BATCH_INSTRUCT,
    format_reactions_for_validation,
    format_repair_candidates,
)
from scripts.ops.raw_reactions.schema import is_valid_reaction_obj


EXPERIMENT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = EXPERIMENT_DIR / "config.yaml"
DEFAULT_DATA_DIR = EXPERIMENT_DIR / "data"

DEFAULT_CONFIG = {
    "seed": 1729,
    "acceptance_threshold": 0.4,
    "audit_rounds": 9,
    "repair_validation_rounds": 9,
    "batch_size": 5,
    "cheap_model": "openai/gpt-oss-120b",
    "strong_model": "google/gemini-2.5-pro",
    "verdict_path": "data/raw_reactions/verdict",
    "samples": {
        "high_accept": {"confidence_min": 0.8, "confidence_max": None, "n": 80},
        "borderline_accept": {"confidence_min": 0.4, "confidence_max": 0.6, "n": 80},
        "borderline_reject": {"confidence_min": 0.2, "confidence_max": 0.4, "n": 80},
        "low_reject": {"confidence_min": 0.0, "confidence_max": 0.1, "n": 40},
    },
    "repair": {"confidence_min": 0.2, "confidence_max": 0.4, "n": 80},
}


def _parse_scalar(value):
    value = value.strip()
    if value in {"", "null", "None", "~"}:
        return None
    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if (
        (value.startswith('"') and value.endswith('"'))
        or (value.startswith("'") and value.endswith("'"))
    ):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _load_simple_yaml(path):
    """Load the small YAML subset used by config.yaml without adding PyYAML."""
    root = {}
    stack = [(-1, root)]
    with open(path) as f:
        for raw_line in f:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            indent = len(raw_line) - len(raw_line.lstrip(" "))
            line = raw_line.strip()
            if ":" not in line:
                raise ValueError(f"Unsupported config line: {raw_line.rstrip()}")
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()

            while indent <= stack[-1][0]:
                stack.pop()

            parent = stack[-1][1]
            if value == "":
                node = {}
                parent[key] = node
                stack.append((indent, node))
            else:
                parent[key] = _parse_scalar(value)
    return root


def deep_merge(base, override):
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(path=DEFAULT_CONFIG_PATH):
    if path is None:
        return deepcopy(DEFAULT_CONFIG)
    path = Path(path)
    if not path.exists():
        return deepcopy(DEFAULT_CONFIG)
    loaded = _load_simple_yaml(path)
    return deep_merge(DEFAULT_CONFIG, loaded)


def resolve_path(path, *, root=None):
    path = Path(path)
    if path.is_absolute():
        return path
    root = Path.cwd() if root is None else Path(root)
    return root / path


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_jsonl(path):
    store = JsonlStore()
    path = str(path)
    if os.path.isdir(path):
        store.register_vault(path, "verdict_")
    return store.load_jsonl(path)


def write_jsonl(entries, path):
    ensure_parent(path)
    with open(path, "w") as f:
        for entry in entries:
            f.write(json.dumps(entry) + "\n")


def append_jsonl(entry, path):
    ensure_parent(path)
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def file_sha256(path):
    path = Path(path)
    h = hashlib.sha256()
    if path.is_dir():
        for child in sorted(path.glob("*.jsonl")):
            h.update(child.name.encode("utf-8"))
            with open(child, "rb") as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(chunk)
    else:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
    return h.hexdigest()


def canonical_hash(value, prefix):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return f"{prefix}:{hashlib.sha1(payload.encode('utf-8')).hexdigest()}"


def get_sample_id(entry):
    rid = entry.get("rid")
    if rid is not None:
        return str(rid)
    reaction = entry.get("reaction")
    if reaction is not None:
        return canonical_hash(reaction, "reaction")
    return canonical_hash(entry, "entry")


def normalize_verdict_entry(entry):
    normalized = entry.copy()
    normalized["sample_id"] = get_sample_id(entry)
    normalized["original_confidence"] = entry.get("confidence", 0.0)
    normalized["original_rounds"] = entry.get("rounds")
    normalized["original_positives"] = entry.get("positives")
    normalized["original_source"] = entry.get("source")
    return normalized


def dedupe_by_sample_id(entries):
    grouped = {}
    for entry in entries:
        normalized = normalize_verdict_entry(entry)
        sample_id = normalized["sample_id"]
        current = grouped.get(sample_id)
        if current is None:
            grouped[sample_id] = normalized
            continue
        current_key = (
            current.get("confidence", 0.0) or 0.0,
            current.get("positives", 0) or 0,
            current.get("rounds", 0) or 0,
        )
        new_key = (
            normalized.get("confidence", 0.0) or 0.0,
            normalized.get("positives", 0) or 0,
            normalized.get("rounds", 0) or 0,
        )
        if new_key > current_key:
            grouped[sample_id] = normalized
    return list(grouped.values())


def filter_by_confidence(entries, confidence_min, confidence_max):
    selected = []
    for entry in entries:
        confidence = entry.get("confidence")
        if confidence is None:
            continue
        if confidence < confidence_min:
            continue
        if confidence_max is not None and confidence >= confidence_max:
            continue
        selected.append(entry)
    return selected


def deterministic_sample(entries, n, seed, salt):
    entries = sorted(entries, key=lambda entry: entry["sample_id"])
    rng = random.Random(f"{seed}:{salt}")
    if n >= len(entries):
        return list(entries)
    return rng.sample(entries, n)


def completed_sample_ids(path):
    return {
        entry["sample_id"]
        for entry in read_jsonl(path)
        if entry.get("sample_id") is not None
    }


def parse_indexed_verdicts(response, expected_count):
    verdicts = []
    for raw_line in response.splitlines():
        line = raw_line.strip().lower()
        if not line:
            continue
        if "invalid" in line:
            verdicts.append(False)
        elif "valid" in line:
            verdicts.append(True)
    if len(verdicts) != expected_count:
        return None
    return verdicts


def parse_repair_response(response, expected_count):
    repairs = []
    for raw_line in response.strip().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.lower() == "null":
            repairs.append(None)
            continue
        try:
            parsed = json.loads(line)
        except (json.JSONDecodeError, TypeError):
            return None
        if not isinstance(parsed, dict):
            return None
        repairs.append(parsed)
    if len(repairs) != expected_count:
        return None
    return repairs


def batched(entries, batch_size):
    for idx in range(0, len(entries), batch_size):
        yield entries[idx:idx + batch_size]


def usage_context(llm_client):
    tracker = getattr(llm_client, "track_usage", None)
    if tracker is None:
        return nullcontext({})
    return tracker()


def fixed_round_validate_entries(
    llm_client,
    entries,
    *,
    model,
    rounds,
    batch_size=5,
    reaction_field="reaction",
):
    votes_by_id = {entry["sample_id"]: [] for entry in entries}
    for _ in range(rounds):
        for batch in batched(entries, batch_size):
            staged = [
                {"cid": entry.get("cid"), "reaction": entry[reaction_field]}
                for entry in batch
            ]
            prompt = f"{VALIDATE_BATCH_INSTRUCT}\n{format_reactions_for_validation(staged)}"
            response = llm_client.fetch_answer_str(prompt, model)
            verdicts = parse_indexed_verdicts(response, len(batch))
            if verdicts is None:
                raise ValueError(
                    f"Model {model!r} returned an invalid verdict count "
                    f"for batch of {len(batch)}"
                )
            for entry, verdict in zip(batch, verdicts):
                votes_by_id[entry["sample_id"]].append(verdict)

    results = []
    for entry in entries:
        votes = votes_by_id[entry["sample_id"]]
        positives = sum(1 for vote in votes if vote)
        confidence = positives / len(votes) if votes else 0.0
        result = {
            "sample_id": entry["sample_id"],
            "rid": entry.get("rid"),
            "cid": entry.get("cid"),
            "reaction": entry[reaction_field],
            "original_confidence": entry.get("original_confidence", entry.get("confidence")),
            "original_rounds": entry.get("original_rounds", entry.get("rounds")),
            "original_positives": entry.get("original_positives", entry.get("positives")),
            "original_source": entry.get("original_source", entry.get("source")),
            "validation_model": model,
            "validation_rounds": len(votes),
            "validation_positives": positives,
            "validation_confidence": confidence,
            "validation_votes": votes,
        }
        results.append(result)
    return results


def build_repair_candidates(entries):
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
        for entry in entries
    ]


def generate_repairs(llm_client, candidates, *, model, batch_size=5):
    results = []
    for batch in batched(candidates, batch_size):
        prompt = f"{REPAIR_BATCH_INSTRUCT}\n\n{format_repair_candidates(batch)}"
        response = llm_client.fetch_answer_str(prompt, model)
        repairs = parse_repair_response(response, len(batch))
        if repairs is None:
            raise ValueError(
                f"Model {model!r} returned invalid repair output "
                f"for batch of {len(batch)}"
            )
        for candidate, repair in zip(batch, repairs):
            if repair is None:
                status = "null_repair"
            elif not is_valid_reaction_obj(repair):
                status = "invalid_schema"
            else:
                status = "ok"
            results.append(
                {
                    "sample_id": candidate["sample_id"],
                    "rid": candidate.get("rid"),
                    "cid": candidate.get("cid"),
                    "original_reaction": candidate["reaction"],
                    "original_confidence": candidate.get("original_confidence"),
                    "original_rounds": candidate.get("original_rounds"),
                    "original_positives": candidate.get("original_positives"),
                    "original_source": candidate.get("original_source"),
                    "fixer_model": model,
                    "fixed_reaction": repair,
                    "repair_status": status,
                }
            )
    return results


def validation_inputs_from_repairs(repairs):
    inputs = []
    for repair in repairs:
        if repair.get("repair_status") != "ok":
            continue
        inputs.append(
            {
                "sample_id": repair["sample_id"],
                "rid": repair.get("rid"),
                "cid": repair.get("cid"),
                "reaction": repair["fixed_reaction"],
                "fixed_reaction": repair["fixed_reaction"],
                "original_reaction": repair.get("original_reaction"),
                "original_confidence": repair.get("original_confidence"),
                "original_rounds": repair.get("original_rounds"),
                "original_positives": repair.get("original_positives"),
                "original_source": repair.get("original_source"),
                "fixer_model": repair.get("fixer_model"),
            }
        )
    return inputs


def add_common_args(parser):
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="Path to model-quality experiment config.",
    )
    parser.add_argument(
        "--data-dir",
        default=str(DEFAULT_DATA_DIR),
        help="Directory for frozen samples and experiment outputs.",
    )
    return parser


def build_arg_parser(description):
    return add_common_args(argparse.ArgumentParser(description=description))
