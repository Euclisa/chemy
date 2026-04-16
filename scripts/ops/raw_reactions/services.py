import hashlib
import json
import time
from collections import Counter

from .presets import POSITION_ANY, POSITION_PRODUCT, POSITION_REAGENT
from .prompts import (
    REPAIR_BATCH_INSTRUCT,
    VALIDATE_BATCH_INSTRUCT,
    build_fetch_prompt,
    build_revalidate_prompt,
    format_reactions_for_validation,
    format_repair_candidates,
)
from .records import GenerationResult, ReactionCandidate, ReactionTask, ValidationResult
from .schema import is_valid_reaction_obj, normalize_reaction_obj


def canonical_payload(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def stable_hash(value, prefix):
    digest = hashlib.sha1(canonical_payload(value).encode("utf-8")).hexdigest()
    return f"{prefix}:{digest}"


def reaction_output_id(reaction):
    return stable_hash(reaction, "reaction")


def generation_run_id(task_id, model):
    return stable_hash({"task_id": task_id, "model": model}, "generation")


def candidate_id(task_id, model, reaction_index, output_id):
    return stable_hash(
        {
            "task_id": task_id,
            "model": model,
            "reaction_index": reaction_index,
            "output_id": output_id,
        },
        "candidate",
    )


def _count(stats, key):
    if stats is not None:
        stats[key] = stats.get(key, 0) + 1


def parse_reactions_jsonl(response, stats=None):
    reactions = []
    for line in (response or "").strip().split("\n"):
        line = line.strip()
        if not line or line.lower() == "null":
            _count(stats, "null_or_empty")
            continue
        try:
            obj = json.loads(line)
        except (json.JSONDecodeError, TypeError):
            _count(stats, "malformed_json")
            continue
        if not is_valid_reaction_obj(obj):
            _count(stats, "invalid_schema")
            continue
        reactions.append(normalize_reaction_obj(obj))
        _count(stats, "accepted")
    return reactions


def parse_indexed_verdicts(response, expected_count):
    verdicts = []
    for raw_line in (response or "").splitlines():
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
    for raw_line in (response or "").strip().splitlines():
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


class LLMCallRunner:
    def __init__(self, llm_client):
        self.llm_client = llm_client

    def call(self, prompt, model, *, reasoning_effort="medium", max_completion_tokens=10000):
        start_input = getattr(self.llm_client, "input_tokens_total", 0)
        start_completion = getattr(self.llm_client, "completion_tokens_total", 0)
        start = time.perf_counter()
        try:
            response = self.llm_client.fetch_answer_str(
                prompt,
                model,
                reasoning_effort=reasoning_effort,
                max_completion_tokens=max_completion_tokens,
            )
            error = None
        except getattr(type(self.llm_client), "CompletionTokensLimitReached", ()):
            raise
        except Exception as exc:
            response = None
            error = str(exc)
        latency_s = time.perf_counter() - start
        usage = {
            "input_tokens": getattr(self.llm_client, "input_tokens_total", 0) - start_input,
            "completion_tokens": (
                getattr(self.llm_client, "completion_tokens_total", 0) - start_completion
            ),
        }
        return {
            "model": model,
            "raw_response": response,
            "usage": usage,
            "latency_s": latency_s,
            "error": error,
        }


class RawReactionTaskBuilder:
    MAX_EXISTING_CONTEXT = 20

    def __init__(self, reactions=None):
        self.reactions = reactions
        self._existing_by_cid = {}

    def set_existing_reactions(self, existing_by_cid):
        self._existing_by_cid = existing_by_cid or {}

    def get_existing_context(self, cid, position):
        by_position = self._existing_by_cid.get(cid, {})
        if position == POSITION_REAGENT:
            reactions = by_position.get(POSITION_REAGENT, [])
        elif position == POSITION_PRODUCT:
            reactions = by_position.get(POSITION_PRODUCT, [])
        elif position == POSITION_ANY:
            reactions = []
            seen = set()
            for side in (POSITION_REAGENT, POSITION_PRODUCT):
                for reaction in by_position.get(side, []):
                    if reaction in seen:
                        continue
                    seen.add(reaction)
                    reactions.append(reaction)
        else:
            raise ValueError(f"Unknown position: {position!r}")
        return reactions[:self.MAX_EXISTING_CONTEXT]

    def build_existing_reactions_context(self, parsed_reactions):
        by_cid = {}
        if self.reactions is None:
            return by_cid

        for react in parsed_reactions:
            try:
                reaction_str = self.reactions.get_reaction_as_str(react)
            except (KeyError, TypeError):
                continue
            for entry in react["reagents"]:
                by_cid.setdefault(entry["cid"], {}).setdefault(
                    POSITION_REAGENT, []
                ).append(reaction_str)
            for entry in react["products"]:
                by_cid.setdefault(entry["cid"], {}).setdefault(
                    POSITION_PRODUCT, []
                ).append(reaction_str)
        return by_cid

    def build_task(self, chem, preset, *, complexity_band=None, mode="cold_start"):
        existing = tuple(self.get_existing_context(chem["cid"], preset.position))
        prompt = build_fetch_prompt(
            chem["cmpdname"],
            preset.position,
            preset.scope,
            existing,
        )
        task_payload = {
            "cid": chem["cid"],
            "compound_name": chem["cmpdname"],
            "preset": preset.name,
            "position": preset.position,
            "scope": preset.scope,
            "complexity_band": complexity_band,
            "mode": mode,
            "existing_reactions": list(existing),
        }
        return ReactionTask(
            task_id=stable_hash(task_payload, "task"),
            cid=chem["cid"],
            compound_name=chem["cmpdname"],
            preset=preset.name,
            position=preset.position,
            scope=preset.scope,
            prompt=prompt,
            existing_reactions=existing,
            complexity_band=complexity_band,
            mode=mode,
            metadata={
                "bertz_complexity": chem.get("bertz_complexity"),
                "heavy_count": chem.get("heavy_count"),
                "organic": chem.get("organic"),
                "wiki": chem.get("wiki"),
            },
        )


class RawReactionGenerationService:
    def __init__(self, llm_client):
        self.runner = LLMCallRunner(llm_client)

    def generate(self, task, model, *, self_review=True):
        task = task if isinstance(task, ReactionTask) else ReactionTask.from_dict(task)
        call = self.runner.call(task.prompt, model)
        if call["error"] is not None:
            return GenerationResult(
                task_id=task.task_id,
                model=model,
                raw_response=call["raw_response"],
                reactions=(),
                parse_stats={},
                usage=call["usage"],
                latency_s=call["latency_s"],
                error=call["error"],
            )

        stats = {}
        reactions = parse_reactions_jsonl(call["raw_response"], stats)
        reviewed_raw = None
        reviewed = []
        reviewed_stats = {}
        total_usage = call["usage"].copy()
        total_latency = call["latency_s"]

        if self_review and reactions:
            reactions_jsonl_str = "\n".join(json.dumps(r) for r in reactions)
            review_prompt = build_revalidate_prompt(reactions_jsonl_str)
            review_call = self.runner.call(review_prompt, model)
            reviewed_raw = review_call["raw_response"]
            total_latency += review_call["latency_s"]
            total_usage["input_tokens"] += review_call["usage"]["input_tokens"]
            total_usage["completion_tokens"] += review_call["usage"]["completion_tokens"]
            if review_call["error"] is None:
                reviewed = parse_reactions_jsonl(reviewed_raw, reviewed_stats)
                if reviewed:
                    reactions = reviewed

        return GenerationResult(
            task_id=task.task_id,
            model=model,
            raw_response=call["raw_response"],
            reactions=tuple(reactions),
            parse_stats=stats,
            reviewed_raw_response=reviewed_raw,
            reviewed_reactions=tuple(reviewed),
            reviewed_parse_stats=reviewed_stats,
            usage=total_usage,
            latency_s=total_latency,
        )


class RawReactionEvaluator:
    def __init__(self, parser, reactions=None):
        self.parser = parser
        self.reactions = reactions or getattr(parser, "reactions", None)

    def _target_flags(self, parsed, task):
        if not parsed:
            return False, False
        cid = task.cid if isinstance(task, ReactionTask) else task["cid"]
        position = task.position if isinstance(task, ReactionTask) else task["position"]
        reagent_cids = {entry["cid"] for entry in parsed["reagents"]}
        product_cids = {entry["cid"] for entry in parsed["products"]}
        present = cid in reagent_cids or cid in product_cids
        if position == POSITION_REAGENT:
            correct = cid in reagent_cids
        elif position == POSITION_PRODUCT:
            correct = cid in product_cids
        else:
            correct = present
        return present, correct

    def _balance_success(self, parsed):
        if not parsed or self.reactions is None:
            return False
        try:
            return self.reactions.balance_reaction(parsed) is not None
        except Exception:
            return False

    def evaluate_reaction(self, task, model, reaction, reaction_index):
        task = task if isinstance(task, ReactionTask) else ReactionTask.from_dict(task)
        schema_valid = is_valid_reaction_obj(reaction)
        parsed = None
        unmapped = set()
        if schema_valid:
            parsed, unmapped = self.parser.parse_structured_reaction(reaction)
        output_id = reaction_output_id(reaction)
        target_present, target_correct = self._target_flags(parsed, task)
        rid = parsed["rid"] if parsed else None
        return ReactionCandidate(
            candidate_id=candidate_id(task.task_id, model, reaction_index, output_id),
            output_id=output_id,
            task_id=task.task_id,
            model=model,
            reaction_index=reaction_index,
            reaction=reaction,
            schema_valid=schema_valid,
            parsed=parsed is not None,
            rid=rid,
            cid=task.cid,
            unmapped_names=tuple(sorted(unmapped)),
            target_compound_present=target_present,
            target_position_correct=target_correct,
            balance_success=self._balance_success(parsed),
            duplicate_key=rid or output_id,
        )

    def evaluate_generation(self, task, generation):
        generation_data = generation.to_dict() if hasattr(generation, "to_dict") else generation
        return [
            self.evaluate_reaction(task, generation_data["model"], reaction, idx).to_dict()
            for idx, reaction in enumerate(generation_data.get("reactions") or [])
        ]


class RawReactionValidationService:
    def __init__(self, llm_client, parser=None):
        self.llm_client = llm_client
        self.parser = parser

    def _get_verdicts_safe(self, batch, model, max_retries=3):
        for _ in range(max_retries):
            prompt = f"{VALIDATE_BATCH_INSTRUCT}\n{format_reactions_for_validation(batch)}"
            response = self.llm_client.fetch_answer_str(prompt, model)
            verdicts = parse_indexed_verdicts(response, len(batch))
            if verdicts is not None:
                return verdicts
        return None

    def _rid_for(self, reaction):
        if self.parser is None:
            return None
        if not is_valid_reaction_obj(reaction):
            return None
        parsed, _ = self.parser.parse_structured_reaction(reaction)
        return parsed["rid"] if parsed else None

    def _build_result(self, entry, positives, rounds_done, model):
        confidence = positives / rounds_done if rounds_done > 0 else 0.0
        reaction = entry.get("reaction")
        sample_id = entry.get("sample_id") or entry.get("candidate_id")
        if sample_id is None:
            sample_id = reaction_output_id(reaction)
        result = ValidationResult(
            sample_id=sample_id,
            reaction=reaction,
            validation_model=model,
            validation_rounds=rounds_done,
            validation_positives=positives,
            validation_confidence=confidence,
            validation_votes=tuple(entry.get("_votes", [])),
            rid=entry.get("rid") or self._rid_for(reaction),
            cid=entry.get("cid"),
            task_id=entry.get("task_id"),
            model=entry.get("model"),
        ).to_dict()
        result.update({
            "confidence": confidence,
            "rounds": rounds_done,
            "positives": positives,
            "source": model,
        })
        return result

    def validate_entries(
        self,
        entries,
        *,
        model,
        max_rounds=9,
        acceptance_threshold=0.4,
        early_stop=True,
        batch_size=None,
    ):
        batch_size = batch_size or len(entries) or 1
        active = list(range(len(entries)))
        positives = [0] * len(entries)
        rounds_done = [0] * len(entries)
        votes = [[] for _ in entries]
        results = [None] * len(entries)

        for round_i in range(max_rounds):
            current_active = active if early_stop else list(range(len(entries)))
            if not current_active:
                break

            verdict_map = {}
            for active_batch in batched(current_active, batch_size):
                batch = [entries[i] for i in active_batch]
                verdicts = self._get_verdicts_safe(batch, model)
                if verdicts is None:
                    return None
                for idx, verdict in zip(active_batch, verdicts):
                    verdict_map[idx] = verdict

            remaining = max_rounds - round_i - 1
            next_active = []
            for idx in current_active:
                verdict = verdict_map[idx]
                votes[idx].append(verdict)
                rounds_done[idx] += 1
                positives[idx] += int(verdict)

                if not early_stop:
                    continue

                best_possible = (positives[idx] + remaining) / max_rounds
                worst_possible = positives[idx] / max_rounds
                if (
                    best_possible < acceptance_threshold
                    or worst_possible >= acceptance_threshold
                ):
                    entry_with_votes = {**entries[idx], "_votes": votes[idx]}
                    results[idx] = self._build_result(
                        entry_with_votes, positives[idx], rounds_done[idx], model,
                    )
                else:
                    next_active.append(idx)
            if early_stop:
                active = next_active

        unresolved = active if early_stop else list(range(len(entries)))
        for idx in unresolved:
            if results[idx] is not None:
                continue
            entry_with_votes = {**entries[idx], "_votes": votes[idx]}
            results[idx] = self._build_result(
                entry_with_votes, positives[idx], rounds_done[idx], model,
            )
        return results


class RawReactionRepairService:
    def __init__(self, llm_client):
        self.llm_client = llm_client

    def generate_repairs(self, candidates, *, model, batch_size=5, max_retries=2):
        results = []
        for batch in batched(candidates, batch_size):
            repairs = None
            for _ in range(max_retries):
                prompt = f"{REPAIR_BATCH_INSTRUCT}\n\n{format_repair_candidates(batch)}"
                response = self.llm_client.fetch_answer_str(prompt, model)
                repairs = parse_repair_response(response, len(batch))
                if repairs is not None:
                    break
            if repairs is None:
                return None
            for candidate, repair in zip(batch, repairs):
                if repair is None:
                    status = "null_repair"
                elif not is_valid_reaction_obj(repair):
                    status = "invalid_schema"
                else:
                    status = "ok"
                results.append({
                    "sample_id": candidate.get("sample_id") or candidate.get("rid"),
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
                })
        return results


def aggregate_gold_votes(votes):
    if not votes:
        return {}
    bool_fields = [
        "valid_chemistry",
        "documented_or_standard",
        "target_compound_present",
        "target_position_correct",
        "phase_ok",
        "roles_ok",
        "solvent_catalyst_ok",
    ]
    result = {}
    for field in bool_fields:
        vals = [vote.get(field) for vote in votes if isinstance(vote.get(field), bool)]
        result[field] = (sum(vals) / len(vals) >= 0.5) if vals else None
    positives = [vote.get("valid_chemistry") is True for vote in votes]
    result["gold_rounds"] = len(votes)
    result["gold_positives"] = sum(1 for item in positives if item)
    result["gold_confidence"] = result["gold_positives"] / len(votes)
    categories = Counter(
        vote.get("error_category") or "none"
        for vote in votes
        if isinstance(vote, dict)
    )
    result["error_category"] = categories.most_common(1)[0][0] if categories else None
    return result
