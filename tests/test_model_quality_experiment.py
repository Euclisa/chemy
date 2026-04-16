import json
import threading

from experiments.model_quality.helpers import read_jsonl, write_jsonl
from experiments.model_quality.run_gold_audit import run_gold_audit
from experiments.model_quality.run_generation import run_generation
from experiments.model_quality.sample import split_complexity_bands
from experiments.model_quality.score import build_summary


def reaction(name="Water", product="Ethanol"):
    return {
        "reagents": [{"name": name, "phase": "l"}],
        "products": [{"name": product, "phase": "l"}],
        "solvent": None,
        "catalyst": None,
    }


class FakeLLMClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []
        self.input_tokens_total = 0
        self.completion_tokens_total = 0
        self.lock = threading.Lock()

    def fetch_answer_str(self, message, model, **kwargs):
        with self.lock:
            self.calls.append({"message": message, "model": model})
            if not self.responses:
                raise AssertionError("Unexpected LLM call")
            self.input_tokens_total += 10
            self.completion_tokens_total += 5
            return self.responses.pop(0)

    def get_future_result(self, future, executor):
        return future.result()


class FakeParser:
    def parse_structured_reaction(self, reaction_obj):
        name_to_cid = {"Water": 1, "Ethanol": 2, "Methanol": 3}
        parsed = {"reagents": [], "products": [], "rid": "rid:water-ethanol"}
        for side in ("reagents", "products"):
            for compound in reaction_obj[side]:
                cid = name_to_cid.get(compound["name"])
                if cid is None:
                    return None, {(compound["name"].lower(), compound["name"])}
                parsed[side].append({
                    "cid": cid,
                    "original_name": compound["name"],
                    "norm_name": compound["name"].lower(),
                    "phase": compound["phase"],
                })
        return parsed, set()


class FakeReactions:
    def balance_reaction(self, parsed):
        return {"rid": parsed["rid"], "balance": []}


class FakeChemy:
    def __init__(self, data_dir):
        self.reaction_llm = FakeParser()
        self.reactions = FakeReactions()


def test_complexity_band_assignment_uses_low_middle_high_compounds():
    chems = [
        {"cid": 1, "cmpdname": "a", "bertz_complexity": 1, "heavy_count": 1},
        {"cid": 2, "cmpdname": "b", "bertz_complexity": 2, "heavy_count": 2},
        {"cid": 3, "cmpdname": "c", "bertz_complexity": 3, "heavy_count": 3},
        {"cid": 4, "cmpdname": "d", "bertz_complexity": 4, "heavy_count": 4},
        {"cid": 5, "cmpdname": "e", "bertz_complexity": 5, "heavy_count": 5},
        {"cid": 6, "cmpdname": "f", "bertz_complexity": 6, "heavy_count": 6},
    ]

    bands = split_complexity_bands(chems)

    assert [chem["cid"] for chem in bands["simple"]] == [1, 2]
    assert [chem["cid"] for chem in bands["medium"]] == [3, 4]
    assert [chem["cid"] for chem in bands["hard"]] == [5, 6]


def test_generation_bakeoff_keeps_models_separately_attributed(tmp_path, monkeypatch):
    monkeypatch.setattr("experiments.model_quality.run_generation.Chemy", FakeChemy)
    data_dir = tmp_path / "bakeoff"
    task = {
        "task_id": "task:1",
        "cid": 1,
        "compound_name": "Water",
        "preset": "default_rp",
        "position": "any",
        "scope": "documented",
        "prompt": "generate water reactions",
        "existing_reactions": [],
        "complexity_band": "simple",
        "mode": "cold_start",
    }
    write_jsonl([task], data_dir / "tasks.jsonl")
    client = FakeLLMClient([
        json.dumps(reaction()),
        json.dumps(reaction(product="Methanol")),
    ])
    config = {
        "chemy_data_dir": "data",
        "candidate_models": ["model-a", "model-b"],
        "generation_self_review": False,
    }

    manifest = run_generation(
        config,
        data_dir=data_dir,
        llm_client=client,
        max_workers=2,
    )
    candidates = read_jsonl(data_dir / "candidates.jsonl")

    assert manifest["ran"] == 2
    assert {call["model"] for call in client.calls} == {"model-a", "model-b"}
    assert {candidate["model"] for candidate in candidates} == {"model-a", "model-b"}
    assert {candidate["task_id"] for candidate in candidates} == {"task:1"}
    assert manifest["max_workers"] == 2


def test_gold_audit_parallel_batches_aggregate_one_row_per_candidate(tmp_path):
    data_dir = tmp_path / "bakeoff"
    candidates = [
        {
            "candidate_id": "c1",
            "task_id": "t1",
            "model": "model-a",
            "rid": "r1",
            "output_id": "o1",
            "compound_name": "Water",
            "position": "any",
            "reaction": reaction(),
        },
        {
            "candidate_id": "c2",
            "task_id": "t2",
            "model": "model-a",
            "rid": "r2",
            "output_id": "o2",
            "compound_name": "Water",
            "position": "any",
            "reaction": reaction(product="Methanol"),
        },
    ]
    write_jsonl(candidates, data_dir / "candidates.jsonl")
    response = "\n".join([
        json.dumps({
            "index": 1,
            "valid_chemistry": True,
            "documented_or_standard": True,
            "target_compound_present": True,
            "target_position_correct": True,
            "phase_ok": True,
            "roles_ok": True,
            "solvent_catalyst_ok": True,
            "error_category": None,
            "confidence": 0.9,
        }),
        json.dumps({
            "index": 2,
            "valid_chemistry": False,
            "documented_or_standard": False,
            "target_compound_present": True,
            "target_position_correct": True,
            "phase_ok": True,
            "roles_ok": False,
            "solvent_catalyst_ok": True,
            "error_category": "wrong_chemistry",
            "confidence": 0.2,
        }),
    ])
    client = FakeLLMClient([response, response])
    config = {
        "strong_model": "strong",
        "gold_rounds": 2,
        "batch_size": 2,
        "max_workers": 2,
    }

    manifest = run_gold_audit(
        config,
        data_dir=data_dir,
        llm_client=client,
        max_workers=2,
    )
    audits = read_jsonl(data_dir / "gold_audits.jsonl")

    assert manifest["audited"] == 2
    assert manifest["max_workers"] == 2
    assert len(client.calls) == 2
    assert len(audits) == 2
    assert {entry["candidate_id"] for entry in audits} == {"c1", "c2"}
    assert all(entry["gold_rounds"] == 2 for entry in audits)


def test_score_summary_handles_failed_duplicate_and_gold_rejected_outputs(tmp_path):
    data_dir = tmp_path / "bakeoff"
    write_jsonl(
        [
            {
                "generation_id": "g1",
                "task_id": "t1",
                "model": "model-a",
                "parse_stats": {"accepted": 2, "invalid_schema": 1},
                "usage": {"input_tokens": 100, "completion_tokens": 50},
            }
        ],
        data_dir / "generations.jsonl",
    )
    write_jsonl(
        [
            {
                "candidate_id": "c1",
                "task_id": "t1",
                "model": "model-a",
                "rid": "r1",
                "output_id": "o1",
                "parsed": True,
                "target_compound_present": True,
                "target_position_correct": True,
                "balance_success": True,
                "preset": "default_rp",
                "complexity_band": "simple",
            },
            {
                "candidate_id": "c2",
                "task_id": "t1",
                "model": "model-a",
                "rid": "r1",
                "output_id": "o2",
                "parsed": False,
                "target_compound_present": False,
                "target_position_correct": False,
                "balance_success": False,
                "preset": "default_rp",
                "complexity_band": "simple",
            },
        ],
        data_dir / "candidates.jsonl",
    )
    write_jsonl(
        [
            {"candidate_id": "c1", "gold_confidence": 1.0, "model": "model-a"},
            {"candidate_id": "c2", "gold_confidence": 0.0, "model": "model-a"},
        ],
        data_dir / "gold_audits.jsonl",
    )

    summary = build_summary(
        {"acceptance_threshold": 0.4},
        data_dir=data_dir,
    )

    model = summary["models"]["model-a"]
    assert model["schema_success_rate"] == 2 / 3
    assert model["parsed_rate"] == 0.5
    assert model["duplicate_rid_rate"] == 0.5
    assert model["gold_precision"] == 0.5
    assert model["gold_accepted_unique_rid_n"] == 1
