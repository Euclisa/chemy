import json

from experiments.model_quality.helpers import (
    build_repair_candidates,
    fixed_round_validate_entries,
    generate_repairs,
    read_jsonl,
    validation_inputs_from_repairs,
    write_jsonl,
)
from experiments.model_quality.run_audit import run_audit_group
from experiments.model_quality.run_repairs import run_validation_route
from experiments.model_quality.sample import build_samples


def reaction(name="Water", product="Hydrogen"):
    return {
        "reagents": [{"name": name, "phase": "l"}],
        "products": [{"name": product, "phase": "g"}],
        "solvent": None,
        "catalyst": None,
    }


class FakeLLMClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def fetch_answer_str(self, message, model):
        self.calls.append({"message": message, "model": model})
        if not self.responses:
            raise AssertionError("Unexpected LLM call")
        return self.responses.pop(0)


def test_sampling_is_deterministic_and_dedupes_by_rid(tmp_path):
    verdict_path = tmp_path / "verdict.jsonl"
    rows = [
        {
            "rid": "r1",
            "cid": 1,
            "reaction": reaction("A", "B"),
            "confidence": 0.25,
            "rounds": 4,
            "positives": 1,
        },
        {
            "rid": "r1",
            "cid": 1,
            "reaction": reaction("A", "C"),
            "confidence": 0.375,
            "rounds": 8,
            "positives": 3,
        },
        {
            "rid": "r2",
            "cid": 2,
            "reaction": reaction("D", "E"),
            "confidence": 0.3333333333,
            "rounds": 6,
            "positives": 2,
        },
        {
            "rid": "r3",
            "cid": 3,
            "reaction": reaction("F", "G"),
            "confidence": 0.9,
            "rounds": 5,
            "positives": 5,
        },
    ]
    write_jsonl(rows, verdict_path)

    config = {
        "seed": 7,
        "verdict_path": str(verdict_path),
        "samples": {
            "borderline_reject": {
                "confidence_min": 0.2,
                "confidence_max": 0.4,
                "n": 2,
            }
        },
        "repair": {"confidence_min": 0.2, "confidence_max": 0.4, "n": 2},
    }

    data_dir_a = tmp_path / "run_a"
    data_dir_b = tmp_path / "run_b"
    build_samples(config, data_dir=data_dir_a)
    build_samples(config, data_dir=data_dir_b)

    sample_a = read_jsonl(data_dir_a / "samples" / "audit_borderline_reject.jsonl")
    sample_b = read_jsonl(data_dir_b / "samples" / "audit_borderline_reject.jsonl")

    assert [entry["sample_id"] for entry in sample_a] == [
        entry["sample_id"] for entry in sample_b
    ]
    assert {entry["sample_id"] for entry in sample_a} == {"r1", "r2"}
    r1 = next(entry for entry in sample_a if entry["sample_id"] == "r1")
    assert r1["confidence"] == 0.375


def test_fixed_round_validation_does_not_early_stop():
    entries = [
        {
            "sample_id": "r1",
            "cid": 1,
            "reaction": reaction(),
            "confidence": 1.0,
        }
    ]
    client = FakeLLMClient(["1. Valid", "1. Valid", "1. Valid"])

    results = fixed_round_validate_entries(
        client,
        entries,
        model="strong",
        rounds=3,
        batch_size=1,
    )

    assert len(client.calls) == 3
    assert results[0]["validation_rounds"] == 3
    assert results[0]["validation_positives"] == 3
    assert results[0]["validation_confidence"] == 1.0


def test_frozen_repairs_feed_the_same_object_to_multiple_validators(tmp_path):
    fixed = reaction("Fixed reagent", "Fixed product")
    candidate = {
        "sample_id": "r1",
        "rid": "r1",
        "cid": 1,
        "reaction": reaction("Bad reagent", "Bad product"),
        "confidence": 0.25,
    }
    repair_client = FakeLLMClient([json.dumps(fixed)])

    repairs = generate_repairs(
        repair_client,
        build_repair_candidates([candidate]),
        model="cheap",
        batch_size=1,
    )
    inputs = validation_inputs_from_repairs(repairs)
    assert inputs[0]["reaction"] == fixed

    cheap_client = FakeLLMClient(["1. Valid"])
    strong_client = FakeLLMClient(["1. Valid"])
    output_a = tmp_path / "cheap.jsonl"
    output_b = tmp_path / "strong.jsonl"

    run_validation_route(
        cheap_client,
        route="cheap_fix_cheap_validator",
        inputs=inputs,
        output_path=output_a,
        model="cheap",
        rounds=1,
        batch_size=1,
        acceptance_threshold=0.4,
    )
    run_validation_route(
        strong_client,
        route="cheap_fix_strong_validator",
        inputs=inputs,
        output_path=output_b,
        model="strong",
        rounds=1,
        batch_size=1,
        acceptance_threshold=0.4,
    )

    assert "Fixed reagent" in cheap_client.calls[0]["message"]
    assert "Fixed reagent" in strong_client.calls[0]["message"]
    assert read_jsonl(output_a)[0]["fixed_reaction"] == fixed
    assert read_jsonl(output_b)[0]["fixed_reaction"] == fixed


def test_audit_resume_skips_completed_rows(tmp_path):
    sample_path = tmp_path / "sample.jsonl"
    output_path = tmp_path / "audit.jsonl"
    write_jsonl(
        [
            {"sample_id": "done", "cid": 1, "reaction": reaction(), "confidence": 0.25},
            {"sample_id": "todo", "cid": 2, "reaction": reaction("A", "B"), "confidence": 0.25},
        ],
        sample_path,
    )
    write_jsonl(
        [
            {
                "sample_id": "done",
                "validation_confidence": 1.0,
            }
        ],
        output_path,
    )
    client = FakeLLMClient(["1. Invalid"])

    summary = run_audit_group(
        client,
        group_name="borderline_reject",
        sample_path=sample_path,
        output_path=output_path,
        model="strong",
        rounds=1,
        batch_size=5,
    )

    assert summary["ran"] == 1
    assert len(client.calls) == 1
    output = read_jsonl(output_path)
    assert [entry["sample_id"] for entry in output] == ["done", "todo"]
    assert output[-1]["validation_confidence"] == 0.0
