import csv
import json
from pathlib import Path

from scripts.ops.raw_reactions.fetcher import parse_reactions_jsonl
from tests.conftest import make_reaction
from scripts.ops.raw_reactions.presets import (
    DEFAULT_PRESET_NAMES,
    POSITION_ANY,
    POSITION_PRODUCT,
    POSITION_REAGENT,
    get_preset,
)
from scripts.ops.raw_reactions.prompts import build_fetch_prompt, format_structured_reaction
from scripts.ops.raw_reactions.validator import ACCEPT_THR


def test_bigsol_parser_groups_valid_rows_and_converts_units(ops_context):
    with open(ops_context.bigsol.big_sol_fn, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["PubChem_CID", "SMILES_Solvent", "Solubility(mol/L)", "Temperature_K"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "PubChem_CID": ops_context.ethanol["cid"],
                "SMILES_Solvent": ops_context.water["smiles"],
                "Solubility(mol/L)": "2.0",
                "Temperature_K": "298.15",
            }
        )
        writer.writerow(
            {
                "PubChem_CID": "9999",
                "SMILES_Solvent": ops_context.water["smiles"],
                "Solubility(mol/L)": "1.0",
                "Temperature_K": "300",
            }
        )
        writer.writerow(
            {
                "PubChem_CID": ops_context.ethanol["cid"],
                "SMILES_Solvent": "unknown",
                "Solubility(mol/L)": "1.0",
                "Temperature_K": "300",
            }
        )

    ops_context.bigsol.parse_bigsol()
    entries = ops_context.store.load_jsonl(ops_context.bigsol.big_sol_parsed_fn)

    assert len(entries) == 1
    assert entries[0]["cid"] == ops_context.ethanol["cid"]
    water_cid = ops_context.water["cid"]
    expected = 2.0 * ops_context.ethanol["mw"]
    assert entries[0]["solubility"][str(water_cid)] == {"298.15": expected}


def test_solubility_generator_merges_crc_and_bigsol_sources(ops_context):
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.ethanol["cid"],
                "solubility": [{"cid": ops_context.water["cid"], "solubility": "miscible"}],
            }
        ],
        ops_context.crc.crc_inorganic_constants_fn,
        backup=False,
    )
    ops_context.store.write_jsonl([], ops_context.crc.crc_organic_constants_fn, backup=False)
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.methanol["cid"],
                "name": "Methanol",
                "solubility": {
                    str(ops_context.water["cid"]): {
                        280.0: 11.2,
                        298.15: 17.3,
                        320.0: 22.4,
                    }
                },
            }
        ],
        ops_context.bigsol.big_sol_parsed_fn,
        backup=False,
    )

    ops_context.solubility_ops.generate_chems_solubility()
    entries = {
        entry["cid"]: entry["solubility"]
        for entry in ops_context.store.load_jsonl(ops_context.solubility_ops.chems_solubility_fn)
    }

    assert entries[ops_context.ethanol["cid"]] == [
        {
            "solvent_cid": ops_context.water["cid"],
            "solvent_name": "Water",
            "value": "miscible",
        }
    ]
    assert entries[ops_context.methanol["cid"]] == [
        {
            "solvent_cid": ops_context.water["cid"],
            "solvent_name": "Water",
            "value": "11.2 g/L (280.0 K); 22.4 g/L (320.0 K)",
        }
    ]


def test_hazard_assembler_prefers_wiki_and_fills_from_llm(ops_context):
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.ethanol["cid"],
                "nfpa": {"health": 1, "flammability": 3, "instability": 0},
                "pictograms": ["GHS02"],
            }
        ],
        ops_context.compounds.chems_hazards_wiki_fn,
        backup=False,
    )
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.ethanol["cid"],
                "categories": ["acute_toxic", "flammable"],
                "source": "openai/gpt-oss-120b",
            },
            {
                "cid": ops_context.methanol["cid"],
                "categories": ["acute_toxic", "flammable"],
                "source": "openai/gpt-oss-120b",
            },
        ],
        ops_context.hazards_ops.chems_hazard_categories_llm_fn,
        backup=False,
    )
    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.methanol["cid"],
                "nfpa": {"health": 2, "flammability": 3, "instability": 0},
                "source": "openai/gpt-oss-120b",
            }
        ],
        ops_context.hazards_ops.chems_nfpa_llm_fn,
        backup=False,
    )

    ops_context.hazards_ops.assemble_chems_hazards()
    hazards = {
        entry["cid"]: entry
        for entry in ops_context.store.load_jsonl(ops_context.hazards_ops.chems_hazards_fn)
    }

    assert hazards[ops_context.ethanol["cid"]]["nfpa"]["source"] == "wikipedia"
    assert hazards[ops_context.ethanol["cid"]]["pictograms"]["value"] == ["GHS02"]
    assert hazards[ops_context.methanol["cid"]]["nfpa"]["value"]["health"] == 2
    assert hazards[ops_context.methanol["cid"]]["pictograms"]["value"] == ["GHS06", "GHS02"]


def test_reaction_llm_parser_parses_control_words_and_unmapped_names(ops_context):
    reaction, unmapped = ops_context.reaction_llm.parse_reaction_scheme("Water + heat -> Ethanol")
    assert reaction is not None
    assert [entry["cid"] for entry in reaction["reagents"]] == [ops_context.water["cid"]]
    assert [entry["cid"] for entry in reaction["products"]] == [ops_context.ethanol["cid"]]
    assert unmapped == set()

    reaction, unmapped = ops_context.reaction_llm.parse_reaction_scheme(
        "Water + Mystery Reagent -> Methanol"
    )
    assert reaction is None
    assert ("mysteryreagent", "Mystery Reagent") in unmapped


def test_raw_reactions_jsonl_parser_requires_phase_and_allows_optional_note():
    valid_without_note = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
        "solvent": None,
        "catalyst": None,
    }
    valid_with_note = {
        "reagents": [{"name": "Hydrogen sulfide", "phase": "g", "note": "bubbled"}],
        "products": [{"name": "Copper(II) sulfide", "phase": "s"}],
    }
    missing_phase = {
        "reagents": [{"name": "Water"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    invalid_phase = {
        "reagents": [{"name": "Water", "phase": "liquid"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    invalid_note = {
        "reagents": [{"name": "Water", "phase": "l", "note": 12}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }

    stats = {}
    parsed = parse_reactions_jsonl(
        "\n".join(
            json.dumps(obj)
            for obj in [
                valid_without_note,
                valid_with_note,
                missing_phase,
                invalid_phase,
                invalid_note,
            ]
        ),
        stats,
    )

    assert parsed == [valid_without_note, {**valid_with_note, "solvent": None, "catalyst": None}]
    assert stats == {"accepted": 2, "invalid_schema": 3}


def test_format_structured_reaction_includes_phase_note_and_top_level_metadata():
    reaction = {
        "reagents": [
            {"name": "hydrogen sulfide", "phase": "g", "note": "bubbled"},
            {"name": "copper(II) sulfate", "phase": "aq"},
        ],
        "products": [{"name": "copper(II) sulfide", "phase": "s"}],
        "solvent": "water",
        "catalyst": None,
    }

    formatted = format_structured_reaction(reaction)

    assert formatted == (
        "hydrogen sulfide(g) {bubbled} + copper(II) sulfate(aq) -> "
        "copper(II) sulfide(s) [solvent: water]"
    )


def test_raw_reactions_preset_names_parse_base_and_position():
    any_preset = get_preset("wiki_crc_rp")
    product_preset = get_preset("wiki_crc_p")
    reagent_preset = get_preset("default_r")

    assert any_preset.base_name == "wiki_crc"
    assert any_preset.name == "wiki_crc_rp"
    assert any_preset.position == POSITION_ANY
    assert product_preset.position == POSITION_PRODUCT
    assert reagent_preset.base_name == "default"
    assert reagent_preset.position == POSITION_REAGENT
    assert "default_rp" in DEFAULT_PRESET_NAMES


def test_raw_reactions_preset_rejects_missing_suffix():
    try:
        get_preset("default")
    except ValueError as e:
        assert "Expected '<base>_<suffix>'" in str(e)
    else:
        raise AssertionError("get_preset accepted an unsuffixed preset name")


def test_fetch_prompt_uses_position_wording_and_context():
    prompt = build_fetch_prompt(
        "Water",
        POSITION_PRODUCT,
        "documented",
        ["Hydrogen + Oxygen -> Water"],
    )

    assert "where Water appears as a product" in prompt
    assert "Generate ADDITIONAL reactions where Water appears as a product" in prompt
    assert "Hydrogen + Oxygen -> Water" in prompt


def test_reaction_llm_parser_structured_parsing_preserves_phase_and_note(ops_context):
    reaction_obj = {
        "reagents": [{"name": "Water", "phase": "l", "note": "neat"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }

    reaction, unmapped = ops_context.reaction_llm.parse_structured_reaction(reaction_obj)

    assert reaction is not None
    assert unmapped == set()
    assert reaction["reagents"][0]["cid"] == ops_context.water["cid"]
    assert reaction["reagents"][0]["phase"] == "l"
    assert reaction["reagents"][0]["note"] == "neat"
    assert reaction["products"][0]["cid"] == ops_context.ethanol["cid"]
    assert reaction["products"][0]["phase"] == "l"
    assert "note" not in reaction["products"][0]


def test_reaction_llm_parser_raw_parsing_filters_invalid_and_deduplicates(ops_context):
    valid_reaction = {
        "reagents": [{"name": "Water", "phase": "l", "note": "neat"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    duplicate_reaction = {
        "reagents": [{"name": "Water", "phase": "aq"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    rejected_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }

    ops_context.store.write_jsonl(
        [
            {
                "reaction": valid_reaction,
                "confidence": ACCEPT_THR,
                "source": "openai/gpt-oss-120b",
            },
            {
                "reaction": duplicate_reaction,
                "confidence": 0.8,
                "source": "qwen/qwen3-235b-a22b",
            },
            {
                "reaction": rejected_reaction,
                "confidence": ACCEPT_THR - 0.01,
                "source": "openai/gpt-oss-120b",
            },
        ],
        ops_context.raw_reactions.layout.verdict(),
        backup=False,
        no_sort=True,
    )

    ops_context.raw_reactions.parse()
    parsed = ops_context.store.load_jsonl(ops_context.reaction_llm.reactions_parsed_llm_fn)

    assert len(parsed) == 1
    assert parsed[0]["source"] == "openai/gpt-oss-120b"
    assert parsed[0]["confidence"] == ACCEPT_THR
    assert "valid" not in parsed[0]
    assert parsed[0]["reagents"][0]["phase"] == "l"
    assert parsed[0]["reagents"][0]["note"] == "neat"


def test_raw_reaction_validator_results_do_not_persist_valid_flag(ops_context, monkeypatch):
    reaction = {
        "cid": ops_context.water["cid"],
        "reaction": {
            "reagents": [{"name": "Water", "phase": "l"}],
            "products": [{"name": "Ethanol", "phase": "l"}],
        },
    }

    monkeypatch.setattr(
        ops_context.llm_client,
        "fetch_answer_str",
        lambda *args, **kwargs: "Valid",
    )

    results = ops_context.raw_reactions._validator.validate_batch(
        [reaction],
        model=ops_context.reaction_llm.gpt_oss,
        max_rounds=1,
        acceptance_threshold=ACCEPT_THR,
    )

    assert results[0]["confidence"] == 1.0
    assert "valid" not in results[0]


def test_checked_in_raw_reaction_verdicts_do_not_persist_valid_flag():
    for path in Path("data/raw_reactions/verdict").glob("*.jsonl"):
        for line in path.read_text().splitlines():
            assert '"valid"' not in line


def test_global_repair_candidates_group_by_rid_and_skip_accepted_or_repaired(ops_context):
    repairer = ops_context.raw_reactions._repairer
    low_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    accepted_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    repaired_reaction = {
        "reagents": [{"name": "Ethanol", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    repaired_rid, _, _ = repairer._get_reaction_id(repaired_reaction)

    ops_context.store.write_jsonl(
        [
            {"cid": ops_context.water["cid"], "reaction": low_reaction, "confidence": 0.2},
            {
                "cid": ops_context.water["cid"],
                "reaction": accepted_reaction,
                "confidence": ACCEPT_THR,
            },
            {
                "cid": ops_context.ethanol["cid"],
                "reaction": repaired_reaction,
                "confidence": 0.3,
            },
            {"cid": ops_context.water["cid"], "reaction": low_reaction, "confidence": 0.35},
        ],
        ops_context.raw_reactions.layout.verdict(),
        backup=False,
        no_sort=True,
    )
    ops_context.store.write_jsonl(
        [{"rid": repaired_rid, "fixed_reaction": None, "reason": "already_done"}],
        ops_context.raw_reactions.layout.repairs(),
        backup=False,
        no_sort=True,
    )

    candidates = repairer._build_candidate_groups(
        lower_bound=0.1,
        acceptance_threshold=ACCEPT_THR,
    )

    assert len(candidates) == 1
    assert candidates[0]["candidate_count"] == 2
    assert candidates[0]["original_confidence"] == 0.35


def test_repair_batch_outputs_success_null_and_validation_failure(ops_context, monkeypatch):
    repairer = ops_context.raw_reactions._repairer

    def make_group(reaction, cid):
        rid, _, _ = repairer._get_reaction_id(reaction)
        return {
            "rid": rid,
            "reaction": reaction,
            "cid": cid,
            "candidate_count": 1,
            "original_confidence": 0.2,
        }

    success_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    rejected_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    null_reaction = {
        "reagents": [{"name": "Ethanol", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    groups = [
        make_group(success_reaction, ops_context.water["cid"]),
        make_group(rejected_reaction, ops_context.water["cid"]),
        make_group(null_reaction, ops_context.ethanol["cid"]),
    ]

    repair_rows = [success_reaction, rejected_reaction, None]

    def fake_repair_fetch(message, *args, **kwargs):
        body = message.split("\n\n", 1)[1]
        assert "RID:" not in body
        assert body.splitlines() == [json.dumps(group["reaction"]) for group in groups]
        return "\n".join(
            "null" if row is None else json.dumps(row)
            for row in repair_rows
        )

    monkeypatch.setattr(
        ops_context.llm_client,
        "fetch_answer_str",
        fake_repair_fetch,
    )
    monkeypatch.setattr(
        repairer.validator,
        "validate_batch",
        lambda entries, **kwargs: [
            {
                **entries[0],
                "confidence": ACCEPT_THR,
                "source": "validator-model",
                "rounds": 1,
                "positives": 1,
            },
            {
                **entries[1],
                "confidence": ACCEPT_THR - 0.01,
                "source": "validator-model",
                "rounds": 1,
                "positives": 0,
            },
        ],
    )

    results = repairer.repair_batch(
        groups,
        repair_model=ops_context.reaction_llm.gpt_oss,
        validation_model=ops_context.reaction_llm.qwen,
        acceptance_threshold=ACCEPT_THR,
    )
    by_rid = {entry["rid"]: entry for entry in results}

    assert by_rid[groups[0]["rid"]]["fixed_reaction"] == success_reaction
    assert by_rid[groups[0]["rid"]]["validation_confidence"] == ACCEPT_THR
    assert "valid" not in by_rid[groups[0]["rid"]]
    assert by_rid[groups[1]["rid"]]["fixed_reaction"] is None
    assert by_rid[groups[1]["rid"]]["reason"] == "validation_rejected"
    assert by_rid[groups[2]["rid"]]["fixed_reaction"] is None
    assert by_rid[groups[2]["rid"]]["reason"] == "null_repair"


def test_parse_uses_global_repairs_only_for_selected_rejected_rids(ops_context):
    selected_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
    }
    unselected_reaction = {
        "reagents": [{"name": "Ethanol", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    failed_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    }
    selected_rid, selected_parsed, _ = ops_context.raw_reactions._repairer._get_reaction_id(
        selected_reaction
    )
    unselected_rid, _, _ = ops_context.raw_reactions._repairer._get_reaction_id(
        unselected_reaction
    )
    failed_rid, _, _ = ops_context.raw_reactions._repairer._get_reaction_id(
        failed_reaction
    )

    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.water["cid"],
                "reaction": selected_reaction,
                "confidence": ACCEPT_THR - 0.01,
                "source": "first-pass",
            },
            {
                "cid": ops_context.water["cid"],
                "reaction": failed_reaction,
                "confidence": ACCEPT_THR - 0.01,
                "source": "first-pass",
            },
        ],
        ops_context.raw_reactions.layout.verdict(),
        backup=False,
        no_sort=True,
    )
    ops_context.store.write_jsonl(
        [
            {
                "rid": selected_rid,
                "fixed_reaction": selected_reaction,
                "fixed_rid": selected_parsed["rid"],
                "candidate_count": 1,
                "original_confidence": ACCEPT_THR - 0.01,
                "repair_source": "repair-model",
                "validation_source": "validator-model",
                "validation_confidence": ACCEPT_THR,
                "rounds": 1,
                "positives": 1,
            },
            {
                "rid": failed_rid,
                "fixed_reaction": None,
                "reason": "validation_rejected",
                "candidate_count": 1,
                "original_confidence": ACCEPT_THR - 0.01,
            },
            {
                "rid": unselected_rid,
                "fixed_reaction": unselected_reaction,
                "fixed_rid": unselected_rid,
                "candidate_count": 1,
                "original_confidence": ACCEPT_THR - 0.01,
                "repair_source": "repair-model",
                "validation_source": "validator-model",
                "validation_confidence": ACCEPT_THR,
                "rounds": 1,
                "positives": 1,
            },
        ],
        ops_context.raw_reactions.layout.repairs(),
        backup=False,
        no_sort=True,
    )

    ops_context.raw_reactions.parse()
    parsed = ops_context.store.load_jsonl(ops_context.reaction_llm.reactions_parsed_llm_fn)

    assert len(parsed) == 1
    assert parsed[0]["rid"] == selected_rid
    assert parsed[0]["confidence"] == ACCEPT_THR
    assert parsed[0]["source"] == "repair:repair-model|validate:validator-model"


def test_llm_fetch_helpers_and_fallback_schedule(ops_context, monkeypatch, tmp_path):
    out_fn = tmp_path / "processed.jsonl"
    ops_context.store.write_jsonl(
        [{"cid": 1, "name": "Water"}, {"cid": 2, "name": "Ethanol"}],
        str(out_fn),
        backup=False,
        no_sort=True,
    )

    assert ops_context.llm_ops._get_processed_entries(str(out_fn), "cid") == {1, 2}
    assert ops_context.llm_ops._get_processed_entries(str(out_fn), lambda x: x.get("name")) == {
        "Water",
        "Ethanol",
    }
    assert ops_context.llm_ops._get_reactions_from_response("A -> B\nnote\nC + D -> E") == [
        "A -> B",
        "C + D -> E",
    ]
    assert ops_context.llm_ops._str_verdict_to_bool("Valid")
    assert not ops_context.llm_ops._str_verdict_to_bool("Invalid")

    calls = []
    structured_reaction = {
        "reagents": [{"name": "Water", "phase": "l"}],
        "products": [{"name": "Ethanol", "phase": "l"}],
        "solvent": None,
        "catalyst": None,
    }

    def fake_fetch(message, model, **kwargs):
        calls.append(model)
        if model == ops_context.llm_ops.gpt_oss:
            return "No reactions available"
        return json.dumps(structured_reaction)

    monkeypatch.setattr(ops_context.llm_client, "fetch_answer_str", fake_fetch)

    result = ops_context.raw_reactions._fetcher._fetch_one_with_fallback(
        ops_context.water, get_preset("default_rp")
    )

    assert result == {"cid": ops_context.water["cid"], "reactions": [structured_reaction]}
    assert calls == [
        ops_context.reaction_llm.gpt_oss,
        ops_context.reaction_llm.qwen,
        ops_context.reaction_llm.qwen,
    ]


def test_existing_reactions_context_is_filtered_by_requested_position(ops_context):
    water_to_ethanol = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.ethanol],
        source="llm",
    )
    ethanol_to_water = make_reaction(
        ops_context.reactions,
        [ops_context.ethanol],
        [ops_context.water],
        source="llm",
    )
    ops_context.store.write_jsonl(
        [water_to_ethanol, ethanol_to_water],
        ops_context.reaction_llm.reactions_parsed_llm_fn,
        backup=False,
    )

    existing = ops_context.raw_reactions._build_existing_reactions_context()
    ops_context.raw_reactions._fetcher.set_existing_reactions(existing)

    assert ops_context.raw_reactions._fetcher._get_existing_context(
        ops_context.water["cid"], POSITION_REAGENT
    ) == ["Water -> Ethanol"]
    assert ops_context.raw_reactions._fetcher._get_existing_context(
        ops_context.water["cid"], POSITION_PRODUCT
    ) == ["Ethanol -> Water"]
    assert ops_context.raw_reactions._fetcher._get_existing_context(
        ops_context.water["cid"], POSITION_ANY
    ) == ["Water -> Ethanol", "Ethanol -> Water"]


def test_llm_fetch_submit_entries_delegates_arguments(ops_context, monkeypatch):
    captured = {}

    def fake_run_batch(
        llm_client,
        entries,
        routine,
        out_fn,
        logger,
        max_workers=1,
        routine_args=None,
        batch_size=None,
        description="Generation",
        store=None,
    ):
        captured.update(
            {
                "llm_client": llm_client,
                "out_fn": out_fn,
                "entries": entries,
                "max_workers": max_workers,
                "routine": routine,
                "logger": logger,
                "routine_args": routine_args,
                "batch_size": batch_size,
                "description": description,
                "store": store,
            }
        )

    monkeypatch.setattr("scripts.ops.llm_fetch.run_batch", fake_run_batch)

    ops_context.llm_ops._submit_entries_to_llm(
        "out.jsonl",
        [1, 2],
        3,
        lambda value: value,
        routine_args=["x"],
        batch_size=5,
        description="Demo",
    )

    assert captured["out_fn"] == "out.jsonl"
    assert captured["llm_client"] == ops_context.llm_client
    assert captured["entries"] == [1, 2]
    assert captured["max_workers"] == 3
    assert captured["routine_args"] == ["x"]
    assert captured["batch_size"] == 5
    assert captured["description"] == "Demo"
