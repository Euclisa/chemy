import csv
import json
from pathlib import Path


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


def test_reaction_llm_parser_raw_parsing_applies_threshold_and_deduplicates(ops_context):
    valid_reaction = "Water -> Ethanol"
    duplicate_reaction = "Water -> Ethanol"
    low_confidence = "Water -> Methanol"

    with open(ops_context.reaction_llm.raw_reactions_verdict_fn, "w") as f:
        for entry in [
            {"reaction": valid_reaction, "confidence": 0.9, "source": "openai/gpt-oss-120b"},
            {"reaction": duplicate_reaction, "confidence": 0.8, "source": "qwen/qwen3-235b-a22b"},
            {"reaction": low_confidence, "confidence": 0.1, "source": "openai/gpt-oss-120b"},
        ]:
            f.write(json.dumps(entry) + "\n")

    ops_context.reaction_llm.parse_raw_llm_reactions()
    parsed = ops_context.store.load_jsonl(ops_context.reaction_llm.reactions_parsed_llm_fn)

    assert len(parsed) == 1
    assert parsed[0]["source"] == "openai/gpt-oss-120b"
    assert parsed[0]["confidence"] == 0.9


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

    def fake_fetch(message, model, **kwargs):
        calls.append(model)
        if model == ops_context.llm_ops.gpt_oss:
            return "No reactions available"
        return "Water -> Ethanol"

    def fake_refine(messages, model, **kwargs):
        calls.append(f"revalidate:{model}")
        return "Water -> Ethanol"

    monkeypatch.setattr(ops_context.llm_ops, "fetch_llm_answer_str", fake_fetch)
    monkeypatch.setattr(ops_context.llm_client, "_fetch_answer", fake_refine)

    result = ops_context.llm_ops._fetch_raw_reactions(ops_context.water)

    assert result == {"cid": ops_context.water["cid"], "reactions": ["Water -> Ethanol"]}
    assert calls == [
        ops_context.llm_ops.gpt_oss,
        ops_context.llm_ops.qwen,
        f"revalidate:{ops_context.llm_ops.qwen}",
    ]


def test_llm_fetch_submit_entries_delegates_arguments(ops_context):
    captured = {}

    def fake_submit(out_fn, entries, max_workers, routine, logger, routine_args=None, batch_size=None, description="Generation"):
        captured.update(
            {
                "out_fn": out_fn,
                "entries": entries,
                "max_workers": max_workers,
                "routine": routine,
                "logger": logger,
                "routine_args": routine_args,
                "batch_size": batch_size,
                "description": description,
            }
        )

    ops_context.llm_client.submit_entries_to_llm = fake_submit

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
    assert captured["entries"] == [1, 2]
    assert captured["max_workers"] == 3
    assert captured["routine_args"] == ["x"]
    assert captured["batch_size"] == 5
    assert captured["description"] == "Demo"
