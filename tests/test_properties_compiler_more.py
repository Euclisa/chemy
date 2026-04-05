from tests.conftest import make_chem, make_reaction


def test_compile_chems_properties_includes_identity_solubility_and_crc(compiler_context):
    compiler_context.store.write_jsonl(
        [
            {
                "cid": compiler_context.ethanol["cid"],
                "solubility": [
                    {
                        "solvent_cid": compiler_context.water["cid"],
                        "solvent_name": "Water",
                        "value": "miscible",
                    }
                ],
            }
        ],
        compiler_context.solubility.chems_solubility_fn,
        backup=False,
    )
    compiler_context.store.write_jsonl(
        [
            {
                "cid": compiler_context.water["cid"],
                "name": "Water",
                "physical_form": "liquid",
                "mp": {"value": 0.0},
                "bp": {"value": 100.0},
                "density": 1.0,
                "refractive_index": 1.33,
            }
        ],
        compiler_context.crc.crc_inorganic_constants_fn,
        backup=False,
    )
    compiler_context.store.write_jsonl([], compiler_context.crc.crc_organic_constants_fn, backup=False)
    compiler_context.store.write_jsonl(
        [
            {
                "cid": compiler_context.ethanol["cid"],
                "flash_point": {"value": 13.0},
                "flash_limits": "3.3-19%",
                "ignition_temp": {"value": 363.0},
            }
        ],
        compiler_context.crc.crc_flammability_fn,
        backup=False,
    )

    compiled = {
        entry["cid"]: {prop["property"] for prop in entry["properties"]}
        for entry in compiler_context.compiler.compile_chems_properties()
    }

    assert "PubChem CID" in compiled[compiler_context.water["cid"]]
    assert "Molecular formula" in compiled[compiler_context.water["cid"]]
    assert "Appearence" in compiled[compiler_context.water["cid"]]
    assert "Boiling point" in compiled[compiler_context.water["cid"]]
    assert "Solubility" in compiled[compiler_context.ethanol["cid"]]
    assert "Flash point" in compiled[compiler_context.ethanol["cid"]]


def test_generate_curiosity_index_normalizes_values(compiler_context):
    reactions = [
        make_reaction(
            compiler_context.reactions,
            [compiler_context.water],
            [compiler_context.ethanol],
            source="ord",
        ),
        make_reaction(
            compiler_context.reactions,
            [compiler_context.water],
            [compiler_context.methanol],
            source="openai/gpt-oss-120b",
        ),
    ]
    compiler_context.compiler.write_parsed_reactions(reactions)
    compiler_context.store.write_jsonl(
        [
            {
                "cid": compiler_context.ethanol["cid"],
                "nfpa": {"value": {"health": 2, "flammability": 3, "instability": 0}},
                "pictograms": {"value": ["GHS02", "GHS06"]},
            }
        ],
        compiler_context.hazards.chems_hazards_fn,
        backup=False,
    )

    compiler_context.compiler.generate_curiosity_index()
    curiosity = compiler_context.store.load_jsonl(compiler_context.compiler.chems_curiosity_fn)

    assert len(curiosity) == len(compiler_context.compounds.chems)
    assert max(entry["curiosity"] for entry in curiosity) == 1.0
    assert all(0.0 <= entry["curiosity"] <= 1.0 for entry in curiosity)


def test_discard_orphaned_unmapped_chems_keeps_connected_and_removes_orphans(compiler_context):
    orphan = make_chem(-2, "Orphan", "CCC", synonyms=["Orphan"], complexity=5)
    compiler_context.compounds.update_chems(list(compiler_context.compounds.chems) + [orphan])

    connected_reaction = make_reaction(
        compiler_context.reactions,
        [compiler_context.unknown],
        [compiler_context.ethanol],
        source="ord",
    )
    compiler_context.compiler.write_parsed_reactions([connected_reaction])
    compiler_context.store.write_jsonl([], compiler_context.compiler.reactions_details_ord_fn, backup=False)
    compiler_context.store.write_jsonl([], compiler_context.compiler.reactions_details_llm_fn, backup=False)

    compiler_context.compiler.discard_orphaned_unmapped_chems()

    remaining_unmapped = {chem["cid"] for chem in compiler_context.compounds.chems_unmapped}
    assert compiler_context.unknown["cid"] in remaining_unmapped
    assert -2 not in remaining_unmapped
