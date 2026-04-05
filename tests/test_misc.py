import json

from tests.conftest import make_chem, make_reaction


def test_misc_sortings_and_background_substances(ops_context):
    reaction_one = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.ethanol],
        source="ord",
    )
    reaction_two = make_reaction(
        ops_context.reactions,
        [ops_context.water, ops_context.ethanol],
        [ops_context.methanol],
        source="openai/gpt-oss-120b",
    )
    ops_context.compiler.write_parsed_reactions([reaction_one, reaction_two])
    ops_context.store.write_jsonl(
        [
            {"cid": ops_context.water["cid"], "curiosity": 0.8},
            {"cid": ops_context.ethanol["cid"], "curiosity": 0.1},
            {"cid": ops_context.methanol["cid"], "curiosity": 0.4},
            {"cid": ops_context.unknown["cid"], "curiosity": 0.0},
        ],
        ops_context.compiler.chems_curiosity_fn,
        backup=False,
    )

    ops_context.misc_ops.get_background_substances(2)

    with open(ops_context.compounds.background_cids_fn) as f:
        background = json.load(f)

    assert background == [ops_context.water["cid"], ops_context.ethanol["cid"]]
    commonness_sort = ops_context.misc_ops.get_commonness_chems_sorting()
    expected_commonness = sorted(
        [chem["cid"] for chem in ops_context.compounds.chems],
        key=lambda cid: ops_context.compiler.get_chems_reactions_occurence()[cid],
    )
    assert commonness_sort == expected_commonness
    assert ops_context.misc_ops.get_curiosity_chems_sorting() == [
        ops_context.unknown["cid"],
        ops_context.ethanol["cid"],
        ops_context.methanol["cid"],
        ops_context.water["cid"],
    ]

    complexity_sort = ops_context.misc_ops.get_complexity_chems_sorting()
    expected = sorted(
        [chem["cid"] for chem in ops_context.compounds.chems],
        key=lambda cid: ops_context.compounds.cid_chem_map[cid]["bertz_complexity"],
    )
    assert complexity_sort == expected


def test_misc_radicals_cids_unbalancing_and_svg_generation(ops_context, tmp_path):
    radical = make_chem(
        77,
        "Radical",
        "C",
        synonyms=["Methyl(.) radical"],
        complexity=5,
    )
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [radical])

    balanced = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.ethanol],
        source="ord",
    )
    unbalanced = make_reaction(
        ops_context.reactions,
        [ops_context.methanol],
        [ops_context.unknown],
        source="openai/gpt-oss-120b",
    )
    ops_context.compiler.write_parsed_reactions([balanced, unbalanced])
    ops_context.store.write_jsonl(
        [
            {
                "rid": balanced["rid"],
                "balance": [
                    {"cid": ops_context.water["cid"], "coeff": 1},
                    {"cid": ops_context.ethanol["cid"], "coeff": 1},
                ],
            }
        ],
        ops_context.reactions.reactions_balance_fn,
        backup=False,
    )
    ops_context.reactions.clear_cached_property("reactions_balance")

    radicals_fn = tmp_path / "radicals.jsonl"
    ops_context.misc_ops.extract_radicals_list(str(radicals_fn))
    radicals = ops_context.store.load_jsonl(str(radicals_fn))
    assert radicals == [{"cid": 77, "name": "Radical"}]

    cids_fn = tmp_path / "cids.txt"
    ops_context.misc_ops.get_chems_cids(str(cids_fn))
    cids = {int(line) for line in cids_fn.read_text().splitlines()}
    assert cids == {chem["cid"] for chem in ops_context.compounds.chems}

    ops_context.misc_ops.find_unbalancing_chems()
    unbalancing = ops_context.store.load_jsonl(ops_context.misc_ops.unbalancing_cids_fn)
    assert {tuple(sorted(entry.items())) for entry in unbalancing} == {
        tuple(sorted({"cid": ops_context.methanol["cid"], "name": "Methanol", "count": 1}.items())),
        tuple(
            sorted(
                {
                    "cid": ops_context.unknown["cid"],
                    "name": "Mystery Solvent",
                    "count": 1,
                }.items()
            )
        ),
    }

    ops_context.misc_ops.generate_chems_structures_svg()
    assert (ops_context.data_dir / "structures" / f"{ops_context.water['cid']}.svg").exists()
