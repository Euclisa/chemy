import json

import pytest

from scripts.ops.thermo_solver import solve_formation_thermo
from tests.conftest import make_reaction


def test_compute_formation_value_for_known_compound(ops_context):
    value = ops_context.thermo_ops.compute_formation_value(
        ops_context.water["cid"],
        20.0,
        {
            101: 2.0,
            102: 8.0,
            103: 1.0,
        },
    )

    assert value == 14.0


def test_llm_thermo_batch_parses_numbered_lines_and_retries_invalid(ops_context, monkeypatch):
    reaction = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.methanol],
        source="openai/gpt-oss-120b",
    )

    responses = iter(["not parseable", "1. -10.5, -5e0"])

    def fake_fetch(*args, **kwargs):
        return next(responses)

    monkeypatch.setattr(ops_context.llm_client, "fetch_answer_str", fake_fetch)

    result = ops_context.thermo_llm_ops._get_reactions_thermo_batch([reaction])

    assert result == [
        {
            "rid": reaction["rid"],
            "estimates": {
                "dH": -10.5,
                "dG": -5.0,
                "source": ops_context.reaction_llm.gpt_oss,
                "prompt": "thermo_v2",
            },
        }
    ]


def test_llm_thermo_line_parser_accepts_float_forms(ops_context):
    assert ops_context.thermo_llm_ops._parse_thermo_line("-.25, +1.2e+2") == (-0.25, 120.0)
    assert ops_context.thermo_llm_ops._parse_thermo_line("2) 10., .5") == (10.0, 0.5)


def test_thermo_llm_has_no_legacy_solve_shim(ops_context):
    assert not hasattr(ops_context.thermo_llm_ops, "process_llm_thermo_estimates")
    assert not hasattr(ops_context.thermo_llm_ops, "solve_formation_estimates")
    assert not hasattr(ops_context.thermo_llm_ops, "validate_solve_against_references")
    assert not hasattr(ops_context.thermo_llm_ops, "chems_formation_thermo_llm_fn")
    assert not hasattr(ops_context.thermo_llm_ops, "corrected_reactions_thermo_llm_fn")


def test_llm_thermo_batch_rejects_partial_parse(ops_context, monkeypatch):
    reaction_one = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.methanol],
        source="openai/gpt-oss-120b",
    )
    reaction_two = make_reaction(
        ops_context.reactions,
        [ops_context.methanol],
        [ops_context.ethanol],
        source="openai/gpt-oss-120b",
    )

    monkeypatch.setattr(
        ops_context.llm_client,
        "fetch_answer_str",
        lambda *args, **kwargs: "-10, -5\nbad row",
    )

    assert ops_context.thermo_llm_ops._get_reactions_thermo_batch([reaction_one, reaction_two]) is None


def test_direct_formation_solver_recovers_exact_synthetic_network():
    h2 = 1
    o2 = 2
    h2o = 3
    carbon = 4
    ch4 = 5
    co2 = 6

    reactions = [
        {
            "rid": "water",
            "reagents": [{"cid": h2}, {"cid": o2}],
            "products": [{"cid": h2o}],
        },
        {
            "rid": "methane",
            "reagents": [{"cid": carbon}, {"cid": h2}],
            "products": [{"cid": ch4}],
        },
        {
            "rid": "co2",
            "reagents": [{"cid": carbon}, {"cid": o2}],
            "products": [{"cid": co2}],
        },
    ]
    balances = {
        "water": {h2: 2, o2: 1, h2o: 2},
        "methane": {carbon: 1, h2: 2, ch4: 1},
        "co2": {carbon: 1, o2: 1, co2: 1},
    }
    reaction_estimates = [
        {"rid": "water", "dH": -115.6, "dG": -109.2, "sigma_dH": 1.0, "sigma_dG": 1.0},
        {"rid": "methane", "dH": -17.8, "dG": -12.1, "sigma_dH": 1.0, "sigma_dG": 1.0},
        {"rid": "co2", "dH": -94.1, "dG": -94.3, "sigma_dH": 1.0, "sigma_dG": 1.0},
    ]
    formation_references = [
        {"cid": h2, "reference_type": "element", "dHf": 0.0, "dGf": 0.0},
        {"cid": o2, "reference_type": "element", "dHf": 0.0, "dGf": 0.0},
        {"cid": carbon, "reference_type": "element", "dHf": 0.0, "dGf": 0.0},
    ]

    result = solve_formation_thermo(reactions, balances, reaction_estimates, formation_references)
    formation = {entry["cid"]: entry for entry in result.formation}

    assert formation[h2o]["dHf"] == pytest.approx(-57.8, abs=0.02)
    assert formation[h2o]["dGf"] == pytest.approx(-54.6, abs=0.02)
    assert formation[ch4]["dHf"] == pytest.approx(-17.8, abs=0.02)
    assert formation[ch4]["dGf"] == pytest.approx(-12.1, abs=0.02)
    assert formation[co2]["dHf"] == pytest.approx(-94.1, abs=0.02)
    assert formation[co2]["dGf"] == pytest.approx(-94.3, abs=0.02)


def test_burcat_extracts_cas_matched_formation_reference(ops_context, tmp_path):
    burcat = tmp_path / "BURCAT.THR.txt"
    burcat.write_text(
        """
H2   Hydrogen REFERENCE ELEMENT
H2 REF ELEMENT    T 2/17H  2.   0.   0.   0.G   200.000  6000.000  A   2.01588 1
 0.29328305E+01 0.82659802E-03-0.14640057E-06 0.15409851E-10-0.68879615E-15    2
-0.81305582E+03-0.10243160E+01 0.29840380E+01 0.14841487E-02-0.57900000E-06    3
 0.10364500E-09-0.76727200E-14-0.91793518E+03-0.15835015E+01 0.00000000E+00    4
O2  CALCULATED FROM ORIGINAL VALUES  HF298=0 KJ
O2 REF ELEMENT    TPIS89O  2.   0.   0.   0.G   200.000  6000.000  A  31.99880 1
 0.36122140E+01 0.74834900E-03-0.19820690E-06 0.30184800E-10-0.18745700E-14    2
-0.11977320E+04 0.36150970E+01 0.36255980E+01-0.18782180E-02 0.70554500E-05    3
-0.67635130E-08 0.21555990E-11-0.10475250E+04 0.43052710E+01 0.00000000E+00    4
7732-18-5
H2O REF=Wooley J. RES. NBS 92 (1987), 35. HF298=-241.826+/-0.04 KJ
H2O               L 5/89H   2O   1    0    0G   200.000  6000.000  A  18.01528 1
 0.26770389E+01 0.29731816E-02-0.77376889E-06 0.94433514E-10-0.42689991E-14    2
-0.29885894E+05 0.68825500E+01 0.41986352E+01-0.20364017E-02 0.65203416E-05    3
-0.54879269E-08 0.17719680E-11-0.30293726E+05-0.84900901E+00-0.29084817E+05    4
""".strip(),
        encoding="utf-8",
    )

    refs = ops_context.thermo_burcat_ops.extract_formation_references(str(burcat))
    water_ref = next(entry for entry in refs if entry["cid"] == ops_context.water["cid"])

    assert water_ref["match"] == "cas"
    assert water_ref["burcat_species"] == "H2O"
    assert water_ref["dHf"] == pytest.approx(-57.8, abs=0.1)
    assert water_ref["dGf"] == pytest.approx(-54.6, abs=1.0)


def test_compute_formation_value_returns_none_for_unknown_or_missing_elements(ops_context):
    assert ops_context.thermo_ops.compute_formation_value(9999, 10.0, {101: 1.0}) is None
    assert (
        ops_context.thermo_ops.compute_formation_value(
            ops_context.water["cid"],
            10.0,
            {101: 1.0},
        )
        is None
    )


def test_thermo_reference_policies_select_expected_train_test_refs(ops_context):
    ops = ops_context.thermo_experiments_ops
    refs = [
        {"cid": 1, "source": "burcat", "reference_type": "burcat", "dHf": -1.0, "dGf": -2.0},
        {"cid": 2, "source": "burcat", "reference_type": "burcat", "dHf": -3.0, "dGf": -4.0},
        {"cid": 3, "source": "burcat", "reference_type": "burcat", "dHf": -5.0, "dGf": -6.0},
        {"cid": 4, "source": "burcat", "reference_type": "burcat", "dHf": -7.0, "dGf": -8.0},
    ]
    ops_context.store.write_jsonl(refs, ops.formation_references_fn, backup=False, no_sort=True)

    train, test = ops.select_formation_references("elements-only")
    assert all(entry["reference_type"] == "element" for entry in train)
    assert [entry["cid"] for entry in test] == [1, 2, 3, 4]

    train, test = ops.select_formation_references("burcat-all")
    assert test == []
    assert {entry["cid"] for entry in refs}.issubset({entry["cid"] for entry in train})

    first_train, first_test = ops.select_formation_references("burcat-split", train_fraction=0.5, seed=42)
    second_train, second_test = ops.select_formation_references("burcat-split", train_fraction=0.5, seed=42)
    assert [entry["cid"] for entry in first_train] == [entry["cid"] for entry in second_train]
    assert [entry["cid"] for entry in first_test] == [entry["cid"] for entry in second_test]
    assert all(entry["reference_type"] == "element" for entry in first_train[:3])
    assert all(entry["reference_type"] != "element" for entry in first_test)
    assert len([entry for entry in first_train if entry["reference_type"] != "element"]) == 2
    assert len(first_test) == 2
    assert (
        {entry["cid"] for entry in first_train if entry["reference_type"] != "element"}
        & {entry["cid"] for entry in first_test}
    ) == set()


def test_thermo_validation_metrics_cover_errors_and_missing_values(ops_context):
    rows = [
        {
            "subset": "train",
            "dHf_ref": 10.0,
            "dHf_error": 2.0,
            "dGf_ref": 20.0,
            "dGf_error": None,
        },
        {
            "subset": "train",
            "dHf_ref": 12.0,
            "dHf_error": -4.0,
            "dGf_ref": None,
            "dGf_error": None,
        },
        {
            "subset": "test",
            "dHf_ref": 30.0,
            "dHf_error": None,
            "dGf_ref": 40.0,
            "dGf_error": 5.0,
        },
    ]

    metrics = ops_context.thermo_experiments_ops.validation_metrics(rows)

    assert metrics["train"]["dHf"]["count"] == 2
    assert metrics["train"]["dHf"]["missing_count"] == 0
    assert metrics["train"]["dHf"]["bias"] == pytest.approx(-1.0)
    assert metrics["train"]["dHf"]["mae"] == pytest.approx(3.0)
    assert metrics["test"]["dHf"]["count"] == 0
    assert metrics["test"]["dHf"]["missing_count"] == 1
    assert metrics["overall"]["dGf"]["count"] == 1
    assert metrics["overall"]["dGf"]["missing_count"] == 1


def test_thermo_experiment_run_writes_outputs(ops_context, monkeypatch):
    ops = ops_context.thermo_experiments_ops
    reaction_one = make_reaction(
        ops_context.reactions,
        [ops_context.water],
        [ops_context.methanol],
        source="openai/gpt-oss-120b",
    )
    reaction_two = make_reaction(
        ops_context.reactions,
        [ops_context.methanol],
        [ops_context.ethanol],
        source="openai/gpt-oss-120b",
    )
    ops_context.store.write_jsonl(
        [
            {
                "rid": reaction_one["rid"],
                "balance": [
                    {"cid": ops_context.water["cid"], "coeff": 1},
                    {"cid": ops_context.methanol["cid"], "coeff": 1},
                ],
            },
            {
                "rid": reaction_two["rid"],
                "balance": [
                    {"cid": ops_context.methanol["cid"], "coeff": 1},
                    {"cid": ops_context.ethanol["cid"], "coeff": 1},
                ],
            },
        ],
        ops_context.reactions.reactions_balance_fn,
        backup=False,
    )
    ops_context.reactions.clear_cached_property("reactions_balance")
    monkeypatch.setattr(ops, "select_solve_reactions", lambda _estimates: [reaction_one, reaction_two])

    references = ops._element_formation_references() + [
        {
            "cid": ops_context.water["cid"],
            "source": "burcat",
            "reference_type": "burcat",
            "dHf": -57.8,
            "dGf": -54.6,
            "sigma_dHf": 1.0,
            "sigma_dGf": 1.0,
        },
        {
            "cid": ops_context.methanol["cid"],
            "source": "burcat",
            "reference_type": "burcat",
            "dHf": -48.1,
            "dGf": -38.9,
            "sigma_dHf": 1.0,
            "sigma_dGf": 1.0,
        },
    ]
    ops_context.store.write_jsonl(references, ops.formation_references_fn, backup=False, no_sort=True)
    reaction_estimates = [
        {"rid": reaction_one["rid"], "dH": 9.7, "dG": 15.7, "sigma_dH": 1.0, "sigma_dG": 1.0},
        {"rid": reaction_two["rid"], "dH": 10.0, "dG": 11.0, "sigma_dH": 1.0, "sigma_dG": 1.0},
    ]

    result = ops.run_policy(
        "burcat-split",
        train_fraction=0.5,
        seed=7,
        name="tiny_split",
        reaction_estimates=reaction_estimates,
    )

    experiment_dir = ops_context.data_dir / "thermo" / "experiments" / "tiny_split"
    assert result["name"] == "tiny_split"
    assert (experiment_dir / "config.json").exists()
    assert (experiment_dir / "chems_formation.jsonl").exists()
    assert (experiment_dir / "validation.jsonl").exists()
    assert (experiment_dir / "metrics.json").exists()

    with open(experiment_dir / "config.json") as f:
        config = json.load(f)
    assert config["policy"] == "burcat-split"
    assert config["formation_references_test"] == 1

    validation = ops_context.store.load_jsonl(str(experiment_dir / "validation.jsonl"))
    assert {row["subset"] for row in validation} == {"train", "test"}
