import pytest

from tests.conftest import make_reaction


def test_reaction_hash_is_order_invariant_within_parts(core_context):
    reactions = core_context.reactions

    reaction_a = {
        "reagents": [{"cid": 1}, {"cid": 2}],
        "products": [{"cid": 3}],
    }
    reaction_b = {
        "reagents": [{"cid": 2}, {"cid": 1}],
        "products": [{"cid": 3}],
    }

    assert reactions.get_reaction_hash(reaction_a) == reactions.get_reaction_hash(reaction_b)


def test_validate_reaction_rejects_duplicates_and_bad_ids(core_context):
    reactions = core_context.reactions

    valid = make_reaction(
        reactions,
        [core_context.water, core_context.ethanol],
        [core_context.methanol],
        source="ord",
    )
    duplicate = {
        **valid,
        "reagents": [valid["reagents"][0], valid["reagents"][0]],
    }
    bad_rid = {**valid, "rid": "not-the-real-hash"}

    assert reactions.validate_reaction(valid)
    assert not reactions.validate_reaction(duplicate)
    assert not reactions.validate_reaction(bad_rid)


def test_assemble_reaction_sets_hash_and_complexity(core_context):
    reactions = core_context.reactions

    reaction = reactions.assemble_reaction(
        {
            "reagents": [
                {"cid": core_context.water["cid"], "original_name": "Water", "norm_name": "water"},
            ],
            "products": [
                {
                    "cid": core_context.ethanol["cid"],
                    "original_name": "Ethanol",
                    "norm_name": "ethanol",
                },
            ],
        }
    )

    expected_complexity = (
        core_context.water["bertz_complexity"] + core_context.ethanol["bertz_complexity"]
    ) / 2

    assert reaction["rid"] == reactions.get_reaction_hash(reaction)
    assert reaction["complexity"] == expected_complexity


def test_convert_details_to_canonic_sets_defaults(core_context):
    reactions = core_context.reactions

    canonical = reactions.convert_details_to_canonic({"rid": "abc", "source": "ord"})

    assert canonical == {
        "rid": "abc",
        "source": "ord",
        "confidence": None,
        "description": None,
        "catalysts": [],
        "solvents": [],
        "provenance": {"doi": None, "patent": None},
    }

    with pytest.raises(Exception, match="Field 'source' must be present"):
        reactions.convert_details_to_canonic({"rid": "abc"})
