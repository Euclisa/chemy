from tests.conftest import make_chem


def test_resolve_conflicting_synonyms_updates_filtered_synonyms(ops_context, monkeypatch):
    chem_a = make_chem(11, "Alpha", "CCN", synonyms=["Alpha", "Shared Name"], complexity=15)
    chem_b = make_chem(12, "Beta", "CCCl", synonyms=["Beta", "Shared Name"], complexity=15)
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [chem_a, chem_b])

    monkeypatch.setattr(
        ops_context.compiler,
        "get_parsed_reactions_participants_norm_names",
        lambda: {"sharedname"},
    )
    monkeypatch.setattr(ops_context.compounds, "display_compare_table", lambda *args, **kwargs: None)
    decisions = iter(["s1"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(decisions))

    ops_context.conflicts_ops.resolve_conflicting_synonyms()

    chem_b_after = ops_context.compounds.cid_chem_map[12]
    assert chem_b_after["cmpdsynonym"] == ["Beta"]
    assert ops_context.compounds.filtered_synonyms_by_cid[12] == {"sharedname"}
    decisions = ops_context.compounds.synonym_filter_decisions
    assert decisions == (
        {
            "cid": 12,
            "norm_synonym": "sharedname",
            "raw_synonyms": ["Shared Name"],
            "decision": "remove_synonym",
            "reason": "manual_conflict",
            "source": "manual_conflict_resolver",
            "other_cids": [11],
        },
    )


def test_find_synonym_conflicts_is_deterministic(ops_context, monkeypatch):
    chem_a = make_chem(31, "Alpha", "CCN", synonyms=["Alpha", "Shared Name"], complexity=15)
    chem_b = make_chem(32, "Beta", "CCCl", synonyms=["Beta", "Shared Name"], complexity=15)
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [chem_a, chem_b])

    monkeypatch.setattr(
        ops_context.compiler,
        "get_parsed_reactions_participants_norm_names",
        lambda: {"sharedname"},
    )

    assert ops_context.conflicts_ops.find_synonym_conflicts() == {"sharedname": (31, 32)}


def test_resolve_conflicting_inchi_merges_and_blacklists_duplicate(ops_context):
    chem_a = make_chem(21, "Dup One", "CCN", synonyms=["Dup One"], complexity=10)
    chem_b = make_chem(22, "Dup Two", "CCN", synonyms=["Dup Two", "Shared Dup"], complexity=10)
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [chem_a, chem_b])

    ops_context.conflicts_ops.resolve_conflicting_inchi()

    assert 22 not in ops_context.compounds.cid_chem_map
    assert 22 in ops_context.compounds.cids_blacklist
    assert "Shared Dup" in ops_context.compounds.cid_chem_map[21]["cmpdsynonym"]
