from tests.conftest import make_chem
from scripts.data import chem_names


def test_chem_names_classifies_and_cleans_synonyms():
    assert chem_names.clean_synonym("Acetone (USP)") == "Acetone"
    assert chem_names.normalize_chem_name("\u03b3 Butyrolactone") == "gammabutyrolactone"

    assert chem_names.classify_name("Acetone").accepted
    assert chem_names.classify_name("7732-18-5").reason == "cas"
    assert chem_names.classify_name("UNII-ABC123DEF4").reason == "identifier"
    assert chem_names.classify_name("III").reason == "generic_word"
    assert chem_names.classify_name("solid").reason == "generic_word"
    assert chem_names.classify_name("A").reason == "generic_word"
    assert chem_names.classify_name("Shared", manual_filtered={"shared"}).reason == "manual_decision"


def test_clean_and_normalize_chem_name(core_context):
    compounds = core_context.compounds

    assert compounds.clean_chem_name("  Sodium\u2011chloride  ") == "Sodium-chloride"
    assert compounds.normalize_chem_name("Aluminum chloride") == "aluminiumchloride"
    assert compounds.normalize_chem_name("  Glacial aluminum solution  ") == "glacialaluminiumsolution"
    assert compounds.normalize_chem_name(compounds.unknown_name_ph) == compounds.unknown_name_ph


def test_good_name_criteria_filters_bad_names(core_context):
    compounds = core_context.compounds

    assert compounds.good_name_criteria("Acetone")
    assert compounds.good_name_criteria(compounds.unknown_name_ph)
    assert not compounds.good_name_criteria("solid")
    assert not compounds.good_name_criteria("12345")
    assert not compounds.good_name_criteria("ABC-123")
    assert not compounds.good_name_criteria("III")


def test_name_resolution_maps_are_explicit_and_fail_closed(core_context):
    compounds = core_context.compounds
    alpha = make_chem(99, "Alpha", "CCN", synonyms=["Alpha", "Shared"], complexity=30)
    beta = make_chem(100, "Beta", "CCCl", synonyms=["Beta", "Shared"], complexity=30)
    compounds.update_chems(list(compounds.chems) + [alpha, beta])

    assert compounds.name_cids_map["shared"] == (99, 100)
    assert "shared" not in compounds.unique_name_cid_map
    assert compounds.ambiguous_name_cids_map["shared"] == (99, 100)
    assert compounds.unique_name_cid_map["alpha"] == 99


def test_synonym_filter_decisions_remove_names_and_replace_cmpdname(core_context):
    compounds = core_context.compounds
    compounds.write_synonym_filter_decisions([
        {
            "cid": 101,
            "norm_synonym": "sharedname",
            "raw_synonyms": ["Shared Name"],
            "decision": "remove_synonym",
            "reason": "manual_conflict",
            "source": "test",
            "other_cids": [102],
        }
    ])
    chem = make_chem(
        101,
        "Shared Name",
        "CCN",
        synonyms=["Shared Name", "Specific Name"],
        complexity=30,
    )

    processed = compounds.process_chem_single(chem)

    assert processed["cmpdname"] == "Specific Name"
    assert processed["cmpdsynonym"] == ["Specific Name"]


def test_update_chems_clears_cached_views(core_context):
    compounds = core_context.compounds
    cached_chems = compounds.chems
    assert "chems" in compounds.__dict__
    assert len(cached_chems) == 4

    acetone = make_chem(
        99,
        "Acetone",
        "CC(=O)C",
        synonyms=["Acetone", "Propanone"],
        cas=["67-64-1"],
        wiki="Acetone",
        iupacname="propan-2-one",
        complexity=30,
    )
    compounds.update_chems(list(cached_chems) + [acetone])

    assert "chems" not in compounds.__dict__
    assert compounds.cid_chem_map[99]["cmpdname"] == "Acetone"
    assert 99 in compounds.cids_mapped


def test_update_cids_blacklist_appends_only_new_positive_ids(core_context):
    compounds = core_context.compounds

    compounds.update_cids_blacklist({core_context.water["cid"], -5, core_context.water["cid"]})
    compounds.update_cids_blacklist({core_context.ethanol["cid"]})

    entries = core_context.store.load_jsonl(compounds.cids_blacklist_fn)
    assert {entry["cid"] for entry in entries} == {
        core_context.water["cid"],
        core_context.ethanol["cid"],
    }


def test_extract_cas_strip_inchi_and_build_chem(core_context):
    compounds = core_context.compounds

    cas_numbers = compounds.extract_chem_cas_numbers(
        {"cmpdsynonym": ["cas-64-17-5", "invalid", "7732-18-5"]}
    )
    assert sorted(cas_numbers) == ["64-17-5", "7732-18-5"]

    assert compounds.strip_inchi_layers("InChI=1S/H2O/h1H2/p-1/q-1") == "InChI=1S/H2O/h1H2"

    built = compounds.build_chem(-10, "Acetone", "CC(=O)C")
    assert built is not None
    assert built["cid"] == -10
    assert built["cmpdname"] == "Acetone"
    assert built["inchi_snone"]
    assert built["inchikey_snone"]


def test_fingerprint_similarity_of_identical_fp_is_one(core_context):
    compounds = core_context.compounds
    fp = core_context.water["ECFP4_fp"]

    assert compounds.get_fp_tanimoto(fp, fp) == 1.0
