from tests.conftest import make_chem


def test_clean_and_normalize_chem_name(core_context):
    compounds = core_context.compounds

    assert compounds.clean_chem_name("  12 Sodium\u2011chloride.  ") == "Sodium-chloride"
    assert compounds.normalize_chem_name(" Glacial aluminum solution ") == "aluminium"
    assert compounds.normalize_chem_name(compounds.unknown_name_ph) == compounds.unknown_name_ph


def test_good_name_criteria_filters_bad_names(core_context):
    compounds = core_context.compounds

    assert compounds.good_name_criteria("Acetone")
    assert compounds.good_name_criteria(compounds.unknown_name_ph)
    assert not compounds.good_name_criteria("solid")
    assert not compounds.good_name_criteria("12345")
    assert not compounds.good_name_criteria("ABC-123")


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
