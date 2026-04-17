import json
import sys
from contextlib import contextmanager
from types import ModuleType, SimpleNamespace

from tests.conftest import make_chem


def _long_description():
    return (
        "This reaction record is intentionally verbose so it passes the ORD parsing threshold and "
        "resembles a real synthetic note with enough detail about setup, transformation, workup, "
        "and outcomes to be useful for downstream indexing and export into the chemistry graph."
    )


def test_ord_parser_converts_curated_raw_records_into_graph_artifacts(ops_context):
    raw_fn = ops_context.data_dir / "ord_cleaned.jsonl"
    description = _long_description()

    valid_entry = {
        "reagents": [
            {"smiles": ops_context.water["smiles"], "name": "Water"},
            {"smiles": ops_context.ethanol["smiles"], "name": "Ethanol"},
        ],
        "products": [
            {"smiles": ops_context.methanol["smiles"], "name": "Methanol"},
        ],
        "solvents": [
            {"smiles": ops_context.water["smiles"], "name": "Water"},
            {"smiles": ops_context.water["smiles"], "name": "Water"},
        ],
        "catalysts": [
            {"smiles": ops_context.ethanol["smiles"], "name": "Ethanol"},
        ],
        "description": description,
        "provenance": {"doi": "10.1000/example", "patent": "US-123"},
    }
    unmapped_entry = {
        "reagents": [{"smiles": "N", "name": "Ammonia"}],
        "products": [{"smiles": ops_context.water["smiles"], "name": "Water"}],
        "solvents": None,
        "catalysts": None,
        "description": description,
        "provenance": {"doi": "10.1000/missing", "patent": "US-999"},
    }

    with open(raw_fn, "w") as f:
        for entry in [valid_entry, unmapped_entry]:
            f.write(json.dumps(entry) + "\n")

    ops_context.ord_ops.parse_raw_ord_reactions(str(raw_fn))

    parsed_reactions = ops_context.store.load_jsonl(ops_context.ord_ops.reactions_parsed_ord_fn)
    parsed_details = ops_context.store.load_jsonl(ops_context.ord_ops.reactions_details_ord_fn)
    unmapped = ops_context.store.load_jsonl(ops_context.ord_ops.unmapped_smiles_fn)

    assert len(parsed_reactions) == 1
    assert parsed_reactions[0]["source"] == "ord"
    assert {entry["cid"] for entry in parsed_reactions[0]["reagents"]} == {
        ops_context.water["cid"],
        ops_context.ethanol["cid"],
    }
    assert [entry["cid"] for entry in parsed_reactions[0]["products"]] == [ops_context.methanol["cid"]]

    assert len(parsed_details) == 1
    assert parsed_details[0]["provenance"] == {"doi": "10.1000/example", "patent": "US-123"}
    assert [entry["cid"] for entry in parsed_details[0]["solvents"]] == [ops_context.water["cid"]]
    assert [entry["cid"] for entry in parsed_details[0]["catalysts"]] == [ops_context.ethanol["cid"]]

    assert unmapped == [{"smiles": "N", "name": "Ammonia", "count": 1}]


def test_pubchem_dump_ingestion_adds_new_mapped_compounds_without_overwriting_existing(ops_context):
    dump_fn = ops_context.data_dir / "pubchem_dump.jsonl"
    incoming = [
        {
            "cid": 2,
            "cmpdname": "Duplicate Ethanol",
            "cmpdsynonym": ["Duplicate Ethanol"],
            "mf": "C2H6O",
            "mw": 46.07,
            "complexity": 25.0,
            "smiles": "CCO",
            "inchi": ops_context.ethanol["inchi"],
            "inchikey": ops_context.ethanol["inchikey"],
            "charge": 0,
            "annotation": "duplicate",
            "iupacname": "ethanol",
        },
        {
            "cid": 500,
            "cmpdname": "Acetone",
            "cmpdsynonym": ["Acetone", "Propanone"],
            "mf": "C3H6O",
            "mw": 58.08,
            "complexity": 35.0,
            "smiles": "CC(=O)C",
            "inchi": "InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
            "inchikey": "CSCPPACGZOOCGX-UHFFFAOYSA-N",
            "charge": 0,
            "annotation": "solvent",
            "iupacname": "propan-2-one",
        },
    ]

    with open(dump_fn, "w") as f:
        for obj in incoming:
            for key, value in obj.items():
                f.write(f'  "{key}": {json.dumps(value)},\n')

    ops_context.pubchem_ops.extract_pubchem_dump_to_chems(str(dump_fn))

    assert 500 in ops_context.compounds.cids_mapped
    assert ops_context.compounds.cid_chem_map[500]["cmpdname"] == "Acetone"
    assert ops_context.compounds.cid_chem_map[2]["cmpdname"] == ops_context.ethanol["cmpdname"]


def test_pubchem_direct_fetch_processes_synonyms_before_writing(ops_context, tmp_path, monkeypatch):
    acetone = SimpleNamespace(
        cid=500,
        iupac_name="Acetone",
        synonyms=["Acetone", "solid", "7732-18-5"],
        molecular_formula="C3H6O",
        molecular_weight=58.08,
        charge=0,
        smiles="CC(=O)C",
        inchi="InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
        inchikey="CSCPPACGZOOCGX-UHFFFAOYSA-N",
        complexity=35.0,
    )
    fake_pubchempy = ModuleType("pubchempy")
    fake_pubchempy.get_compounds = lambda cid, namespace: [acetone]
    monkeypatch.setitem(sys.modules, "pubchempy", fake_pubchempy)

    out_fn = tmp_path / "pubchem.jsonl"
    ops_context.pubchem_ops.fetch_chems_cids_from_pubchem([500], str(out_fn))

    entries = ops_context.store.load_jsonl(str(out_fn))
    assert len(entries) == 1
    assert entries[0]["cmpdsynonym"] == ["Acetone"]


def test_reaction_name_lookup_fails_closed_on_ambiguous_synonyms(ops_context):
    alpha = make_chem(501, "Alpha", "CCN", synonyms=["Alpha", "Shared"], complexity=20)
    beta = make_chem(502, "Beta", "CCCl", synonyms=["Beta", "Shared"], complexity=20)
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [alpha, beta])

    reaction, unmapped = ops_context.reaction_llm.parse_structured_reaction({
        "reagents": [{"name": "Shared", "phase": "s"}, {"name": "Water", "phase": "l"}],
        "products": [{"name": "Methanol", "phase": "l"}],
    })

    assert reaction is None
    assert ("shared", "Shared") in unmapped


def test_sql_exporter_emits_relational_rows_from_compiled_pipeline_state(ops_context, monkeypatch):
    reaction = {
        "reagents": [
            {
                "cid": ops_context.water["cid"],
                "original_name": "Water",
                "norm_name": "water",
            }
        ],
        "products": [
            {
                "cid": ops_context.ethanol["cid"],
                "original_name": "Ethanol",
                "norm_name": "ethanol",
            }
        ],
    }
    reaction = ops_context.reactions.assemble_reaction(reaction)
    reaction["source"] = "ord"
    reaction["confidence"] = 0.9

    ops_context.compiler.write_parsed_reactions([reaction])
    ops_context.store.write_jsonl(
        [
            {
                "rid": reaction["rid"],
                "source": "ord",
                "description": "Hydration-style example",
                "confidence": 0.9,
                "catalysts": [{"cid": ops_context.methanol["cid"]}],
                "solvents": [{"cid": ops_context.water["cid"]}],
                "provenance": {"doi": "10.1000/example", "patent": "US-123"},
            }
        ],
        ops_context.compiler.reactions_details_ord_fn,
        backup=False,
    )
    ops_context.store.write_jsonl([], ops_context.compiler.reactions_details_llm_fn, backup=False)
    ops_context.store.write_jsonl(
        [
            {
                "rid": reaction["rid"],
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
    ops_context.compiler.clear_cached_property("reactions_details")
    ops_context.compiler.clear_cached_property("parsed_reactions")
    ops_context.compiler.clear_cached_property("parsed_reactions_balanced")

    alpha = make_chem(501, "Alpha", "CCN", synonyms=["Alpha", "Shared"], complexity=20)
    beta = make_chem(502, "Beta", "CCCl", synonyms=["Beta", "Shared"], complexity=20)
    ops_context.compounds.update_chems(list(ops_context.compounds.chems) + [alpha, beta])

    ops_context.store.write_jsonl(
        [
            {
                "cid": ops_context.ethanol["cid"],
                "nfpa": {"value": {"health": 1, "flammability": 3, "instability": 0}},
                "pictograms": {"value": ["GHS02"]},
            }
        ],
        ops_context.hazards_ops.chems_hazards_fn,
        backup=False,
    )
    ops_context.store.write_jsonl(
        [
            {"code": "solv", "name": "Solvent", "domain": "role"},
        ],
        ops_context.compounds.categories_fn,
        backup=False,
    )
    ops_context.store.write_jsonl(
        [
            {"cid": ops_context.water["cid"], "categories": ["solv"]},
        ],
        ops_context.compounds.chems_categories_fn,
        backup=False,
    )
    ops_context.store.write_jsonl([], ops_context.solubility.chems_solubility_fn, backup=False)
    ops_context.store.write_jsonl([], ops_context.crc.crc_inorganic_constants_fn, backup=False)
    ops_context.store.write_jsonl([], ops_context.crc.crc_organic_constants_fn, backup=False)
    ops_context.store.write_jsonl([], ops_context.crc.crc_flammability_fn, backup=False)
    ops_context.store.write_jsonl(
        [{"cid": ops_context.water["cid"], "description": "Water description"}],
        ops_context.reaction_llm.chems_descriptions_fn,
        backup=False,
    )
    ops_context.store.write_jsonl(
        [
            {"cid": ops_context.water["cid"], "curiosity": 0.7},
            {"cid": ops_context.ethanol["cid"], "curiosity": 0.3},
            {"cid": ops_context.methanol["cid"], "curiosity": 0.1},
            {"cid": ops_context.unknown["cid"], "curiosity": 0.0},
        ],
        ops_context.compiler.chems_curiosity_fn,
        backup=False,
    )

    executed = []

    class FakeCursor:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class FakeConnection:
        def __init__(self):
            self.commits = 0

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def cursor(self):
            return FakeCursor()

        def commit(self):
            self.commits += 1

    fake_conn = FakeConnection()

    fake_psycopg2 = ModuleType("psycopg2")
    fake_psycopg2.connect = lambda database: fake_conn
    fake_extras = ModuleType("psycopg2.extras")
    fake_extras.execute_values = lambda cursor, sql, data: executed.append((sql, list(data)))

    monkeypatch.setitem(sys.modules, "psycopg2", fake_psycopg2)
    monkeypatch.setitem(sys.modules, "psycopg2.extras", fake_extras)

    @contextmanager
    def fake_progress(*args, **kwargs):
        class ProgressStub:
            def add_task(self, description, total):
                return 1

            def advance(self, task):
                return None

            def update(self, task, advance=1):
                return None

            def refresh(self):
                return None

        yield ProgressStub()

    monkeypatch.setattr(ops_context.logger, "progress", fake_progress)

    ops_context.sql_ops.populate_db()

    assert len(executed) == 20

    sql_to_rows = {sql: rows for sql, rows in executed}
    compounds_rows = next(rows for sql, rows in executed if "INSERT INTO compounds" in sql)
    synonyms_rows = next(rows for sql, rows in executed if "INSERT INTO compound_synonyms" in sql)
    reactions_rows = next(rows for sql, rows in executed if "INSERT INTO reactions" in sql)
    details_rows = next(rows for sql, rows in executed if "INSERT INTO reaction_details" in sql)
    solvents_rows = next(rows for sql, rows in executed if "INSERT INTO reaction_solvents" in sql)
    catalysts_rows = next(rows for sql, rows in executed if "INSERT INTO reaction_catalysts" in sql)

    assert any(row[0] == ops_context.water["cid"] for row in compounds_rows)
    assert reactions_rows == [
        (
            reaction["rid"],
            reaction["complexity"],
            "ord",
            True,
            0.9,
        )
    ]
    assert details_rows == [
        (
            reaction["rid"],
            "10.1000/example",
            "US-123",
            "Hydration-style example",
            "ord",
            0.9,
        )
    ]
    assert solvents_rows == [(ops_context.water["cid"], reaction["rid"])]
    assert catalysts_rows == [(ops_context.methanol["cid"], reaction["rid"])]
    assert (501, "Shared") not in synonyms_rows
    assert (502, "Shared") not in synonyms_rows
    assert (501, "Alpha") in synonyms_rows
    assert (502, "Beta") in synonyms_rows
    assert fake_conn.commits == 1
