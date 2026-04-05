import re
from pathlib import Path


SQL_PATH = Path(__file__).resolve().parents[1] / "sql" / "create_db.sql"


def _load_sql():
    return SQL_PATH.read_text()


def test_schema_declares_expected_tables():
    sql = _load_sql()
    tables = re.findall(r"CREATE TABLE\s+([a-z_]+)\s*\(", sql, flags=re.IGNORECASE)

    assert set(tables) == {
        "compounds",
        "compound_cas",
        "compound_fingerprints",
        "compound_synonyms",
        "compound_wiki",
        "compound_nfpa",
        "categories",
        "compound_categories",
        "compound_hazard_pictograms",
        "compound_properties",
        "compound_descriptions",
        "compound_commonness_sorting",
        "compound_complexity_sorting",
        "compound_curiosity_sorting",
        "reactions",
        "reaction_reactants",
        "reaction_products",
        "reaction_solvents",
        "reaction_catalysts",
        "reaction_details",
    }


def test_schema_contains_core_compound_and_reaction_constraints():
    sql = _load_sql()

    assert "cid INTEGER PRIMARY KEY" in sql
    assert "inchikey CHAR(27) UNIQUE NOT NULL" in sql
    assert "bits BIGINT[16] NOT NULL" in sql
    assert "rid CHAR(24) PRIMARY KEY" in sql
    assert "balanced BOOLEAN NOT NULL" in sql
    assert "coeff INTEGER CHECK (coeff > 0 OR coeff IS NULL)" in sql
    assert "reaction_details (\n    rid CHAR(24) PRIMARY KEY REFERENCES reactions(rid)," in sql


def test_schema_uses_foreign_keys_for_cross_table_links():
    sql = _load_sql()

    expected_refs = [
        "REFERENCES compounds(cid)",
        "REFERENCES categories(code)",
        "REFERENCES reactions(rid)",
    ]

    for ref in expected_refs:
        assert ref in sql
