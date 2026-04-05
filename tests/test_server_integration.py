import json
import os
import socket
import subprocess
import time
import uuid
from pathlib import Path
from types import SimpleNamespace

import psycopg2
import pytest
import requests
from psycopg2 import sql
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, GraphDescriptors, inchi, rdMolDescriptors

from scripts.data import CompoundStore, ReactionStore
from scripts.infra import JsonlStore, Logger
from scripts.ops.misc import MiscOps
from scripts.ops.properties_compiler import PropertiesCompiler
from scripts.ops.sql_export import SqlExporter


ROOT = Path(__file__).resolve().parents[1]
SCHEMA_SQL = ROOT / "sql" / "create_db.sql"
CHEM_DEBUG_BIN = ROOT / "build" / "chem_debug"


def make_chem(
    cid,
    name,
    smiles,
    *,
    synonyms=None,
    cas=None,
    wiki=None,
    iupacname=None,
    complexity=50,
):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for test fixture: {smiles}")

    smiles = Chem.MolToSmiles(mol, canonical=True)
    inchi_val = Chem.MolToInchi(mol)
    inchi_snone = inchi.MolToInchi(mol, options="/SNon")
    inchikey = Chem.InchiToInchiKey(inchi_val)
    inchikey_snone = inchi.MolToInchiKey(mol, options="/SNon")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)

    bits = []
    for idx in range(0, 1024, 64):
        chunk = fp.ToBitString()[idx:idx + 64]
        bits.append(int(chunk, 2))

    return {
        "cid": cid,
        "cmpdname": name,
        "cmpdsynonym": list(synonyms or [name]),
        "cas": list(cas or []),
        "wiki": wiki,
        "iupacname": iupacname,
        "smiles": smiles,
        "inchi": inchi_val,
        "inchikey": inchikey,
        "inchi_snone": inchi_snone,
        "inchikey_snone": inchikey_snone,
        "mf": rdMolDescriptors.CalcMolFormula(mol),
        "mw": Descriptors.MolWt(mol),
        "complexity": complexity,
        "bertz_complexity": GraphDescriptors.BertzCT(mol),
        "heavy_count": mol.GetNumHeavyAtoms(),
        "organic": True,
        "charge": Chem.GetFormalCharge(mol),
        "ECFP4_fp": {"bits": bits, "popcount": fp.GetNumOnBits()},
    }


def _pg_conn_kwargs(dbname=None):
    kwargs = {"dbname": dbname or os.getenv("PGDATABASE", "postgres")}

    env_map = {
        "PGHOST": "host",
        "PGPORT": "port",
        "PGUSER": "user",
        "PGPASSWORD": "password",
    }
    for env_name, kw_name in env_map.items():
        value = os.getenv(env_name)
        if value:
            kwargs[kw_name] = int(value) if kw_name == "port" else value

    return kwargs


def _connect_or_skip(dbname=None):
    try:
        return psycopg2.connect(**_pg_conn_kwargs(dbname))
    except Exception as exc:
        pytest.skip(f"PostgreSQL not available for server integration tests: {exc}")


def _create_database(db_name):
    with _connect_or_skip() as conn:
        conn.autocommit = True
        with conn.cursor() as cur:
            cur.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier(db_name)))

    with _connect_or_skip(db_name) as conn:
        with conn.cursor() as cur:
            cur.execute(SCHEMA_SQL.read_text())
        conn.commit()


def _drop_database(db_name):
    with _connect_or_skip() as conn:
        conn.autocommit = True
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT pg_terminate_backend(pid)
                FROM pg_stat_activity
                WHERE datname = %s AND pid <> pg_backend_pid()
                """,
                (db_name,),
            )
            cur.execute(sql.SQL("DROP DATABASE IF EXISTS {}").format(sql.Identifier(db_name)))


def _find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _build_server_binary():
    subprocess.run(
        ["cmake", "--build", "build", "--target", "chem_debug", "-j2"],
        cwd=ROOT,
        check=True,
    )


def _seed_fixture_data(data_dir, db_name):
    store = JsonlStore()
    logger = Logger()
    compounds = CompoundStore(str(data_dir), store, logger)
    reactions = ReactionStore(str(data_dir), compounds, store, logger)

    water = make_chem(
        1,
        "Water",
        "O",
        synonyms=["Water", "Oxidane"],
        cas=["7732-18-5"],
        wiki="Water",
        iupacname="oxidane",
        complexity=10,
    )
    ethanol = make_chem(
        2,
        "Ethanol",
        "CCO",
        synonyms=["Ethanol", "Ethyl alcohol"],
        cas=["64-17-5"],
        wiki="Ethanol",
        iupacname="ethanol",
        complexity=25,
    )
    methanol = make_chem(
        3,
        "Methanol",
        "CO",
        synonyms=["Methanol", "Methyl alcohol"],
        cas=["67-56-1"],
        wiki="Methanol",
        iupacname="methanol",
        complexity=20,
    )
    compounds.update_chems([water, ethanol, methanol])

    crc = SimpleNamespace(
        crc_inorganic_constants_fn=str(data_dir / "chems_properties" / "crc_inorganic.jsonl"),
        crc_organic_constants_fn=str(data_dir / "chems_properties" / "crc_organic.jsonl"),
        crc_flammability_fn=str(data_dir / "chems_properties" / "crc_flammability.jsonl"),
    )
    solubility = SimpleNamespace(
        chems_solubility_fn=str(data_dir / "chems_properties" / "chems_solubility.jsonl"),
    )
    hazards = SimpleNamespace(
        chems_hazards_fn=str(data_dir / "chems_properties" / "chems_hazards.jsonl"),
    )

    compiler = PropertiesCompiler(
        str(data_dir),
        compounds,
        reactions,
        store,
        logger,
        crc,
        solubility,
        hazards,
    )
    store.register_sorting(compiler.reactions_parsed_ord_fn, "rid")
    store.register_sorting(compiler.reactions_parsed_llm_fn, "rid")
    store.register_sorting(compiler.reactions_details_ord_fn, "rid")
    store.register_sorting(compiler.reactions_details_llm_fn, "rid")
    store.register_vault(compiler.reactions_parsed_ord_fn, "reactions_")
    store.register_vault(compiler.reactions_parsed_llm_fn, "reactions_")
    store.register_vault(compiler.reactions_details_ord_fn, "reactions_details_")
    store.register_vault(compiler.reactions_details_llm_fn, "reactions_details_")

    store.write_jsonl(
        [
            {"code": "solvent", "name": "Solvents", "domain": "role"},
            {"code": "alcohol", "name": "Alcohols", "domain": "functional_group"},
        ],
        compounds.categories_fn,
        backup=False,
    )
    store.write_jsonl(
        [
            {"cid": water["cid"], "categories": ["solvent"]},
            {"cid": ethanol["cid"], "categories": ["alcohol", "solvent"]},
            {"cid": methanol["cid"], "categories": ["alcohol", "solvent"]},
        ],
        compounds.chems_categories_fn,
        backup=False,
    )
    store.write_jsonl(
        [
            {
                "cid": ethanol["cid"],
                "nfpa": {
                    "value": {
                        "health": "2",
                        "flammability": "3",
                        "instability": "0",
                    }
                },
                "pictograms": {"value": ["flammable"]},
            }
        ],
        hazards.chems_hazards_fn,
        backup=False,
    )
    store.write_jsonl(
        [
            {"cid": water["cid"], "description": "Water fixture description."},
            {"cid": ethanol["cid"], "description": "Ethanol fixture description."},
            {"cid": methanol["cid"], "description": "Methanol fixture description."},
        ],
        str(data_dir / "chems_properties" / "llm" / "chems_descriptions.jsonl"),
        backup=False,
    )

    reaction = reactions.assemble_reaction(
        {
            "reagents": [
                {
                    "cid": ethanol["cid"],
                    "original_name": ethanol["cmpdname"],
                    "norm_name": ethanol["cmpdname"].lower(),
                }
            ],
            "products": [
                {
                    "cid": methanol["cid"],
                    "original_name": methanol["cmpdname"],
                    "norm_name": methanol["cmpdname"].lower(),
                }
            ],
        }
    )
    reaction["source"] = "ord"
    reaction["confidence"] = 0.91
    compiler.write_parsed_reactions([reaction])
    store.write_jsonl(
        [
            {
                "rid": reaction["rid"],
                "source": "ord",
                "confidence": 0.91,
                "description": "Fixture reaction description.",
                "solvents": [{"cid": water["cid"]}],
                "catalysts": [{"cid": water["cid"]}],
                "provenance": {"doi": "10.1000/test-doi", "patent": "US-TEST"},
            }
        ],
        compiler.reactions_details_ord_fn,
        backup=False,
    )
    store.write_jsonl(
        [
            {
                "rid": reaction["rid"],
                "balance": [
                    {"cid": ethanol["cid"], "coeff": 1},
                    {"cid": methanol["cid"], "coeff": 1},
                ],
            }
        ],
        reactions.reactions_balance_fn,
        backup=False,
    )
    reactions.clear_cached_property("reactions_balance")

    misc = MiscOps(str(data_dir), compounds, reactions, store, logger, compiler)
    exporter = SqlExporter(
        compounds,
        reactions,
        store,
        logger,
        compiler,
        misc,
        hazards,
        db_name=db_name,
    )
    exporter.populate_db()

    db_info = {"dbname": db_name}
    conn_params = _pg_conn_kwargs(db_name)
    if "host" in conn_params:
        db_info["host"] = conn_params["host"]
    if "port" in conn_params:
        db_info["port"] = conn_params["port"]
    if "user" in conn_params:
        db_info["user"] = conn_params["user"]
    if "password" in conn_params:
        db_info["password"] = conn_params["password"]

    (data_dir / "db.json").write_text(json.dumps(db_info))
    (data_dir / "ui" / "index.html").write_text("<!doctype html><html><body>fixture</body></html>")

    return SimpleNamespace(
        data_dir=data_dir,
        water=water,
        ethanol=ethanol,
        methanol=methanol,
        reaction=reaction,
    )


@pytest.fixture(scope="module")
def server_env(tmp_path_factory):
    _connect_or_skip().close()
    _build_server_binary()

    db_name = f"chemy_server_test_{uuid.uuid4().hex[:8]}"
    _create_database(db_name)

    data_dir = tmp_path_factory.mktemp("server-fixture-data")
    for rel in [
        "ui",
        "structures",
        "chems",
        "misc",
        "chems_properties",
        "chems_properties/llm",
        "chems_properties/wiki",
        "assets",
        "assets/structures",
        "reactions_parsed",
        "reactions_details",
    ]:
        (data_dir / rel).mkdir(parents=True, exist_ok=True)

    fixture = _seed_fixture_data(data_dir, db_name)

    port = _find_free_port()
    process = subprocess.Popen(
        [str(CHEM_DEBUG_BIN), str(data_dir), "127.0.0.1", str(port)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    base_url = f"http://127.0.0.1:{port}"
    try:
        deadline = time.time() + 20
        while time.time() < deadline:
            if process.poll() is not None:
                output = process.stdout.read() if process.stdout else ""
                raise RuntimeError(f"Server exited early:\n{output}")
            try:
                response = requests.get(base_url + "/", timeout=0.5)
                if response.status_code in {200, 404}:
                    break
            except requests.RequestException:
                time.sleep(0.2)
        else:
            raise RuntimeError("Timed out waiting for the server to start")

        yield SimpleNamespace(
            base_url=base_url,
            process=process,
            data_dir=data_dir,
            water=fixture.water,
            ethanol=fixture.ethanol,
            methanol=fixture.methanol,
            reaction=fixture.reaction,
        )
    finally:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=10)
        _drop_database(db_name)


def test_compound_info_endpoint_accepts_integer_arrays(server_env):
    response = requests.post(
        server_env.base_url + "/api/compound_info",
        json=[server_env.ethanol["cid"], server_env.methanol["cid"]],
        timeout=5,
    )

    assert response.status_code == 200
    payload = response.json()
    assert [entry["cid"] for entry in payload] == [server_env.ethanol["cid"], server_env.methanol["cid"]]
    assert payload[0]["name"] == "Ethanol"
    assert payload[0]["description"] == "Ethanol fixture description."


def test_compound_and_reaction_info_reject_malformed_bodies(server_env):
    compound_response = requests.post(
        server_env.base_url + "/api/compound_info",
        json={"cid": server_env.ethanol["cid"]},
        timeout=5,
    )
    reaction_response = requests.post(
        server_env.base_url + "/api/reaction_info",
        json={"rid": server_env.reaction["rid"]},
        timeout=5,
    )

    assert compound_response.status_code == 400
    assert "Expected JSON array" in compound_response.json()["error"]

    assert reaction_response.status_code == 400
    assert "Expected JSON array" in reaction_response.json()["error"]


def test_reaction_info_endpoint_returns_reaction_metadata(server_env):
    response = requests.post(
        server_env.base_url + "/api/reaction_info",
        json=[server_env.reaction["rid"]],
        timeout=5,
    )

    assert response.status_code == 200
    payload = response.json()
    assert len(payload) == 1
    assert payload[0]["rid"] == server_env.reaction["rid"]
    assert payload[0]["description"] == "Fixture reaction description."
    assert payload[0]["reactants"][0]["cid"] == server_env.ethanol["cid"]
    assert payload[0]["products"][0]["cid"] == server_env.methanol["cid"]


def test_build_graph_validates_fields_and_returns_expected_edge(server_env):
    invalid_response = requests.post(
        server_env.base_url + "/api/build_graph",
        json={
            "sources": [server_env.ethanol["cid"]],
            "targets": [server_env.methanol["cid"]],
            "primary_only": False,
            "max_paths": -1,
            "max_cost": 1,
        },
        timeout=5,
    )
    valid_response = requests.post(
        server_env.base_url + "/api/build_graph",
        json={
            "sources": [server_env.ethanol["cid"]],
            "targets": [server_env.methanol["cid"]],
            "primary_only": False,
            "max_paths": 5,
            "max_cost": 1,
        },
        timeout=5,
    )

    assert invalid_response.status_code == 400
    assert "max_paths" in invalid_response.json()["error"]

    assert valid_response.status_code == 200
    graph = valid_response.json()["graph"]
    ethanol_entry = next(entry for entry in graph if entry["cid"] == server_env.ethanol["cid"])
    assert any(
        target["cid"] == server_env.methanol["cid"] and server_env.reaction["rid"] in target["reactions"]
        for target in ethanol_entry["targets"]
    )


def test_search_validates_parameters_and_returns_exact_and_fuzzy_matches(server_env):
    invalid_sort = requests.get(
        server_env.base_url + "/api/search",
        params={"sorting_order": "bogus"},
        timeout=5,
    )
    invalid_direction = requests.get(
        server_env.base_url + "/api/search",
        params={"sorting_direction": "sideways"},
        timeout=5,
    )
    invalid_page = requests.get(
        server_env.base_url + "/api/search",
        params={"page": 99},
        timeout=5,
    )
    exact_match = requests.get(
        server_env.base_url + "/api/search",
        params={"q": "Ethyl alcohol"},
        timeout=5,
    )
    fuzzy_match = requests.get(
        server_env.base_url + "/api/search",
        params={"q": "ethanl"},
        timeout=5,
    )

    assert invalid_sort.status_code == 400
    assert "Invalid sorting order" in invalid_sort.json()["error"]

    assert invalid_direction.status_code == 400
    assert "sorting_direction" in invalid_direction.json()["error"]

    assert invalid_page.status_code == 400
    assert "Invalid offset" in invalid_page.json()["error"]

    assert exact_match.status_code == 200
    assert exact_match.json()["cids"][0] == server_env.ethanol["cid"]

    assert fuzzy_match.status_code == 200
    assert server_env.ethanol["cid"] in fuzzy_match.json()["cids"]


def test_adjacent_edges_requires_cid_and_returns_real_neighbors(server_env):
    missing_cid = requests.get(server_env.base_url + "/api/adjacent_edges", timeout=5)
    response = requests.get(
        server_env.base_url + "/api/adjacent_edges",
        params={"cid": server_env.ethanol["cid"]},
        timeout=5,
    )

    assert missing_cid.status_code == 400
    assert "cid" in missing_cid.json()["error"]

    assert response.status_code == 200
    payload = response.json()
    ethanol_entry = next(entry for entry in payload if entry["cid"] == server_env.ethanol["cid"])
    assert [target["cid"] for target in ethanol_entry["targets"]] == [server_env.methanol["cid"]]


def test_structure_route_generates_and_reuses_svg(server_env):
    svg_path = server_env.data_dir / "structures" / f"{server_env.ethanol['cid']}.svg"
    if svg_path.exists():
        svg_path.unlink()

    first = requests.get(
        server_env.base_url + f"/assets/structures/{server_env.ethanol['cid']}.svg",
        timeout=5,
    )
    second = requests.get(
        server_env.base_url + f"/assets/structures/{server_env.ethanol['cid']}.svg",
        timeout=5,
    )

    assert first.status_code == 200
    assert "svg" in first.headers["Content-Type"]
    assert svg_path.exists()

    assert second.status_code == 200
    assert first.text == second.text
