from pathlib import Path
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, GraphDescriptors, inchi, rdMolDescriptors

from scripts.data import CompoundStore, ReactionStore
from scripts.infra import JsonlStore, Logger
from scripts.ops.bigsol_parser import BigSolParser
from scripts.ops.conflicts import ConflictResolver
from scripts.ops.hazards import HazardAssembler
from scripts.ops.llm_fetch import LLMFetchOps
from scripts.ops.misc import MiscOps
from scripts.ops.ord_parser import OrdParser
from scripts.ops.properties_compiler import PropertiesCompiler
from scripts.ops.pubchem_fetch import PubChemFetcher
from scripts.ops.reaction_llm_parser import ReactionLLMParser
from scripts.ops.sql_export import SqlExporter
from scripts.ops.solubility import SolubilityGenerator
from scripts.ops.thermo import ThermoOps


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


def make_reaction(store, reagents, products, *, source, confidence=0.5):
    reaction = {
        "reagents": [
            {
                "cid": chem["cid"],
                "original_name": chem["cmpdname"],
                "norm_name": chem["cmpdname"].lower(),
            }
            for chem in reagents
        ],
        "products": [
            {
                "cid": chem["cid"],
                "original_name": chem["cmpdname"],
                "norm_name": chem["cmpdname"].lower(),
            }
            for chem in products
        ],
    }
    reaction = store.assemble_reaction(reaction)
    reaction["source"] = source
    reaction["confidence"] = confidence
    return reaction


@pytest.fixture
def sample_data_dir(tmp_path):
    data_dir = tmp_path / "data"
    dirs = [
        "assets/chems_properties",
        "chems",
        "chems_properties",
        "chems_properties/llm",
        "chems_properties/wiki",
        "misc",
        "reactions_details",
        "reactions_parsed",
        "structures",
        "thermo",
    ]
    for rel in dirs:
        (data_dir / rel).mkdir(parents=True, exist_ok=True)
    return data_dir


@pytest.fixture
def core_context(sample_data_dir):
    store = JsonlStore()
    logger = Logger()
    compounds = CompoundStore(str(sample_data_dir), store, logger)
    reactions = ReactionStore(str(sample_data_dir), compounds, store, logger)

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
    unknown = make_chem(
        -1,
        "Mystery Solvent",
        "CC",
        synonyms=["Mystery Solvent"],
        complexity=5,
    )

    compounds.update_chems([water, ethanol, methanol, unknown])
    store.write_jsonl(
        [
            {"cid": 101, "symbol": "H", "atom_count": 2},
            {"cid": 102, "symbol": "O", "atom_count": 2},
            {"cid": 103, "symbol": "C", "atom_count": 1},
        ],
        compounds.elements_fn,
        backup=False,
    )

    return SimpleNamespace(
        data_dir=sample_data_dir,
        store=store,
        logger=logger,
        compounds=compounds,
        reactions=reactions,
        water=water,
        ethanol=ethanol,
        methanol=methanol,
        unknown=unknown,
    )


@pytest.fixture
def compiler_context(core_context):
    data_dir = Path(core_context.data_dir)
    store = core_context.store

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
        core_context.compounds,
        core_context.reactions,
        store,
        core_context.logger,
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

    return SimpleNamespace(
        **core_context.__dict__,
        compiler=compiler,
        crc=crc,
        solubility=solubility,
        hazards=hazards,
    )


@pytest.fixture
def ops_context(compiler_context):
    data_dir = str(compiler_context.data_dir)

    bigsol = BigSolParser(
        data_dir,
        compiler_context.compounds,
        compiler_context.store,
        compiler_context.logger,
    )
    hazards = HazardAssembler(
        data_dir,
        compiler_context.compounds,
        compiler_context.store,
        compiler_context.logger,
    )
    reaction_llm = ReactionLLMParser(
        data_dir,
        compiler_context.compounds,
        compiler_context.reactions,
        compiler_context.store,
        compiler_context.logger,
    )
    solubility = SolubilityGenerator(
        data_dir,
        compiler_context.compounds,
        compiler_context.store,
        compiler_context.logger,
        compiler_context.crc,
        bigsol,
    )
    llm_client = SimpleNamespace(
        completion_tokens_total=0,
        fetch_answer_str=lambda *args, **kwargs: "",
        _fetch_answer=lambda *args, **kwargs: "",
        restrict_ctt=lambda max_ctt: None,
        submit_entries_to_llm=lambda *args, **kwargs: None,
    )
    llm = LLMFetchOps(
        compiler_context.compounds,
        compiler_context.reactions,
        compiler_context.store,
        compiler_context.logger,
        compiler_context.compiler,
        llm_client,
        reaction_llm,
    )
    misc = MiscOps(
        data_dir,
        compiler_context.compounds,
        compiler_context.reactions,
        compiler_context.store,
        compiler_context.logger,
        compiler_context.compiler,
    )
    conflicts = ConflictResolver(
        compiler_context.compounds,
        compiler_context.logger,
        compiler_context.store,
        compiler_context.compiler,
    )
    thermo = ThermoOps(data_dir, compiler_context.compounds)
    ord_parser = OrdParser(
        data_dir,
        compiler_context.compounds,
        compiler_context.reactions,
        compiler_context.store,
        compiler_context.logger,
    )
    pubchem = PubChemFetcher(
        data_dir,
        compiler_context.compounds,
        compiler_context.store,
        compiler_context.logger,
    )
    sql = SqlExporter(
        compiler_context.compounds,
        compiler_context.reactions,
        compiler_context.store,
        compiler_context.logger,
        compiler_context.compiler,
        misc,
        hazards,
        db_name="chemy_test",
    )

    return SimpleNamespace(
        **compiler_context.__dict__,
        bigsol=bigsol,
        conflicts_ops=conflicts,
        hazards_ops=hazards,
        reaction_llm=reaction_llm,
        solubility_ops=solubility,
        llm_client=llm_client,
        llm_ops=llm,
        misc_ops=misc,
        ord_ops=ord_parser,
        pubchem_ops=pubchem,
        sql_ops=sql,
        thermo_ops=thermo,
    )
