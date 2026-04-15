from .data import CompoundStore, ReactionStore
from .infra import JsonlStore, Logger, LLMClient
from .ops.bigsol_parser import BigSolParser
from .ops.raw_reactions import RawReactionsPipeline
from .ops.conflicts import ConflictResolver
from .ops.crc_parser import CRCParser
from .ops.hazards import HazardAssembler
from .ops.llm_fetch import LLMFetchOps
from .ops.misc import MiscOps
from .ops.ord_parser import OrdParser
from .ops.properties_compiler import PropertiesCompiler
from .ops.pubchem_fetch import PubChemFetcher
from .ops.reaction_llm_parser import ReactionLLMParser
from .ops.sql_export import SqlExporter
from .ops.solubility import SolubilityGenerator
from .ops.thermo import ThermoOps
from .ops.thermo_llm import ThermoLLMOps
from .ops.thermo_xtb import ThermoXtbOps


class Chemy:
    def __init__(self, data_dir, api_key=None, db_name=None):
        self.data_dir = data_dir
        self.api_key = api_key
        self.db_name = db_name

        self.store = JsonlStore()
        self.logger = Logger()
        self.compounds = CompoundStore(data_dir, self.store, self.logger)
        self.reactions = ReactionStore(data_dir, self.compounds, self.store, self.logger)
        self.llm_client = LLMClient(api_key, self.logger)

        self.bigsol = BigSolParser(data_dir, self.compounds, self.store, self.logger)
        self.crc = CRCParser(data_dir, self.compounds, self.store, self.logger)
        self.pubchem = PubChemFetcher(data_dir, self.compounds, self.store, self.logger)
        self.ord = OrdParser(data_dir, self.compounds, self.reactions, self.store, self.logger)
        self.reaction_llm = ReactionLLMParser(
            data_dir,
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
        )
        self.solubility = SolubilityGenerator(
            data_dir,
            self.compounds,
            self.store,
            self.logger,
            self.crc,
            self.bigsol,
        )
        self.hazards = HazardAssembler(data_dir, self.compounds, self.store, self.logger)
        self.properties = PropertiesCompiler(
            data_dir,
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
            self.crc,
            self.solubility,
            self.hazards,
        )
        self.llm = LLMFetchOps(
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
            self.properties,
            self.llm_client,
            self.reaction_llm,
        )
        self.misc = MiscOps(
            data_dir,
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
            self.properties,
        )
        self.conflicts = ConflictResolver(
            self.compounds,
            self.logger,
            self.store,
            self.properties,
        )
        self.sql = SqlExporter(
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
            self.properties,
            self.misc,
            self.hazards,
            db_name=db_name,
        )
        self.raw_reactions = RawReactionsPipeline(
            data_dir,
            self.compounds,
            self.reactions,
            self.store,
            self.logger,
            self.llm_client,
            self.reaction_llm,
        )

        self.thermo = ThermoOps(data_dir, self.compounds)
        self.thermo_xtb = ThermoXtbOps(
            data_dir,
            self.compounds,
            self.thermo,
            self.reactions,
            self.store,
            self.logger,
        )
        self.thermo_llm = ThermoLLMOps(
            data_dir,
            self.reactions,
            self.store,
            self.logger,
            self.properties,
            self.thermo,
            self.llm_client,
            self.reaction_llm,
        )

        self.compiler = self.properties
