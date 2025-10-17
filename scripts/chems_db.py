import threading
import os
import json
import shutil
from rich.console import Console
from rich import traceback


class ChemsDB:

    def __init__(self, data_dir):
        self.print_lock = threading.Lock()

        self.data_dir = data_dir
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        
        self.structures_dir = os.path.join(self.data_dir, 'structures')
        if not os.path.exists(self.structures_dir):
            os.makedirs(self.structures_dir)
        
        self.tmp_dir = os.path.join(data_dir, 'tmp/')

        self.raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions.jsonl")
        self.wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "wiki_raw_reactions.jsonl")
        self.top_rare_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "top_rare_raw_reactions.jsonl")
        self.raw_reactions_verdict_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions_verdict.jsonl")
        self.products_wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', 'wiki_products_raw_reactions.jsonl')

        self.reactions_parsed_fn = os.path.join(self.data_dir, 'reactions_parsed', "reactions_parsed.jsonl")
        self.reactions_parsed_llm_fn = os.path.join(self.data_dir, 'reactions_parsed', "reactions_parsed_llm.jsonl")
        self.reactions_parsed_ord_fn = os.path.join(self.data_dir, 'reactions_parsed', 'reactions_parsed_ord.jsonl')

        self.chems_descriptions_fn = os.path.join(self.data_dir, 'chems', 'chems_descriptions.jsonl')
        self.chems_fn = os.path.join(self.data_dir, 'chems', "chems.jsonl")
        self.chems_categories_fn = os.path.join(self.data_dir, 'chems', "chems_categories.jsonl")
        self.wiki_chems_fn = os.path.join(self.data_dir, 'chems', "wiki_chems.jsonl")
        self.hazards_chems_fn = os.path.join(self.data_dir, 'chems', "hazards_chems.jsonl")
        self.chems_edges_fn = os.path.join(self.data_dir, 'chems', 'chems_edges.jsonl')
        self.elements_fn = os.path.join(self.data_dir, 'chems', 'elements.jsonl')

        self.categories_fn = os.path.join(self.data_dir, 'misc', "categories.jsonl")
        self.background_cids_fn = os.path.join(self.data_dir, 'misc', 'background_cids.json')
        self.commonness_sorted_cids_fn = os.path.join(self.data_dir, 'misc', 'commonness_sorted_cids.json')
        self.cids_blacklist_fn = os.path.join(self.data_dir, 'misc', 'cids_blacklist.jsonl')
        self.cids_filtered_synonyms_fn = os.path.join(self.data_dir, 'misc', 'cids_filtered_synonyms.jsonl')

        self.reactions_details_fn = os.path.join(self.data_dir, 'reactions_details', 'reactions_details.jsonl')
        self.reactions_details_ord_fn = os.path.join(self.data_dir, 'reactions_details', 'reactions_details_ord.jsonl')
        self.reactions_details_llm_fn = os.path.join(self.data_dir, 'reactions_details', 'reactions_details_llm.jsonl')
        self.reactions_descriptions_fn = os.path.join(self.data_dir, 'reactions_details', 'reactions_descriptions.jsonl')

        self.unmapped_names_fn = os.path.join(self.data_dir, "unmapped_names.jsonl")
        self.chem_names_blacklisted_fn = os.path.join(self.data_dir, "unmapped_names_blacklisted.txt")
        self.unbalancing_cids_fn = os.path.join(self.data_dir, "unbalancing_cids.jsonl")
        self.unmapped_smiles_fn = os.path.join(self.data_dir, 'unmapped_smiles.jsonl')
        self.chem_smiles_blacklisted_fn = os.path.join(self.data_dir, 'unmapped_smiles_blacklisted.txt')

        self.chems_thermo_xtb_fn = os.path.join(self.data_dir, 'thermo', 'chems_thermo_xtb.jsonl')
        self.reactions_thermo_xtb_fn = os.path.join(self.data_dir, 'thermo', 'reactions_thermo_xtb.jsonl')

        self.crc_inorganic_constants_fn = os.path.join(self.data_dir, 'crc_handbook', 'inorganic_constants.jsonl')
        self.crc_organic_constants_fn = os.path.join(self.data_dir, 'crc_handbook', 'organic_constants.jsonl')
        self.crc_flammability_fn = os.path.join(self.data_dir, 'crc_handbook', 'flammability.jsonl')

        self.crc_unmapped_names_fn = os.path.join(self.data_dir, 'crc_handbook', 'crc_unmapped_names.txt')

        self._file_sorting_prefs = {
            self.reactions_parsed_fn: 'rid',
            self.reactions_parsed_llm_fn: 'rid',
            self.reactions_parsed_ord_fn: 'rid',

            self.chems_descriptions_fn: 'cid',
            self.chems_fn: 'complexity',
            self.chems_categories_fn: 'cid',
            self.wiki_chems_fn: 'cid',
            self.hazards_chems_fn: 'cid',
            self.chems_edges_fn: 'eid',
            self.elements_fn: None,

            self.categories_fn: 'code',
            self.cids_blacklist_fn: 'cid',
            self.cids_filtered_synonyms_fn: 'cid',

            self.reactions_details_fn: 'rid',
            self.reactions_details_ord_fn: 'rid',
            self.reactions_details_llm_fn: 'rid',
            self.reactions_descriptions_fn: 'rid',

            self.unmapped_names_fn: ('count', True),
            self.unmapped_smiles_fn: ('count', True),
            self.unbalancing_cids_fn: ('count', True),

            self.chems_thermo_xtb_fn: 'cid',
            self.reactions_thermo_xtb_fn: 'rid',

            self.crc_flammability_fn: 'name',
            self.crc_inorganic_constants_fn: 'name',
            self.crc_organic_constants_fn: 'name'
        }

        self.console = Console()
        traceback.install(console=self.console)
    

    def _load_jsonl(self, filename):
        if not os.path.exists(filename):
            return []

        with open(filename) as f:
            content = f.read().strip()
            if not content:
                return []

            return [json.loads(x) for x in content.split('\n')]
    
    def _write_jsonl(self, entries, filename, backup=True, force_natural_order=False):
        staged_entries = entries

        if not force_natural_order:
            if filename not in self._file_sorting_prefs:
                raise Exception(f"Sorting preferences aren't assigned for file '{filename}'")
            
            sorting_prefs = self._file_sorting_prefs[filename]
            if sorting_prefs is not None:
                if isinstance(sorting_prefs, tuple) and len(sorting_prefs) == 2:
                    sorting_field = sorting_prefs[0]
                    sorting_reverse = sorting_prefs[1]
                elif isinstance(sorting_prefs, str):
                    sorting_field = sorting_prefs
                    sorting_reverse = False
                else:
                    raise Exception(f"Invalid format of sorting preferences for file '{filename}': {str(sorting_prefs)}")

                staged_entries = sorted(entries, key=lambda x: x[sorting_field], reverse=sorting_reverse)


        if os.path.exists(filename) and backup:
            shutil.copy(filename, f"{filename}.backup")

        with open(filename, 'w') as f:
            for entry in staged_entries:
                f.write(json.dumps(entry) + '\n')


    def log(self, message=""):
        with self.print_lock:
            self.console.print(message)
    
