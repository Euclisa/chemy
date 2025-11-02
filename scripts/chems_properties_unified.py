import os
import numpy as np

from functools import cached_property

from chems_solubility import ChemsSolubility
from chems_crc import ChemsCRC
from chems_hazards import ChemsHazards
from chems_ord_parse import ChemsOrdParse


class ChemsPropertiesUnified(ChemsOrdParse, ChemsSolubility, ChemsCRC, ChemsHazards):
    
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_properties_compiled_fn = os.path.join(self.chems_properties_dir, 'chems_properties.jsonl')
        self.chems_curiosity_fn = os.path.join(self.chems_properties_dir, 'chems_curiosity.jsonl')
        
        self._file_sorting_prefs[self.chems_properties_compiled_fn] = 'cid'
        self._file_sorting_prefs[self.chems_curiosity_fn] = ('curiosity', True)

        self.property_max_synonyms = 7
        self.property_max_cas = 5
    

    @cached_property
    def parsed_reactions(self):
        ord_reactions = self._load_jsonl(self.reactions_parsed_ord_fn)
        llm_reactions = self._load_jsonl(self.reactions_parsed_llm_fn)
        parsed_reactions = {react['rid']: react for react in llm_reactions + ord_reactions}   # order matters
        return tuple(parsed_reactions.values())
    
    @cached_property
    def parsed_reactions_balanced(self):
        return tuple(react for react in self.parsed_reactions if react['rid'] in self.reactions_balance)
    

    @cached_property
    def reactions_details(self):
        ord_details = self._load_jsonl(self.reactions_details_ord_fn)
        llm_details = self._load_jsonl(self.reactions_details_llm_fn)
        parsed_reactions = {entry['rid']: entry for entry in llm_details + ord_details}   # order matters
        return tuple(parsed_reactions.values())
    

    def _write_parsed_reactions(self, reactions):
        ord_reactions = []
        llm_reactions = []
        for react in reactions:
            if react['source'] == self.ord_source:
                ord_reactions.append(react)
            else:
                llm_reactions.append(react)
        
        self._write_jsonl(ord_reactions, self.reactions_parsed_ord_fn)
        self._write_jsonl(llm_reactions, self.reactions_parsed_llm_fn)

    
    def merge_parsed_files(self, out_fn, *parsed_reactions_files):
        rid_reaction = dict()
        total_reactions = 0
        for fn in parsed_reactions_files:
            reactions = self._load_jsonl(fn)
            
            total_reactions += len(reactions)
            
            for react in reactions:
                rid = react['rid']
                if rid not in rid_reaction:
                    rid_reaction[rid] = react
                else:
                    old_react = rid_reaction[rid]
                    old_source = old_react.get('source')
                    new_source = react.get('source')
                    new_source_priority = self.sources_priority.get(new_source, -1)
                    old_source_priority = self.sources_priority.get(old_source, -1)
                    if new_source_priority > old_source_priority:
                        rid_reaction[rid] = react

        reactions_res = list(rid_reaction.values())
        self._write_jsonl(reactions_res, out_fn, backup=False)
    

    def generate_edges(self):        
        edge_reaction_id_map = dict()
        for react in self._rich_track(self.parsed_reactions, "Generating edges..."):
            react_id = react['rid']
            for r in react['reagents']:
                r_cid = r['cid']
                for p in react['products']:
                    p_cid = p['cid']
                    edge = (r_cid, p_cid)
                    if edge not in edge_reaction_id_map:
                        edge_reaction_id_map[edge] = []
                    edge_reaction_id_map[edge].append(react_id)
        
        def get_eid(edge):
            source = edge[0]
            target = edge[1]
            return ((source+target)*target) % 2**64
        
        edges = []
        for edge in edge_reaction_id_map:
            entry = {'eid': get_eid(edge), 'source': edge[0], 'target': edge[1], 'reactions': edge_reaction_id_map[edge]}
            edges.append(entry)

        self._write_jsonl(edges, self.chems_edges_fn, backup=False)
        self.log(f"Generated {len(edge_reaction_id_map)} edges")
    

    def discard_orphaned_unmapped_chems(self):
        connected_cids = set()
        for reaction in self.parsed_reactions:
            all_cids = self._get_all_reaction_cids(reaction)
            for cid in all_cids:
                connected_cids.add(cid)
        
        for entry in self.reactions_details:
            for cmpd in entry['solvents']+entry['catalysts']:
                connected_cids.add(cmpd['cid'])

        
        res_chems = [chem for chem in self.chems_unmapped if chem['cid'] in connected_cids]
        self._update_unmapped_chems(res_chems)
    

    def balance_parsed_reactions(self, reactions_parsed_fn=None):
        if reactions_parsed_fn is None:
            reactions = self._load_jsonl(self.reactions_parsed_llm_fn)+self._load_jsonl(self.reactions_parsed_ord_fn)
        else:
            reactions = self._load_jsonl(reactions_parsed_fn)
        
        balances = self._load_jsonl(self.reactions_balance_fn)
        balanced_rids = set(entry['rid'] for entry in balances)
        
        balanced_cnt = 0
        for react in self._rich_track(reactions, "Balancing parsed reactions"):
            if react['rid'] not in balanced_rids:
                bal = self._balance_reaction(react)
                if bal:
                    balances.append(bal)
                    balanced_cnt += 1

        self._write_jsonl(balances, self.reactions_balance_fn)        
        self.log(f"Processed {len(reactions)}; balanced {balanced_cnt}; total balanced: {len(balances)}")
    

    def _get_chems_reactions_occurence(self, reactions=None):
        if reactions is None:
            reactions = self.parsed_reactions
        else:
            reactions = reactions

        chem_reactions_occurence = dict()
        for chem in self.chems:
            chem_reactions_occurence[chem['cid']] = 0

        for react in reactions:
            all_cids = self._get_all_reaction_cids(react)
            for cid in all_cids:
                chem_reactions_occurence[cid] += 1
        
        return chem_reactions_occurence
    

    def _compile_chems_properties(self):
        solubility = self._load_jsonl(self.chems_solubility_fn)
        crc_inorganic = self._load_jsonl(self.crc_inorganic_constants_fn)
        crc_organic = self._load_jsonl(self.crc_organic_constants_fn)
        crc_flammability = self._load_jsonl(self.crc_flammability_fn)

        result = dict()

        def add_property(cid, property, value):
            if cid not in result:
                return

            result[cid].append({'property': property, 'value': f"{value}"})

        for chem in self._rich_track(self.chems, 'Adding basic properties', transient=True):
            cid = chem['cid']
            result[cid] = []

            add_property(cid, 'PubChem CID', cid)
            add_property(cid, 'Synonyms', ', '.join(chem['cmpdsynonym'][:self.property_max_synonyms]))
            add_property(cid, 'Molecular formula', chem['mf'])
            add_property(cid, 'Molecular weight', f"{chem['mw']} g/mol")
            if chem['wiki']:
                add_property(cid, 'Wikipedia', chem['wiki'])
            if chem['iupacname']:
                add_property(cid, 'IUPAC name', chem['iupacname'])
            if chem['cas']:
                add_property(cid, 'CAS number', ', '.join(chem['cas'][:self.property_max_cas]))
            add_property(cid, 'SMILES', chem['smiles'])
            add_property(cid, 'InChI', chem['inchi'])
            add_property(cid, 'InChI-key', chem['inchikey'])
            add_property(cid, 'Pubchem complexity', chem['complexity'])
            add_property(cid, 'Bertz complexity', chem['bertz_complexity'])
        
        
        for entry in self._rich_track(solubility, 'Adding solubility', transient=True):
            cid = entry['cid']
            sols = entry['solubility']
            value = '\n'.join(f"{s['solvent_name']}: {s['value']}" for s in sols)

            add_property(cid, 'Solubility', value)
        
        def build_crc_str(entry):
            crc_str = str(entry['value'])
            if entry.get('approx'):
                crc_str = '~ ' + crc_str
            if entry.get('decomposes', False):
                crc_str += ' (decomposes)'
            if entry.get('sublimes', False):
                crc_str += ' (sublimes)'
            
            return crc_str
        
        processed_crc_cids = set()
        for entry in self._rich_track(crc_inorganic+crc_organic, 'Adding CRC constants', transient=True):
            cid = entry['cid']
            if cid in processed_crc_cids:
                continue

            if entry['physical_form']:
                add_property(cid, 'Appearence', entry['physical_form'])
            if entry['mp'] and entry['mp']['value'] is not None:
                add_property(cid, 'Melting point', build_crc_str(entry['mp']))
            if entry['bp'] and entry['bp']['value'] is not None:
                add_property(cid, 'Boiling point', build_crc_str(entry['bp']))
            if entry['density']:
                add_property(cid, 'Density', f"{entry['density']} g/cm^3")
            if entry.get('refractive_index', None) is not None:
                add_property(cid, 'Refractive index', entry['refractive_index'])
            
            processed_crc_cids.add(cid)
            

        for entry in self._rich_track(crc_flammability, 'Adding CRC flammability data', transient=True):
            cid = entry['cid']
            if entry['flash_point']:
                fp = build_crc_str(entry['flash_point'])
                add_property(cid, 'Flash point', f"{fp} °C")
            if entry['flash_limits']:
                add_property(cid, 'Flash limits', entry['flash_limits'])
            if entry['ignition_temp']:
                it = build_crc_str(entry['ignition_temp'])
                add_property(cid, 'Ignition temperature', f"{it} °C")
        

        result_entries = [{'cid': cid, 'properties': props} for cid, props in result.items()]

        return result_entries
    


    def generate_curiosity_index(self):
        hazards = self._load_jsonl(self.chems_hazards_fn)
        
        cid_to_curiosity = dict()
        for chem in self.chems:
            cid = chem['cid']
            cid_to_curiosity[cid] = (chem['complexity'] / self.complexity_thr)
        
        for entry in hazards:
            cid = entry['cid']

            if entry['nfpa'] is None:
                nfpa_coeff = 0
            else:
                nfpa_coeff = sum(entry['nfpa']['value'].values()) / 12
            
            if entry['pictograms'] is None:
                pictograms_coeff = 0
            else:
                pictograms_coeff = len(entry['pictograms']['value']) / len(self.LLM_HAZARD_CATEGORIES)
            
            if cid in cid_to_curiosity:
                cid_to_curiosity[cid] += nfpa_coeff + pictograms_coeff
        
        chems_occurence = self._get_chems_reactions_occurence()

        # Convert occurrences to array for percentile calculation
        occur_values = np.array(list(chems_occurence.values()), dtype=float)
        sorted_vals = np.sort(occur_values)
        n = len(sorted_vals)

        def percentile_rank(x):
            # Rank position divided by number of values
            return (np.searchsorted(sorted_vals, x, side='right') / n)

        for cid, occurence in chems_occurence.items():
            if cid in cid_to_curiosity:
                cid_to_curiosity[cid] += 3 * percentile_rank(occurence)
        
        max_curiosity = max(cid_to_curiosity.values())
        curiosity_entries = [{'cid': cid, 'curiosity': curiosity / max_curiosity} for cid, curiosity in cid_to_curiosity.items()]

        self._write_jsonl(curiosity_entries, self.chems_curiosity_fn)
    

    def _get_parsed_reactions_participants_norm_names(self):
        norm_names = set()
        for react in self.parsed_reactions:
            norm_names.update(entry['norm_name'] for entry in react['reagents']+react['products'])
        
        return norm_names




if __name__ == "__main__":
    unified = ChemsPropertiesUnified('data/')
    unified.discard_orphaned_unmapped_chems()
            
