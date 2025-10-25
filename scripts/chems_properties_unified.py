import os
import numpy as np

from chems_solubility import ChemsSolubility
from chems_crc import ChemsCRC
from chems_hazards import ChemsHazards


class ChemsPropertiesUnified(ChemsSolubility, ChemsCRC, ChemsHazards):
    
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_properties_compiled_fn = os.path.join(self.chems_properties_dir, 'chems_properties.jsonl')
        self.chems_curiosity_fn = os.path.join(self.chems_properties_dir, 'chems_curiosity.jsonl')
        
        self._file_sorting_prefs[self.chems_properties_compiled_fn] = 'cid'
        self._file_sorting_prefs[self.chems_curiosity_fn] = ('curiosity', True)

        self.property_max_synonyms = 7
        self.property_max_cas = 5
    

    def compile_chems_properties(self):
        solubility = self._load_jsonl(self.chems_solubility_fn)
        crc_inorganic = self._load_jsonl(self.crc_inorganic_constants_fn)
        crc_organic = self._load_jsonl(self.crc_organic_constants_fn)
        crc_flammability = self._load_jsonl(self.crc_flammability_fn)

        result = dict()

        def add_property(cid, property, value):
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
            add_property(cid, 'Complexity', chem['complexity'])
        
        
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
        
        for entry in self._rich_track(crc_inorganic+crc_organic, 'Adding CRC constants', transient=True):
            cid = entry['cid']
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
        self._write_jsonl(result_entries, self.chems_properties_compiled_fn)
    


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




if __name__ == "__main__":
    props_compile = ChemsPropertiesUnified('data/')
    props_compile.generate_curiosity_index()
            
