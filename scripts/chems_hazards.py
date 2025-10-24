import os

from chems_pubchem_parse import ChemsParsePubchem
from chems_llm_parse import ChemsLLMParse



class ChemsHazards(ChemsLLMParse):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_hazards_fn = os.path.join(self.data_dir, 'chems', 'chems_hazards.jsonl')

        self._file_sorting_prefs[self.chems_hazards_fn] = 'cid'
    

    def assemble_chems_hazards(self):
        wiki_hazards = self._load_jsonl(self.chems_hazards_wiki_fn)
        llm_hazards_cats = self._load_jsonl(self.chems_hazard_categories_llm_fn)
        llm_nfpa = self._load_jsonl(self.chems_nfpa_llm_fn)

        cid_hazards = dict()
        for entry in wiki_hazards:
            source = 'wikipedia'
            nfpa = {'value': entry['nfpa'], 'source': source}
            picts = {'value': entry['pictograms'], 'source': source}
            cid_hazards[entry['cid']] = {'nfpa': nfpa, 'pictograms': picts}
        

        llm_to_ghs = {
            'explosive': 'GHS01',
            'flammable': 'GHS02',
            'oxidizer': 'GHS03',
            'corrosive': 'GHS05',
            'acute_toxic': 'GHS06',
            'serious_health_hazard': 'GHS08',
            'environment_hazard': 'GHS09'
        }
        
        for entry in llm_hazards_cats:
            cid = entry['cid']
            if cid in cid_hazards:
                continue
            
            picts = {'value': [llm_to_ghs[x] for x in entry['categories']], 'source': entry['source']}
            cid_hazards[cid] = {'pictograms': picts}
        
        for entry in llm_nfpa:
            cid = entry['cid']
            cid_hazards.setdefault(cid, dict())
            if 'nfpa' not in cid_hazards[cid]:
                nfpa = {'value': entry['nfpa'], 'source': entry['source']}
                cid_hazards[cid]['nfpa'] = nfpa
        

        chems_hazards = []
        for cid, hazards in cid_hazards.items():
            entry = {'cid': cid}
            entry['nfpa'] = hazards['nfpa'] if 'nfpa' in hazards else None
            entry['pictograms'] = hazards['pictograms'] if 'pictograms' in hazards else None
            chems_hazards.append(entry)

        self._write_jsonl(chems_hazards, self.chems_hazards_fn)   


if __name__ == "__main__":
    chems_hazards = ChemsHazards('data/')
    chems_hazards.assemble_chems_hazards()