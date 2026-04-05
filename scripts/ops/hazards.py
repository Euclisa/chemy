import os


class HazardAssembler:
    LLM_HAZARD_CATEGORIES = {
        'explosive',
        'acute_toxic',
        'flammable',
        'oxidizer',
        'corrosive',
        'serious_health_hazard',
        'environment_hazard',
    }

    def __init__(self, data_dir, compounds, store, logger):
        self.compounds = compounds
        self.store = store
        self.logger = logger

        chems_properties_dir = os.path.join(data_dir, 'chems_properties')
        chems_properties_llm_dir = os.path.join(chems_properties_dir, 'llm')

        self.chems_hazards_fn = os.path.join(chems_properties_dir, 'chems_hazards.jsonl')
        self.chems_hazard_categories_llm_fn = os.path.join(
            chems_properties_llm_dir,
            'chems_hazard_categories_llm.jsonl',
        )
        self.chems_nfpa_llm_fn = os.path.join(chems_properties_llm_dir, 'chems_nfpa_llm.jsonl')

        store.register_sorting(self.chems_hazards_fn, 'cid')
        store.register_sorting(self.compounds.chems_hazards_wiki_fn, 'cid')

    def assemble_chems_hazards(self):
        wiki_hazards = self.store.load_jsonl(self.compounds.chems_hazards_wiki_fn)
        llm_hazards_cats = self.store.load_jsonl(self.chems_hazard_categories_llm_fn)
        llm_nfpa = self.store.load_jsonl(self.chems_nfpa_llm_fn)

        cid_hazards = dict()
        for entry in wiki_hazards:
            source = 'wikipedia'
            nfpa = {'value': entry['nfpa'], 'source': source}
            picts = {'value': entry['pictograms'], 'source': source} if entry['pictograms'] else None
            cid_hazards[entry['cid']] = {'nfpa': nfpa, 'pictograms': picts}

        llm_to_ghs = {
            'explosive': 'GHS01',
            'flammable': 'GHS02',
            'oxidizer': 'GHS03',
            'corrosive': 'GHS05',
            'acute_toxic': 'GHS06',
            'serious_health_hazard': 'GHS08',
            'environment_hazard': 'GHS09',
        }

        for entry in llm_hazards_cats:
            cid = entry['cid']
            cid_hazards.setdefault(cid, dict())
            if 'pictograms' not in cid_hazards[cid] or cid_hazards[cid]['pictograms'] is None:
                picts = None
                if entry['categories']:
                    picts = {
                        'value': [llm_to_ghs[item] for item in entry['categories']],
                        'source': entry['source'],
                    }
                cid_hazards[cid]['pictograms'] = picts

        for entry in llm_nfpa:
            cid = entry['cid']
            cid_hazards.setdefault(cid, dict())
            if 'nfpa' not in cid_hazards[cid]:
                cid_hazards[cid]['nfpa'] = {'value': entry['nfpa'], 'source': entry['source']}

        chems_hazards = []
        for cid, hazards in cid_hazards.items():
            chems_hazards.append(
                {
                    'cid': cid,
                    'nfpa': hazards['nfpa'] if 'nfpa' in hazards else None,
                    'pictograms': hazards['pictograms'] if 'pictograms' in hazards else None,
                }
            )

        self.store.write_jsonl(chems_hazards, self.chems_hazards_fn)
