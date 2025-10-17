import json
import os

from rdkit import Chem
from rdkit.Chem import Draw

from chems_parse_reactions import ChemsParseReactions


class ChemsMisc(ChemsParseReactions):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.unbalancing_cids_fn = os.path.join(self.data_dir, "unbalancing_cids.jsonl")

        self._file_sorting_prefs[self.unbalancing_cids_fn] = ('count', True)
    

    def get_background_substances(self, k):
        reactions = self._load_jsonl(self.reactions_parsed_fn)
        
        reagents_cids_count = dict()
        for cid in self.cid_chem_map:
            reagents_cids_count[cid] = 0
        
        for react in reactions:
            for cid in [x['cid'] for x in react['reagents']]:
                reagents_cids_count[cid] += 1
        
        background_cids = sorted(list(reagents_cids_count.keys()), key=lambda x: reagents_cids_count[x], reverse=True)[:k]

        with open(self.background_cids_fn, 'w') as f:
            f.write(json.dumps(background_cids, indent=2))
    

    def get_commonnes_chems_sorting(self):
        with open(self.reactions_parsed_fn) as f:
            reactions = [json.loads(x) for x in f.read().strip().split('\n')]
        
        reagents_cids_count = dict()
        for cid in self.cid_chem_map:
            reagents_cids_count[cid] = 0
        
        for react in reactions:
            for cid in [x['cid'] for x in react['reagents']]:
                reagents_cids_count[cid] += 1
        
        sorted_cids = sorted(list(reagents_cids_count.keys()), key=lambda x: reagents_cids_count[x], reverse=True)

        with open(self.commonness_sorted_cids_fn, 'w') as f:
            f.write(json.dumps(sorted_cids, indent=2))

    
    

    def extract_radicals_list(self, out_fn):
        chems = self._load_jsonl(self.chems_fn)
        radicals = []
        for chem in chems:
            all_names = chem['cmpdsynonym'] + [chem['cmpdname']]
            for name in all_names:
                if '(.)' in name:
                    radicals.append({'cid': chem['cid'], 'name': chem['cmpdname']})
                    break
        
        self._write_jsonl(radicals, out_fn, backup=False)
    

    def get_chems_cids(self, out_fn):
        with open(out_fn, 'w') as f:
            cids_str = '\n'.join(str(chem['cid']) for chem in self.chems)
            f.write(cids_str)
    

    def find_unbalancing_chems(self):
        unbalanced_cids = set()
        unbalanced_cids_cnt = dict()
        balanced_cids = set()
        cid_to_name = dict()
        with open(self.reactions_parsed_fn) as f:
            for line in f:
                entry = json.loads(line)
                cids = set()

                for chem in entry['reagents']:
                    cid = chem['cid']
                    name = chem['original_name']
                    cid_to_name[cid] = name
                    cids.add(cid)
                
                for chem in entry['products']:
                    cid = chem['cid']
                    name = chem['original_name']
                    cid_to_name[cid] = name
                    cids.add(cid)

                if entry['balanced']:
                    balanced_cids.update(cids)
                else:
                    unbalanced_cids.update(cids)
                    for cid in cids:
                        if cid not in unbalanced_cids_cnt:
                            unbalanced_cids_cnt[cid] = 0
                        unbalanced_cids_cnt[cid] += 1
        
        res_cids = unbalanced_cids - balanced_cids
        res_entries = []
        for cid in res_cids:
            res_entries.append({'cid': cid, 'name': cid_to_name[cid], 'count': unbalanced_cids_cnt[cid]})

        self._write_jsonl(res_entries, self.unbalancing_cids_fn)
    

    def generate_chems_structures_svg(self, force=False):
        for chem in self._rich_track(self.chems, "Generating structures..."):
            cid = chem['cid']
            name = chem['cmpdname']
            smiles = chem['smiles']

            svg_fn = os.path.join(self.structures_dir, f"{cid}.svg")
            if os.path.exists(svg_fn) and not force:
                continue

            try:
                mol = Chem.MolFromSmiles(smiles)

                drawer = Draw.MolDraw2DSVG(300, 300)  # width, height
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
            except Exception:
                self.log(f"Failed to generate structure for '{name}'")
                continue

            with open(os.path.join(self.structures_dir, f"{cid}.svg"), "w") as f:
                f.write(svg)

    

    def show_bertz_complexity_diff(self):        
        ct_diff_sum = 0
        ct_sum = 0
        ct_bertz_sum = 0
        ct_cnt = 0
        ct_max_cnt = 500
        for i, chem in enumerate(self.chems):
            ct = chem['complexity']
            bertz_ct = chem['bertz_complexity']
            diff = ct - bertz_ct
            ct_diff_sum += diff
            ct_sum += ct
            ct_bertz_sum += bertz_ct
            ct_cnt += 1
            if ct_cnt == ct_max_cnt:
                av_ct_diff = ct_diff_sum / ct_max_cnt
                av_ct = ct_sum / ct_max_cnt
                av_ct_bertz = ct_bertz_sum / ct_max_cnt
                ct_cnt = ct_diff_sum = ct_sum = ct_bertz_sum = 0
                print(f"{i-ct_max_cnt+1}-{i+1}: av_ct_diff = {av_ct_diff}; av_ct = {av_ct}; av_ct_bertz = {av_ct_bertz}")


if __name__ == "__main__":
    chems_misc = ChemsMisc('data/')
    chems_misc.get_commonnes_chems_sorting()