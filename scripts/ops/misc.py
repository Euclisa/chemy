import json
import os

from rdkit import Chem
from rdkit.Chem import Draw


class MiscOps:
    def __init__(self, data_dir, compounds, reactions, store, logger, properties):
        self.data_dir = data_dir
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.properties = properties

        self.unbalancing_cids_fn = os.path.join(self.data_dir, "unbalancing_cids.jsonl")
        store.register_sorting(self.unbalancing_cids_fn, ('count', True))

    def get_background_substances(self, k):
        reagents_cids_count = {cid: 0 for cid in self.compounds.cid_chem_map}
        for react in self.properties.parsed_reactions:
            for cid in [entry['cid'] for entry in react['reagents']]:
                reagents_cids_count[cid] += 1

        background_cids = sorted(
            reagents_cids_count.keys(),
            key=lambda cid: reagents_cids_count[cid],
            reverse=True,
        )[:k]

        with open(self.compounds.background_cids_fn, 'w') as f:
            f.write(json.dumps(background_cids, indent=2))

    def get_commonness_chems_sorting(self):
        occurence = self.properties.get_chems_reactions_occurence()
        return sorted(
            [chem['cid'] for chem in self.compounds.chems],
            key=lambda cid: occurence[cid],
        )

    def get_complexity_chems_sorting(self):
        return sorted(
            [chem['cid'] for chem in self.compounds.chems],
            key=lambda cid: self.compounds.cid_chem_map[cid]['bertz_complexity'],
        )

    def get_curiosity_chems_sorting(self):
        curiosity = self.store.load_jsonl(self.properties.chems_curiosity_fn)
        cid_to_curiosity = {entry['cid']: entry['curiosity'] for entry in curiosity}
        return sorted(
            [chem['cid'] for chem in self.compounds.chems],
            key=lambda cid: cid_to_curiosity.get(cid, 0),
        )

    def extract_radicals_list(self, out_fn):
        radicals = []
        for chem in self.compounds.chems:
            all_names = chem['cmpdsynonym'] + [chem['cmpdname']]
            for name in all_names:
                if '(.)' in name:
                    radicals.append({'cid': chem['cid'], 'name': chem['cmpdname']})
                    break

        self.store.write_jsonl(radicals, out_fn, backup=False)

    def get_chems_cids(self, out_fn):
        with open(out_fn, 'w') as f:
            cids_str = '\n'.join(str(chem['cid']) for chem in self.compounds.chems)
            f.write(cids_str)

    def find_unbalancing_chems(self):
        unbalanced_cids = set()
        unbalanced_cids_cnt = dict()
        balanced_cids = set()
        cid_to_name = dict()

        for entry in self.properties.parsed_reactions:
            cids = set()

            for chem in entry['reagents']:
                cid = chem['cid']
                cid_to_name[cid] = chem['original_name']
                cids.add(cid)

            for chem in entry['products']:
                cid = chem['cid']
                cid_to_name[cid] = chem['original_name']
                cids.add(cid)

            if self.reactions.is_react_balanced(entry):
                balanced_cids.update(cids)
            else:
                unbalanced_cids.update(cids)
                for cid in cids:
                    unbalanced_cids_cnt[cid] = unbalanced_cids_cnt.get(cid, 0) + 1

        res_cids = unbalanced_cids - balanced_cids
        res_entries = [
            {'cid': cid, 'name': cid_to_name[cid], 'count': unbalanced_cids_cnt[cid]}
            for cid in res_cids
        ]

        self.store.write_jsonl(res_entries, self.unbalancing_cids_fn)

    def generate_chems_structures_svg(self, force=False):
        for chem in self.logger.track(self.compounds.chems, "Generating structures..."):
            cid = chem['cid']
            name = chem['cmpdname']
            smiles = chem['smiles']

            svg_fn = os.path.join(self.compounds.structures_dir, f"{cid}.svg")
            if os.path.exists(svg_fn) and not force:
                continue

            try:
                mol = Chem.MolFromSmiles(smiles)
                drawer = Draw.MolDraw2DSVG(300, 300)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
            except Exception:
                self.logger.log(f"Failed to generate structure for '{name}'")
                continue

            with open(svg_fn, "w") as f:
                f.write(svg)

    def show_bertz_complexity_diff(self):
        ct_diff_sum = 0
        ct_sum = 0
        ct_bertz_sum = 0
        ct_cnt = 0
        ct_max_cnt = 500
        for i, chem in enumerate(self.compounds.chems):
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
                print(
                    f"{i-ct_max_cnt+1}-{i+1}: av_ct_diff = {av_ct_diff}; "
                    f"av_ct = {av_ct}; av_ct_bertz = {av_ct_bertz}"
                )
