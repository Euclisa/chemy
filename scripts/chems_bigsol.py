import os
import csv

from chems_properties import ChemsProperties


class ChemsBigSol(ChemsProperties):

    """
    Krasnov, L., Malikov, D., Kiseleva, M. et al.
    BigSolDB 2.0, dataset of solubility values for organic compounds in different solvents at various temperatures.
    Sci Data 12, 1236 (2025). https://doi.org/10.1038/s41597-025-05559-8
    """

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.big_sol_fn = os.path.join(self.chems_properties_assets_dir, 'bigsol.csv')
        self.big_sol_parsed_fn = os.path.join(self.chems_properties_dir, 'big_sol.jsonl')

        self._file_sorting_prefs[self.big_sol_parsed_fn] = 'cid'
    

    def parse_bigsol(self):
        parsed_entries_dict = dict()
        with open(self.big_sol_fn) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    cid = int(row['PubChem_CID'])
                    if cid not in self.cid_chem_map:
                        continue

                    solvent_cid = self.smiles_cid_map.get(row['SMILES_Solvent'])
                    if solvent_cid is None:
                        continue
                    solvent_cid = int(solvent_cid)

                    chem = self.cid_chem_map[cid]
                    solubility = float(row['Solubility(mol/L)']) * chem['mw']
                    temp = float(row['Temperature_K'])

                    parsed_entries_dict.setdefault(cid, dict()).setdefault(solvent_cid, dict())[temp] = solubility
                
                except Exception:
                    continue
        
        parsed_entries = []
        for cid, sols in parsed_entries_dict.items():
            name = self.cid_chem_map[cid]['cmpdname']
            parsed_entries.append({'cid': cid, 'name': name, 'solubility': sols})

        self._write_jsonl(parsed_entries, self.big_sol_parsed_fn)


if __name__ == "__main__":
    bigsol = ChemsBigSol('data/')
    bigsol.parse_bigsol()

                
