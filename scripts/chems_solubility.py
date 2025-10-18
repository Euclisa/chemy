import os

from chems_crc import ChemsCRC
from chems_bigsol import ChemsBigSol


class ChemsSolubility(ChemsBigSol, ChemsCRC):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_solubility_fn = os.path.join(self.chems_properties_dir, 'chems_solubility.jsonl')
        
        self._file_sorting_prefs[self.chems_solubility_fn] = 'cid'
    

    def generate_chems_solubility(self):
        inorganic_crc = self._load_jsonl(self.crc_inorganic_constants_fn)
        organic_crc = self._load_jsonl(self.crc_organic_constants_fn)
        bigsol = self._load_jsonl(self.big_sol_parsed_fn)

        result = dict()
        for entry in inorganic_crc+organic_crc:
            sols = entry.get('solubility')
            if not sols:
                continue

            cid = entry['cid']
            result.setdefault(cid, dict())
            for sol in sols:
                result[cid][sol['cid']] = sol['solubility']
        

        for entry in bigsol:
            cid = entry['cid']
            result.setdefault(cid, dict())
            sols_list = entry['solubility']
            for sol_cid, sol_temps in sols_list.items():
                sol_cid = int(sol_cid)
                sol_temps = sorted(list(sol_temps.items()))

                def build_sol_str(entry):
                    return f"{entry[1]:.1f} g/L ({entry[0]} K)"

                if len(sol_temps) > 1:
                    sol_temps = [sol_temps[0], sol_temps[-1]]
                
                sol_str = '; '.join(build_sol_str(x) for x in sol_temps)
                
                result[cid][sol_cid] = sol_str
        
        result_entries = []
        for cid, sols in result.items():
            sols_list = []
            for sol_cid, sol_str in sols.items():
                sol_name = self.cid_chem_map[sol_cid]['cmpdname']
                sols_list.append({'solvent_cid': sol_cid, 'solvent_name': sol_name, 'value': sol_str})
            entry = {'cid': cid, 'solubility': sols_list}
            result_entries.append(entry)
        
        self._write_jsonl(result_entries, self.chems_solubility_fn)
            


if __name__ == "__main__":
    sols = ChemsSolubility('data/')
    sols.generate_chems_solubility()