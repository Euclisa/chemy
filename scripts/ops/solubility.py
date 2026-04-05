import os


class SolubilityGenerator:
    def __init__(self, data_dir, compounds, store, logger, crc, bigsol):
        self.compounds = compounds
        self.store = store
        self.logger = logger
        self.crc = crc
        self.bigsol = bigsol

        chems_properties_dir = os.path.join(data_dir, 'chems_properties')
        self.chems_solubility_fn = os.path.join(chems_properties_dir, 'chems_solubility.jsonl')

        store.register_sorting(self.chems_solubility_fn, 'cid')

    def generate_chems_solubility(self):
        inorganic_crc = self.store.load_jsonl(self.crc.crc_inorganic_constants_fn)
        organic_crc = self.store.load_jsonl(self.crc.crc_organic_constants_fn)
        bigsol = self.store.load_jsonl(self.bigsol.big_sol_parsed_fn)

        result = dict()
        for entry in inorganic_crc + organic_crc:
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

                def build_sol_str(item):
                    return f"{item[1]:.1f} g/L ({item[0]} K)"

                if len(sol_temps) > 1:
                    sol_temps = [sol_temps[0], sol_temps[-1]]

                sol_str = '; '.join(build_sol_str(item) for item in sol_temps)
                result[cid][sol_cid] = sol_str

        result_entries = []
        for cid, sols in result.items():
            sols_list = []
            for sol_cid, sol_str in sols.items():
                sol_name = self.compounds.cid_chem_map[sol_cid]['cmpdname']
                sols_list.append({'solvent_cid': sol_cid, 'solvent_name': sol_name, 'value': sol_str})
            result_entries.append({'cid': cid, 'solubility': sols_list})

        self.store.write_jsonl(result_entries, self.chems_solubility_fn)
