import os
from rdkit import Chem
import json


class PubChemFetcher:

    def __init__(self, data_dir, compounds, store, logger):
        self.data_dir = data_dir
        self.compounds = compounds
        self.store = store
        self.logger = logger

    def _get_inchikey_from_smiles(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            return Chem.MolToInchiKey(mol)
        except Exception:
            return None

    def fetch_smiles_from_pubchem(self, smiles_fn, out_fn):
        import pubchempy as pcp

        smiles_list = self.store.load_jsonl(smiles_fn)

        blacklist = set()
        if os.path.exists(self.compounds.chem_smiles_blacklisted_fn):
            with open(self.compounds.chem_smiles_blacklisted_fn) as f:
                blacklist = set(f.read().strip().split('\n'))

        with open(out_fn, 'a') as f_out:
            for entry in smiles_list:
                chem_smiles = entry['smiles']

                if chem_smiles in blacklist:
                    continue

                inchikey = self._get_inchikey_from_smiles(chem_smiles)
                if inchikey is None or inchikey in self.compounds.inchikey_cid_map:
                    continue

                try:
                    fetched_chems = pcp.get_compounds(chem_smiles, 'smiles')

                    if not fetched_chems:
                        self.logger.log(f"Failed to fetch pubchem data for smiles: '{chem_smiles}'")
                        continue
                    chem = fetched_chems[0]

                    mol_original = Chem.MolFromSmiles(chem_smiles)
                    mol_fetched = Chem.MolFromSmiles(chem.smiles)
                    inchi_original = Chem.MolToInchi(mol_original)
                    inchi_fetched = Chem.MolToInchi(mol_fetched)

                    if inchi_original != inchi_fetched:
                        self.logger.log(f"Fetched substance's inchi doesn't match with original for '{chem_smiles}'. ({inchi_original} != {inchi_fetched})")
                        continue

                    chem_pc_data = {
                        'cid': chem.cid,
                        'cmpdname': chem.iupac_name,
                        'cmpdsynonym': chem.synonyms,
                        'mf': chem.molecular_formula,
                        'mw': chem.molecular_weight,
                        'charge': chem.charge,
                        'smiles': chem.smiles,
                        'inchi': chem.inchi,
                        'inchikey': chem.inchikey,
                        'complexity': chem.complexity
                    }
                    chem_pc_data = self.compounds.process_chem_single(chem_pc_data, force=True)
                    if not chem_pc_data:
                        self.logger.log(f"Fetched compound for smiles '{chem_smiles}' was rejected during processing")
                        continue

                except Exception as e:
                    self.logger.log(f"Exception during fetching: {e}")
                    continue

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.logger.log(f"Fetched '{chem_pc_data['cmpdname']}' for unmapped smiles '{chem_smiles}'")

    def fetch_names_from_pubchem(self, names_fn, out_fn):
        import pubchempy as pcp

        names_list = self.store.load_jsonl(names_fn)

        blacklist = set()
        if os.path.exists(self.compounds.chem_names_blacklisted_fn):
            with open(self.compounds.chem_names_blacklisted_fn) as f:
                blacklist = set(map(lambda x: self.compounds.normalize_chem_name(x), f.read().strip().split('\n')))

        with open(out_fn, 'a') as f_out, open(self.compounds.chem_names_blacklisted_fn, 'a') as f_out_black:
            for entry in names_list:
                chem_name_norm, chem_name, cnt = entry['norm_name'], entry['name'], entry['count']
                if chem_name_norm in self.compounds.name_cids_map:
                    continue
                if chem_name_norm in blacklist:
                    continue

                fetched_chems = pcp.get_compounds(chem_name, 'name')
                if not fetched_chems:
                    f_out_black.write(chem_name + '\n')
                    f_out_black.flush()
                    self.logger.log(f"Failed to fetch pubchem data for unmapped name '{chem_name}'")
                    continue
                chem = fetched_chems[0]
                chem_pc_data = {
                    'cid': chem.cid,
                    'cmpdname': chem.iupac_name,
                    'cmpdsynonym': chem.synonyms,
                    'mf': chem.molecular_formula,
                    'mw': chem.molecular_weight,
                    'charge': chem.charge,
                    'smiles': chem.smiles,
                    'inchi': chem.inchi,
                    'inchikey': chem.inchikey,
                    'complexity': chem.complexity
                }
                chem_pc_data = self.compounds.process_chem_single(chem_pc_data, force=True)
                if not chem_pc_data:
                    self.logger.log(f"Fetched compound for name '{chem_name}' was rejected during processing")
                    continue

                if chem.inchikey in self.compounds.inchikey_cid_map:
                    curr_cid = self.compounds.inchikey_cid_map[chem.inchikey]
                    self.logger.log(f"Compound with name '{chem_name}' has equivalent compound in chem file (CID: {curr_cid})")
                    continue

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.logger.log(f"Fetched '{chem_pc_data['cmpdname']}' for name '{chem_name}'")

    def fetch_chems_cids_from_pubchem(self, cids, out_fn):
        import pubchempy as pcp

        with open(out_fn, 'a') as f_out:
            for cid in cids:
                if cid in self.compounds.cid_chem_map:
                    self.logger.log(f"CID {cid} is in the black list")
                    continue

                fetched_chems = pcp.get_compounds(cid, 'cid')
                if not fetched_chems:
                    self.logger.log(f"Failed to fetch pubchem data for cid '{cid}'")
                    continue

                chem = fetched_chems[0]

                if chem.inchikey in self.compounds.inchikey_cid_map:
                    curr_cid = self.compounds.inchikey_cid_map[chem.inchikey]
                    self.logger.log(f"Compound with CID {cid} has equivalent compound in chem file (CID: {curr_cid})")
                    continue

                chem_pc_data = {
                    'cid': chem.cid,
                    'cmpdname': chem.iupac_name,
                    'cmpdsynonym': chem.synonyms,
                    'mf': chem.molecular_formula,
                    'mw': chem.molecular_weight,
                    'charge': chem.charge,
                    'smiles': chem.smiles,
                    'inchi': chem.inchi,
                    'inchikey': chem.inchikey,
                    'complexity': chem.complexity
                }
                chem_pc_data = self.compounds.process_chem_single(chem_pc_data, force=True)
                if not chem_pc_data:
                    self.logger.log(f"Fetched compound for CID '{cid}' was rejected during processing")
                    continue

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.logger.log(f"Fetched '{chem_pc_data['cmpdname']}' for CID '{cid}'")

    def fetch_elements(self):
        import periodictable
        import requests

        with open(self.compounds.elements_fn, 'w') as f:
            base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/TXT"
            for i, el in enumerate(periodictable.elements):
                if i >= 118:
                    break
                name = el.name
                response = requests.get(base_url.format(name))
                if response.status_code != 200:
                    print(f"Failed to obtain cid for '{name}'")
                    continue
                cid = int(response.text.strip())
                if cid not in self.compounds.cid_chem_map:
                    print(f"'{name}' not in chems file. Skipping...")
                    continue

                chem = self.compounds.cid_chem_map[cid]
                mol = Chem.MolFromInchi(chem['inchi'])
                mol = Chem.AddHs(mol)
                if not mol:
                    print(f"Failed to build mol object for '{name}'")
                    continue
                atom_count = mol.GetNumAtoms()

                entry = {'cid': cid, 'name': name, 'symbol': el.symbol, 'atom_count': atom_count}
                f.write(json.dumps(entry) + '\n')
                f.flush()

                print(f"Fetched '{name}' ({cid}), atom count: {atom_count}")

    def extract_pubchem_dump_to_chems(self, dump_fn, override=False):
        fields_to_keep = ['cid', 'cmpdname', 'cmpdsynonym', 'mf', 'mw', 'complexity', 'smiles', 'inchi', 'inchikey', 'charge', 'annotation', 'iupacname']
        curr_chem = None
        total_incoming = 0
        res_chems = [] if override else list(self.compounds.chems_mapped)
        unique_cids = set(chem['cid'] for chem in res_chems)
        initial_chems_num = len(res_chems)

        def flush_curr_chem(curr_chem_local):
            nonlocal total_incoming

            if curr_chem_local is None or 'cmpdsynonym' not in curr_chem_local:
                return

            total_incoming += 1
            curr_chem_local['charge'] = int(curr_chem_local['charge'])
            cid = curr_chem_local['cid'] = int(curr_chem_local['cid'])
            curr_chem_local['mw'] = float(curr_chem_local['mw'])
            curr_chem_local['complexity'] = float(curr_chem_local['complexity'])
            if not isinstance(curr_chem_local['cmpdsynonym'], list):
                curr_chem_local['cmpdsynonym'] = [curr_chem_local['cmpdsynonym']]

            if cid in unique_cids:
                return

            curr_chem_local = self.compounds.process_chem_single(curr_chem_local, force=True)
            if curr_chem_local:
                res_chems.append(curr_chem_local)
                unique_cids.add(cid)

        with open(dump_fn) as f_in:
            total = self.compounds.count_file_lines(dump_fn)
            for line in self.logger.track(f_in, "Parsing pubchem dump", total=total):
                for field in fields_to_keep:
                    if f'"{field}":' in line:
                        if field == 'cid':
                            flush_curr_chem(curr_chem)

                            curr_chem = dict()
                        curr_chem[field] = json.loads('{' + line.strip().strip(',') + '}')[field]

        flush_curr_chem(curr_chem)
        self.compounds.update_mapped_chems(res_chems)

        print(f"Initial chems num: {initial_chems_num}; Incoming chems processed: {total_incoming}; Total written: {len(res_chems)}")
