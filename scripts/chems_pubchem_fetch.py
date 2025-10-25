import pubchempy as pcp
import os
from rdkit import Chem
import json
import periodictable
import requests

from chems_properties import ChemsProperties



class ChemsFetchPubchem(ChemsProperties):

    def __init__(self, data_dir):
        super().__init__(data_dir)
    

    def __build_chem_from_fetched(self, fetched):
        chem = {
            'cid': fetched.cid,
            'cmpdname': None,
            'cmpdsynonym': fetched.synonyms,
            'mf': fetched.molecular_formula,
            'mw': fetched.molecular_weight,
            'charge': fetched.charge,
            'smiles': fetched.smiles,
            'inchi': fetched.inchi,
            'inchikey': fetched.inchikey,
            'complexity': fetched.complexity,
            'iupacname': fetched.iupac_name
        }



    def _get_inchikey_from_smiles(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            return Chem.MolToInchiKey(mol)
        except Exception:
            return None

    def fetch_smiles_from_pubchem(self, smiles_fn):
        smiles_list = self._load_jsonl(smiles_fn)
        
        blacklist = set()
        if os.path.exists(self.chem_smiles_blacklisted_fn):
            with open(self.chem_smiles_blacklisted_fn) as f:
                blacklist = set(f.read().strip().split('\n'))
        
        with open(self.chems_fn, 'a') as f_out:
            for entry in smiles_list:
                chem_smiles = entry['smiles']

                if chem_smiles in blacklist:
                    continue
    
                inchikey = self._get_inchikey_from_smiles(chem_smiles)
                if inchikey is None or inchikey in self.inchikey_cid_map:
                    continue

                try:
                    fetched_chems = pcp.get_compounds(chem_smiles, 'smiles')

                    if not fetched_chems:
                        self.log(f"Failed to fetch pubchem data for smiles: '{chem_smiles}'")
                        continue
                    chem = fetched_chems[0]
                    
                    mol_original = Chem.MolFromSmiles(chem_smiles)
                    mol_fetched = Chem.MolFromSmiles(chem.smiles)
                    inchi_original = Chem.MolToInchi(mol_original)
                    inchi_fetched = Chem.MolToInchi(mol_fetched)

                    if inchi_original != inchi_fetched:
                        self.log(f"Fetched substance's inchi doesn't match with original for '{chem_smiles}'. ({inchi_original} != {inchi_fetched})")
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

                except Exception as e:
                    self.log(f"Exception during fetching: {e}")
                    continue

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.log(f"Fetched '{chem_pc_data['cmpdname']}' for unmapped smiles '{chem_smiles}'")
    

    def fetch_chems_cids_from_pubchem_f(self, cids_fn):
        with open(cids_fn) as f:
            cids = [int(x) for x in f.read().strip().split('\n')]
        
        self.fetch_chems_cids_from_pubchem(*cids)
    

    def fetch_names_from_pubchem(self, names_fn):
        names_list = self._load_jsonl(names_fn)
        
        blacklist = set()
        if os.path.exists(self.chem_names_blacklisted_fn):
            with open(self.chem_names_blacklisted_fn) as f:
                blacklist = set(map(lambda x: self._normalize_chem_name(x), f.read().strip().split('\n')))
        
        with open(self.chems_fn, 'a') as f_out, open(self.chem_names_blacklisted_fn, 'a') as f_out_black:
            for entry in names_list:
                chem_name_norm, chem_name, cnt = entry['norm_name'], entry['name'], entry['count']
                if chem_name_norm in self.name_cid_map:
                    continue
                if chem_name_norm in blacklist:
                    continue

                fetched_chems = pcp.get_compounds(chem_name, 'name')
                if not fetched_chems:
                    f_out_black.write(chem_name + '\n')
                    f_out_black.flush()
                    self.log(f"Failed to fetch pubchem data for unmapped name '{chem_name}'")
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

                if chem.inchikey in self.inchikey_cid_map:
                    curr_cid = self.inchikey_cid_map[chem.inchikey]
                    self.log(f"Compound with name '{chem_name}' has equivalent compound in chem file (CID: {curr_cid})")
                    continue

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.log(f"Fetched '{chem_pc_data['cmpdname']}' for name '{chem_name}'")
        
    
    def fetch_chems_cids_from_pubchem(self, cids):
        with open(self.chems_fn, 'a') as f_out:
            for cid in cids:
                if cid in self.cid_chem_map:
                    self.log(f"CID {cid} is in the black list")
                    continue

                fetched_chems = pcp.get_compounds(cid, 'cid')
                if not fetched_chems:
                    self.log(f"Failed to fetch pubchem data for cid '{cid}'")
                    continue
                
                chem = fetched_chems[0]

                if chem.inchikey in self.inchikey_cid_map:
                    curr_cid = self.inchikey_cid_map[chem.inchikey]
                    self.log(f"Compound with CID {cid} has equivalent compound in chem file (CID: {curr_cid})")
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

                f_out.write(json.dumps(chem_pc_data) + '\n')
                f_out.flush()
                self.log(f"Fetched '{chem_pc_data['cmpdname']}' for CID '{cid}'")
    

    def fetch_elements(self):
        with open(self.elements_fn, 'w') as f:
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
                if cid not in self.cid_chem_map:
                    print(f"'{name}' not in chems file. Skipping...")
                    continue

                chem = self.cid_chem_map[cid]
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
        i = 0
        total_incoming = 0
        unique_inchikeys_chems = dict() if override else self._process_chems()
        initial_chems_num = len(unique_inchikeys_chems)
        with open(dump_fn) as f_in:
            for line in f_in:
                for field in fields_to_keep:
                    if f'"{field}":' in line:
                        if field == 'cid':
                            if curr_chem is not None:
                                total_incoming += 1
                                if 'cmpdsynonym' in curr_chem:
                                    curr_chem['charge'] = int(curr_chem['charge'])
                                    curr_chem['cid'] = int(curr_chem['cid'])
                                    curr_chem['mw'] = float(curr_chem['mw'])
                                    curr_chem['complexity'] = float(curr_chem['complexity'])
                                    if not isinstance(curr_chem['cmpdsynonym'], list):
                                        curr_chem['cmpdsynonym'] = [curr_chem['cmpdsynonym']]

                                    self.__process_chem_single(curr_chem, unique_inchikeys_chems, force=True)

                                    i += 1
                                    if i % 1000 == 0:
                                        print(i)

                            curr_chem = dict()
                        curr_chem[field] = json.loads('{' + line.strip().strip(',') + '}')[field]
        
        unique_chems = list(unique_inchikeys_chems.values())
        self._update_chems(unique_chems)

        print(f"Initial chems num: {initial_chems_num}; Incoming chems processed: {total_incoming}; Total written: {len(unique_chems)}")