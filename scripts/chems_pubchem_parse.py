import json
import re
import shutil
import os

from rdkit import Chem
from rdkit.Chem import AllChem, GraphDescriptors, inchi

from functools import cached_property
import inspect

from rich.progress import track, Progress

from chems_db import ChemsDB


class ChemsParsePubchem(ChemsDB):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_descriptions_fn = os.path.join(self.data_dir, 'chems', 'chems_descriptions.jsonl')
        self.chems_fn = os.path.join(self.data_dir, 'chems', "chems.jsonl")
        self.chems_categories_fn = os.path.join(self.data_dir, 'chems', "chems_categories.jsonl")
        self.wiki_chems_fn = os.path.join(self.data_dir, 'chems', "wiki_chems.jsonl")
        self.hazards_chems_fn = os.path.join(self.data_dir, 'chems', "hazards_chems.jsonl")
        self.chems_edges_fn = os.path.join(self.data_dir, 'chems', 'chems_edges.jsonl')
        self.elements_fn = os.path.join(self.data_dir, 'chems', 'elements.jsonl')

        self.categories_fn = os.path.join(self.data_dir, 'misc', "categories.jsonl")
        self.background_cids_fn = os.path.join(self.data_dir, 'misc', 'background_cids.json')
        self.commonness_sorted_cids_fn = os.path.join(self.data_dir, 'misc', 'commonness_sorted_cids.json')
        self.cids_blacklist_fn = os.path.join(self.data_dir, 'misc', 'cids_blacklist.jsonl')
        self.cids_filtered_synonyms_fn = os.path.join(self.data_dir, 'misc', 'cids_filtered_synonyms.jsonl')

        self._file_sorting_prefs[self.chems_descriptions_fn] = 'cid'
        self._file_sorting_prefs[self.chems_fn] = 'complexity'
        self._file_sorting_prefs[self.chems_categories_fn] = 'cid'
        self._file_sorting_prefs[self.wiki_chems_fn] = 'cid'
        self._file_sorting_prefs[self.hazards_chems_fn] = 'cid'
        self._file_sorting_prefs[self.chems_edges_fn] = 'eid'
        self._file_sorting_prefs[self.elements_fn] = None

        self._file_sorting_prefs[self.categories_fn] = 'code'
        self._file_sorting_prefs[self.cids_blacklist_fn] = 'cid'
        self._file_sorting_prefs[self.cids_filtered_synonyms_fn] = 'cid'

        self.complexity_thr = 550
        self.max_synonyms_thr = 150

        self.CAS_PATTERN = r'\d{2,7}-\d{2}-\d'
    

    @cached_property
    def chems(self):
        return self._load_jsonl(self.chems_fn)
    
    @cached_property
    def cids_blacklist(self):
        return set([x['cid'] for x in self._load_jsonl(self.cids_blacklist_fn)])

    @cached_property
    def cid_chem_map(self):
        return {chem['cid']: chem for chem in self.chems}
    
    @cached_property
    def inchikey_cid_map(self):
        return {chem['inchikey_snone']: chem['cid'] for chem in self.chems}
    
    @cached_property
    def smiles_cid_map(self):
        return {chem['smiles']: chem['cid'] for chem in self.chems}

    @cached_property
    def cid_mf_map(self):
        return {chem['cid']: chem['mf'] for chem in self.chems}
    

    @cached_property
    def cid_wiki_map(self):
        wiki_entries = self._load_jsonl(self.wiki_chems_fn)
        return {entry['cid']: entry['wiki'] for entry in wiki_entries}


    @cached_property
    def cids_filtered_synonyms(self):
        entries = self._load_jsonl(self.cids_filtered_synonyms_fn)
        return {x['cid']: set(x['synonyms']) for x in entries}

    @cached_property
    def name_cid_map(self):
        _name_cid_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            _name_cid_map[self._normalize_chem_name(chem['cmpdname'], is_clean=True)] = cid
            for syn in chem['cmpdsynonym']:
                _name_cid_map[self._normalize_chem_name(syn, is_clean=True)] = cid

        return _name_cid_map

    @cached_property
    def cas_cid_map(self):
        _cas_cid_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            cas_list = chem['cas']
            if cas_list:
                for cas in cas_list:
                    _cas_cid_map[cas] = cid
        
        return _cas_cid_map


    def __clear_runtime_chems_properties(self):
        for name, attr in inspect.getmembers(type(self)):
            if isinstance(attr, cached_property) and name in self.__dict__:
                delattr(self, name)

    def _update_chems(self, new_chems):
        new_chems.sort(key=lambda x: x['complexity'])
        self._write_jsonl(new_chems, self.chems_fn)
        self.__clear_runtime_chems_properties()

    def _update_cids_blacklist(self, cids):
        cids = set(cids)
        unique_new_cids = cids - self.cids_blacklist
        if not unique_new_cids:
            return

        with open(self.cids_blacklist_fn, 'a') as f:
            for cid in unique_new_cids:
                name = self.cid_chem_map.get(cid, {}).get('cmpdname', 'Unknown')
                f.write(json.dumps({'cid': cid, 'name': name}) + '\n')
            f.flush()

        self.cids_blacklist.update(unique_new_cids)
    

    def _count_file_lines(self, filename):
        count = 0
        with open(filename, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                count += chunk.count(b"\n")
        return count
    

    def _rich_track(self, iterator, description, total=None, transient=False, auto_refresh=True):
        if total is None:
            try:
                total = len(iterator)
            except TypeError:
                pass

        with self._rich_progress(transient=transient, auto_refresh=auto_refresh) as progress:
            task = progress.add_task(description=description, total=total)
            for item in iterator:
                yield item
                progress.update(task, advance=1)
                progress.refresh()
    

    def _rich_progress(self, transient=False, auto_refresh=True):
        return Progress(console=self._console, auto_refresh=auto_refresh, transient=transient, refresh_per_second=1)
    

    def _get_mol_fingerprint(self, mol):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        bitstring = fp.ToBitString()
        popcount = sum([int(x) for x in bitstring])

        chunks = [bitstring[i:i+32] for i in range(0, 1024, 32)]
        ints32 = [int(c, 2) - 2**32 if int(c, 2) >= 2**31 else int(c, 2) for c in chunks]

        return {'bits': ints32, 'popcount': popcount}


    def generate_chems_fingerprints(self):
        for chem in self.chems:
            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")

            chem['ECFP4_fp'] = self._get_mol_fingerprint(mol)
        
        self._update_chems(self.chems)
    

    def _extract_chem_cas_numbers(self, chem):
        cas_numbers = map(lambda x: re.sub(r'^cas-', '', x.strip()), chem['cmpdsynonym'])
        cas_numbers = filter(lambda x: re.fullmatch(self.CAS_PATTERN, x), cas_numbers)
        cas_numbers = list(set(cas_numbers))

        return cas_numbers
    

    def extract_chems_cas_numbers(self):
        shutil.copy(self.chems_fn, f"{self.chems_fn}.backup")
        
        cnt = 0
        for chem in self.chems:
            cas_numbers = self._extract_chem_cas_numbers(chem)
            chem['cas'] = cas_numbers
            if len(cas_numbers) != 0:
                cnt += 1

        self._update_chems(self.chems)
        
        print(f"Extracted CAS numbers for {cnt} chems")
    

    def _get_mol_bertz_complexity(self, mol):
        return GraphDescriptors.BertzCT(mol)

    def generate_chems_bertz_complexity(self):        
        for chem in self.chems:
            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")
            chem['bertz_complexity'] =  self._get_mol_bertz_complexity(mol)
        
        self._update_chems(self.chems)
    

    def merge_wiki_chems(self):
        for chem in self.chems:
            chem['wiki'] = self.cid_wiki_map.get(chem['cid'], None)
        
        self._update_chems(self.chems)

    def _get_mol_organic_mark(self, mol):
        try:
            mol = Chem.AddHs(mol)

            atoms = mol.GetAtoms()
            bonds = mol.GetBonds()

            has_carbon = any(atom.GetAtomicNum() == 6 for atom in atoms)

            if not has_carbon:
                return False

            for bond in bonds:
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()

                if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                    return True

                if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 1) or \
                        (atom1.GetAtomicNum() == 1 and atom2.GetAtomicNum() == 6):
                    return True

            return False

        except Exception as e:
            self.log("Error generating organic mark. Assuming True")
            # If complex smiles then likely organic
            return True

    def generate_organic_marks_for_chems(self):
        for chem in self.chems:
            smiles = chem['smiles']
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")
            chem['organic'] = self._get_mol_organic_mark(mol)

        self._update_chems(self.chems)

    def __process_chem_synonyms(self, chem):

        filtered_synonyms = self.cids_filtered_synonyms.get(chem['cid'], set())

        def good_name_criteria(name):
            if not name:
                return False

            norm_name = self._normalize_chem_name(name, is_clean=True)
            if norm_name in filtered_synonyms:
                return False

            DISCARD_PATTERNS = [
                r'[:=%<>@/\\_.#&*";]',
                r'[-,.]$',
                self.CAS_PATTERN,
                r'unii-',
                r'\(\d:\d\)',
                r'\d{3,}',
                r'[a-z]{14}-[a-z]{10}-[a-z]',
                r'\s{2,}',
                r'\b[nm]m\b',
                r'-,'
            ]

            DISCARD_PATTERNS_IGNORECASE = [
                r'\bpowder\b',
                r'\bbeads\b',
                r'\bimpurity\b',
                r'\bgrade\b',
                r'\bdry\b'
            ]

            name = name.lower()
            if any(re.search(p, name) for p in DISCARD_PATTERNS) or any(
                    re.search(p, name, flags=re.IGNORECASE) for p in DISCARD_PATTERNS_IGNORECASE):
                return False
            # UNII identifiers
            if re.fullmatch(r'[a-z0-9]{10}', name) and re.search(r'[abdefgijklmqrtuvwxyz]', name) and re.search(r'\d', name):
                return False

            return True

        def clean_synonym(name):

            name = name.strip()

            tags = [
                'USP', 'EP', 'HSDB', 'EP MONOGRAPH', 'WHO-DD', 'USAN', 'MI',
                'USP MONOGRAPH', 'VANDF', 'CZECH', 'ACGIH', 'NDIPA', 'German',
                'VAN', 'INN', '9CI', '8CI', 'JAN', 'French', 'Latin', 'ISO', 'Standard',
                'natural', 'HPUS', 'IARC', 'NF', 'INCI', 'TN', 'Dutch', 'Italian',
                'FCC', 'DCIT', 'BAN', 'ORANGE BOOK', 'Spanish', 'IUPAC', 'FHFI', 'Polish'
            ]

            pattern = r'[\(\[]\s*(?:' + '|'.join(re.escape(tag) for tag in tags) + r')\s*[\)\]]$'
            name = re.sub(pattern, '', name, flags=re.IGNORECASE).strip()

            name = re.sub(r'\bfume\b', '', name, flags=re.IGNORECASE).strip()

            return name

        synonyms = chem['cmpdsynonym']
        if not synonyms:
            return False

        if not isinstance(synonyms, list):
            synonyms = [synonyms]

        synonyms = map(lambda x: clean_synonym(x), synonyms)
        synonyms = list(filter(lambda x: good_name_criteria(x), synonyms))
        if not synonyms:
            return False

        if good_name_criteria(chem['iupacname']):
            synonyms = [chem['iupacname']] + synonyms

        if chem['cmpdname']:
            chem['cmpdname'] = clean_synonym(chem['cmpdname'])

        if chem['cmpdname'] is None or not good_name_criteria(chem['cmpdname']):
            chem['cmpdname'] = synonyms[0].lower()
        else:
            synonyms = [chem['cmpdname']] + synonyms

        synonyms = list(dict.fromkeys(synonyms))
        chem['cmpdsynonym'] = synonyms[:self.max_synonyms_thr]

        return True
    


    def __process_chem_single(self, chem, unique_inchikeys_chems, force):
        cid = chem['cid']

        try:
            if chem['charge'] != 0:
                return False

            if cid in self.cids_blacklist:
                return False

            if chem['complexity'] > self.complexity_thr:
                return False

            if '/i' in chem['inchi']:
                return False

            try:
                mol = Chem.MolFromSmiles(chem['smiles'])
                chem['smiles'] = Chem.MolToSmiles(mol, canonical=True)
            except Exception:
                return False

            # Don't use 'force' since we remove CAS numbers from synonyms anyway
            if 'cas' not in chem:
                chem['cas'] = self._extract_chem_cas_numbers(chem)

            if 'iupacname' not in chem:
                chem['iupacname'] = None

            if not self.__process_chem_synonyms(chem):
                return False

            def is_hydrate_inchi(inchi: str) -> bool:
                try:
                    formula_part = inchi.split("/")[1]
                except IndexError:
                    # No layers, take everything after "InChI=1S/"
                    formula_part = inchi[8:]

                # Split by dot (separate molecules)
                molecules = formula_part.split(".")

                # Check for water molecules
                water_pattern = re.compile(r'^(\d*)H2O$', re.IGNORECASE)
                for mol in molecules:
                    if water_pattern.match(mol):
                        return True

                return False

            # Filter hydrates (two checks for reliability)
            if is_hydrate_inchi(chem['inchi']) and 'hydrate' in chem['cmpdname']:
                return False

            if 'ECFP4_fp' not in chem or force:
                chem['ECFP4_fp'] = self._get_mol_fingerprint(mol)

            if 'bertz_complexity' not in chem or force:
                chem['bertz_complexity'] = self._get_mol_bertz_complexity(mol)

            if 'organic' not in chem or force:
                chem['organic'] = self._get_mol_organic_mark(mol)

            if 'wiki' not in chem or force:
                chem['wiki'] = self.cid_wiki_map.get(cid, None)

            if 'inchi_snone' not in chem or force:
                chem['inchi_snone'] = inchi.MolToInchi(mol, options="/SNon")

            if 'inchikey_snone' not in chem or force:
                chem['inchikey_snone'] = inchi.MolToInchiKey(mol, options="/SNon")

            inchikey = chem['inchikey_snone']
            if inchikey in unique_inchikeys_chems:
                old_chem = unique_inchikeys_chems[inchikey]
                old_inchi = old_chem['inchi']
                curr_inchi = chem['inchi']
                if len(curr_inchi) < len(old_inchi):
                    unique_inchikeys_chems[inchikey] = chem
            else:
                unique_inchikeys_chems[inchikey] = chem

            return True

        except Exception as e:
            self.log(f"Exception during processing compound with CID {cid}: {e}")
            return False

    def _process_chems(self, force=False):
        __unique_inchikeys_chems = dict()
        for chem in self.chems:
            self.__process_chem_single(chem, __unique_inchikeys_chems, force=force)

        return __unique_inchikeys_chems

    def organize_chems_file(self, force=False):
        shutil.copy(self.chems_fn, f"{self.chems_fn}.backup")

        initial_chems_num = len(self.chems)

        unique_chems = list(self._process_chems(force=force).values())
        self._update_chems(unique_chems)

        print(f"Discarded {initial_chems_num - len(unique_chems)} chems")

    def __chem_name_to_ascii(self, chem_name_raw):
        unicode_map = {
            '‑': '-',
            'α': 'alpha',
            'γ': 'gamma,',
            '–': '-',
            '\u2019': "'"
        }
        chem_name_ascii = ""
        for char in chem_name_raw:
            if not char.isascii():
                if char in unicode_map:
                    char = unicode_map[char]
                else:
                    char = ""
            chem_name_ascii += char

        return chem_name_ascii

    def _clean_chem_name(self, chem_name_raw, is_clean=False):
        chem_name = self.__chem_name_to_ascii(chem_name_raw)
        chem_name = chem_name.strip()
        chem_name = re.sub(r'\s+', ' ', chem_name)

        if not is_clean:
            chem_name = chem_name.strip('`\'".,;:')
            chem_name = re.sub(r'^\d+ ', '', chem_name)

        return chem_name

    def _normalize_chem_name(self, chem_name_raw, is_clean=False):
        chem_name = self._clean_chem_name(chem_name_raw, is_clean=is_clean)
        chem_name = chem_name.lower()
        chem_name = chem_name.strip()
        chem_name = chem_name.replace("aluminum", "aluminium")

        if not is_clean:
            chem_name = re.sub(r' \([^\d]+\)$', '', chem_name)
            chem_name = chem_name.replace(' vapor', '')
            chem_name = chem_name.replace(' dust', '')
            chem_name = chem_name.replace('solution', '')
            chem_name = chem_name.replace('concentrated', '')
            chem_name = chem_name.replace('dilute ', '')
            chem_name = chem_name.replace('fuming ', '')
            chem_name = chem_name.replace('solid', '')
            chem_name = chem_name.replace('glacial ', '')
            chem_name = chem_name.replace('elemental', '')
            chem_name = chem_name.replace(' metal', '')
            chem_name = chem_name.replace('aqueous', '')
            chem_name = chem_name.replace(' gas', '')
            chem_name = chem_name.replace('hot ', '')
            chem_name = chem_name.replace('uv light', 'light')
            chem_name = chem_name.replace('blue light', 'light')
            chem_name = chem_name.replace('ultraviolet light', 'light')

            if "catalyst" in chem_name or 'raney nickel' in chem_name:
                chem_name = "catalyst"

        chem_name = re.sub(r'\s+', '', chem_name)

        return chem_name


if __name__ == "__main__":
    parse = ChemsParsePubchem('data/')
    parse.organize_chems_file()