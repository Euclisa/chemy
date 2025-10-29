import json
import re
import glob
import os
import shutil

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GraphDescriptors, inchi, rdMolDescriptors, Descriptors

from functools import cached_property
from types import MappingProxyType
from itertools import chain

import numpy as np
import matplotlib.pyplot as plt

import inspect

from rich.progress import track, Progress

from chems_db import ChemsDB


# Disable all RDKit warnings and info messages
RDLogger.DisableLog('rdApp.*')


class ChemsProperties(ChemsDB):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.chems_dir = os.path.join(self.data_dir, 'chems')
        self.chems_mapped_fn = os.path.join(self.chems_dir, 'mapped')
        self.chems_unmapped_fn = os.path.join(self.chems_dir, 'unmapped')

        os.makedirs(self.chems_mapped_fn, exist_ok=True)
        os.makedirs(self.chems_unmapped_fn, exist_ok=True)

        self.chems_wiki_fn = os.path.join(self.data_dir, 'chems', "chems_wiki.jsonl")
        self.chems_edges_fn = os.path.join(self.data_dir, 'chems', 'chems_edges.jsonl')
        self.elements_fn = os.path.join(self.data_dir, 'chems', 'elements.jsonl')

        self.categories_fn = os.path.join(self.data_dir, 'misc', "categories.jsonl")
        self.background_cids_fn = os.path.join(self.data_dir, 'misc', 'background_cids.json')
        self.cids_blacklist_fn = os.path.join(self.data_dir, 'misc', 'cids_blacklist.jsonl')
        self.cids_filtered_synonyms_fn = os.path.join(self.data_dir, 'misc', 'cids_filtered_synonyms.jsonl')

        self.chems_properties_dir = os.path.join(self.data_dir, 'chems_properties')
        self.chems_properties_assets_dir = os.path.join(self.data_dir, 'assets', 'chems_properties')

        self.chems_categories_fn = os.path.join(self.chems_properties_dir, "chems_categories.jsonl")

        self._file_sorting_prefs[self.chems_mapped_fn] = 'complexity'
        self._file_sorting_prefs[self.chems_unmapped_fn] = 'cid'
        
        self._file_sorting_prefs[self.chems_categories_fn] = 'cid'
        self._file_sorting_prefs[self.chems_wiki_fn] = 'cid'
        self._file_sorting_prefs[self.chems_edges_fn] = 'eid'
        self._file_sorting_prefs[self.elements_fn] = None

        self._file_sorting_prefs[self.categories_fn] = 'code'
        self._file_sorting_prefs[self.cids_blacklist_fn] = 'cid'
        self._file_sorting_prefs[self.cids_filtered_synonyms_fn] = 'cid'

        self._dir_vault_prefs[self.chems_mapped_fn] = 'chems_'
        self._dir_vault_prefs[self.chems_unmapped_fn] = 'chems_'

        self.complexity_thr = 700
        self.max_synonyms_thr = 150

        self.unmapped_bertz_complexity_thr = 1000
        self.unmapped_heavy_count_thr = 15

        self.unknown_name_ph = "<Unknown>"

        self.CAS_PATTERN = r'\d{2,7}-\d{2}-\d'
    

    @cached_property
    def chems(self):
        return tuple(chain(self.chems_mapped, self.chems_unmapped))
    
    @cached_property
    def chems_mapped(self):
        mapped = self._load_jsonl(self.chems_mapped_fn)
        return tuple(mapped)
    
    @cached_property
    def chems_unmapped(self):
        unmapped = self._load_jsonl(self.chems_unmapped_fn)
        return tuple(unmapped)
    
    @cached_property
    def cids_mapped(self):
        return frozenset(chem['cid'] for chem in self.chems_mapped)
    
    @cached_property
    def cids_unmapped(self):
        return frozenset(chem['cid'] for chem in self.chems_unmapped)
    
    @cached_property
    def cids_blacklist(self):
        return frozenset([x['cid'] for x in self._load_jsonl(self.cids_blacklist_fn)])

    @cached_property
    def cid_chem_map(self):
        return MappingProxyType({chem['cid']: chem for chem in self.chems})
    
    @cached_property
    def mapped_cid_chem_map(self):
        return MappingProxyType({chem['cid']: chem for chem in self.chems_mapped})
    
    @cached_property
    def inchikey_cid_map(self):
        return MappingProxyType({chem['inchikey_snone']: chem['cid'] for chem in self.chems})
    
    @cached_property
    def mapped_inchikey_cid_map(self):
        return MappingProxyType({chem['inchikey_snone']: chem['cid'] for chem in self.chems_mapped})
    
    def __get_norm_inchi_cid_map_chems(self, chems):
        _norm_inchi_cid_map = dict()
        for chem in chems:
            norm_inchi = self._strip_inchi_layers(chem['inchi_snone'])
            _norm_inchi_cid_map[chem['cid']] = norm_inchi
        return _norm_inchi_cid_map
    
    @cached_property
    def norm_inchi_cid_map(self):
        return MappingProxyType(self.__get_norm_inchi_cid_map_chems(self.chems))
    
    @cached_property
    def mapped_norm_inchi_cid_map(self):
        return MappingProxyType(self.__get_norm_inchi_cid_map_chems(self.chems_mapped))
    
    @cached_property
    def smiles_cid_map(self):
        return MappingProxyType({chem['smiles']: chem['cid'] for chem in self.chems})
    
    @cached_property
    def mapped_smiles_cid_map(self):
        return MappingProxyType({chem['smiles']: chem['cid'] for chem in self.chems_mapped})

    @cached_property
    def cid_mf_map(self):
        return MappingProxyType({chem['cid']: chem['mf'] for chem in self.chems})
    

    @cached_property
    def cid_wiki_map(self):
        wiki_entries = self._load_jsonl(self.chems_wiki_fn)
        return MappingProxyType({entry['cid']: entry['wiki'] for entry in wiki_entries})


    @cached_property
    def cids_filtered_synonyms(self):
        entries = self._load_jsonl(self.cids_filtered_synonyms_fn)
        return MappingProxyType({x['cid']: set(x['synonyms']) for x in entries})

    @cached_property
    def name_cid_map(self):
        _name_cid_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            _name_cid_map[self._normalize_chem_name(chem['cmpdname'], is_clean=True)] = cid
            for syn in chem['cmpdsynonym']:
                _name_cid_map[self._normalize_chem_name(syn, is_clean=True)] = cid

        return MappingProxyType(_name_cid_map)

    @cached_property
    def cas_cid_map(self):
        _cas_cid_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            cas_list = chem['cas']
            if cas_list:
                for cas in cas_list:
                    _cas_cid_map[cas] = cid
        
        return MappingProxyType(_cas_cid_map)
    
    @cached_property
    def symb_to_el(self):
        return MappingProxyType({el['symbol']: el for el in self._load_jsonl(self.elements_fn)})
    

    def _clear_cached_property(self, name: str):
        cls_attr = getattr(type(self), name, None)
        if isinstance(cls_attr, cached_property):
            if name in self.__dict__:
                del self.__dict__[name]
                return True
            return False
        raise AttributeError(f"'{name}' is not a cached_property of {type(self).__name__}")


    def __clear_runtime_chems_properties(self):
        for name, attr in inspect.getmembers(type(self)):
            if isinstance(attr, cached_property) and name in self.__dict__:
                delattr(self, name)


    def _update_chems(self, new_chems):

        mapped_chems = [chem for chem in new_chems if chem.get('cid', 0) > 0]
        self._write_jsonl(mapped_chems, self.chems_mapped_fn)

        unmapped_chems = [chem for chem in new_chems if chem.get('cid', 0) < 0]
        self._write_jsonl(unmapped_chems, self.chems_unmapped_fn)

        self.__clear_runtime_chems_properties()
    

    def _update_unmapped_chems(self, new_chems):

        unmapped_chems = [chem for chem in new_chems if chem.get('cid', 0) < 0]
        self._write_jsonl(unmapped_chems, self.chems_unmapped_fn)

        self.__clear_runtime_chems_properties()


    def _update_cids_blacklist(self, cids):
        cids = set(filter(lambda cid: cid > 0, cids))
        unique_new_cids = cids - self.cids_blacklist
        if not unique_new_cids:
            return

        with open(self.cids_blacklist_fn, 'a') as f:
            for cid in unique_new_cids:
                name = self.cid_chem_map.get(cid, {}).get('cmpdname', 'Unknown')
                f.write(json.dumps({'cid': cid, 'name': name}) + '\n')
            f.flush()

        self._clear_cached_property('cids_blacklist')
    

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
                if not auto_refresh:
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
    

    def _good_name_criteria(self, name: str):
        if not name:
            return False
        
        if name == self.unknown_name_ph:
            return True
        
        name = name.strip()
        
        DISCARD_WORDS_PART = [
            'oil', 'solid', 'liquid', 'dry', 'powder', 'nanopowder',
            'beads', 'impurity', 'grade', 'intermediate', 'title', 'desired',
            'material', 'solution', 'syrup', 'crystals', 'residue', 'compound',
            'product', 'titled'
        ]

        discard_word_part_pattern = '(' + '|'.join([r"\b"+word+r"\b" for word in DISCARD_WORDS_PART]) + ')'

        DISCARD_WORDS_WHOLE = [
            'acetate', 'acid', 'salt', 'phosphonate', 'ester'
        ]

        discard_word_whole_pattern = '(' + '|'.join([r"^"+word+r"$" for word in DISCARD_WORDS_WHOLE]) + ')'

        DISCARD_PATTERNS = [
            r'[:=%<>@/\\_.#&*";?!]',
            r'[-,.]$',
            self.CAS_PATTERN,
            r'unii-',
            r'\(\d:\d\)',
            r'\d{3,}',
            r'[a-z]{14}-[a-z]{10}-[a-z]',
            r'\s{2,}',
            r'\b[nm]m\b',
            r'-,',
            r'^\d+$',
            r'^[a-z]?[0-9-()\s]+[a-z]?$',
            discard_word_part_pattern,
            discard_word_whole_pattern
        ]

        name = name.lower()
        if any(re.search(p, name) for p in DISCARD_PATTERNS):
            return False
        # UNII identifiers
        if re.fullmatch(r'[a-z0-9]{10}', name) and re.search(r'[abdefgijklmqrtuvwxyz]', name) and re.search(r'\d', name):
            return False

        return True

    def __process_chem_synonyms(self, chem):

        def good_name_criteria(name):
            if not name:
                return False

            filtered_synonyms = self.cids_filtered_synonyms.get(chem['cid'], set())
            norm_name = self._normalize_chem_name(name, is_clean=True)
            if norm_name in filtered_synonyms:
                return False
            
            return self._good_name_criteria(name)

        def clean_synonym(name):

            if name == self.unknown_name_ph:
                return name

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


    def _merge_synonyms(self, chem1, chem2):
        result_synonyms = dict()
        synonyms1 = chem1['cmpdsynonym']
        synonyms2 = chem2['cmpdsynonym']
    
        if len(synonyms1) > len(synonyms2):
            synonyms1, synonyms2 = synonyms2, synonyms1

        for i in range(len(synonyms1)):
            result_synonyms[synonyms1[i]] = 0
            result_synonyms[synonyms2[i]] = 0
        
        for i in range(len(synonyms1), len(synonyms2)):
            result_synonyms[synonyms2[i]] = 0
        
        return list(dict.fromkeys(result_synonyms))



    def _process_chem_single(self, chem, force=False):
        cid = chem['cid']

        try:

            if chem['charge'] != 0:
                return None

            if cid in self.cids_blacklist:
                return None

            if cid > 0 and chem['complexity'] > self.complexity_thr:
                return None

            if '/i' in chem['inchi']:
                return None

            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                return None
            chem['smiles'] = Chem.MolToSmiles(mol, canonical=True)

            # Don't use 'force' since we remove CAS numbers from synonyms anyway
            if 'cas' not in chem:
                chem['cas'] = self._extract_chem_cas_numbers(chem)

            if 'iupacname' not in chem:
                chem['iupacname'] = None

            if cid > 0:
                if not self.__process_chem_synonyms(chem):
                    return None
            else:
                if not self._good_name_criteria(chem['cmpdname']):
                    chem['cmpdname'] = self.unknown_name_ph

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
                for mol_str in molecules:
                    if water_pattern.match(mol_str):
                        return True

                return None

            if cid > 0:
                # Filter hydrates (two checks for reliability)
                if is_hydrate_inchi(chem['inchi']) and 'hydrate' in chem['cmpdname']:
                    return None

            if 'ECFP4_fp' not in chem or force:
                chem['ECFP4_fp'] = self._get_mol_fingerprint(mol)

            if 'bertz_complexity' not in chem or force:
                chem['bertz_complexity'] = self._get_mol_bertz_complexity(mol)
            
            if 'heavy_count' not in chem or force:
                chem['heavy_count'] = mol.GetNumHeavyAtoms()
            
            if cid < 0 and not self._unmapped_complexity_criteria(chem):
                return None

            if 'organic' not in chem or force:
                chem['organic'] = self._get_mol_organic_mark(mol)

            if 'wiki' not in chem or force:
                chem['wiki'] = self.cid_wiki_map.get(cid, None)

            if 'inchi_snone' not in chem or force:
                chem['inchi_snone'] = inchi.MolToInchi(mol, options="/SNon")

            if 'inchikey_snone' not in chem or force:
                chem['inchikey_snone'] = inchi.MolToInchiKey(mol, options="/SNon")

            if not chem['inchi_snone'] or not chem['inchikey_snone']:
                return None

            return chem

        except Exception as e:
            self.log_warn(f"Exception during processing compound with CID {cid}: {e}")
            return None

    def _process_chems(self, force=False):
        processed_chems = []
        for chem in self._rich_track(self.chems, "Processing compounds"):
            if (chem := self._process_chem_single(chem, force=force)):
                processed_chems.append(chem)

        return processed_chems
    

    def _get_unique_inchikeys_chems(self, organize=False):
        if organize:
            return self._process_chems()
        else:
            unique_inchikeys_chems = dict()
            for chem in self.chems:
                inchikey = chem['inchikey_snone']
                if not inchikey:
                    raise Exception(f"Invalid InChI-key found in main compounds file. Run with 'organize=True'")
                if inchikey in unique_inchikeys_chems:
                    raise Exception(f"Duplicates found in main compounds file. Run with 'organize=True'")
                unique_inchikeys_chems[inchikey] = chem

            return unique_inchikeys_chems



    def organize_chems_file(self, force=False):
        initial_chems_num = len(self.chems)

        processed_chems = self._process_chems(force=force)
        self._update_chems(processed_chems)

        print(f"Discarded {initial_chems_num - len(processed_chems)} chems")

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
        if chem_name_raw == self.unknown_name_ph:
            return chem_name_raw

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
    

    def __build_basic_chem(self, cid, name, smiles):
        if not name:
            name = "<Unknown>"
        
        if cid >= 0:
            raise ValueError(f"Invalid CID: {cid}")

        chem = dict()
        chem['cid'] = cid
        chem['cmpdname'] = name

        try:
            mol = Chem.MolFromSmiles(smiles)
            inchi = Chem.MolToInchi(mol)
            if not inchi:
                raise Exception("Failed to obtain inchi")
            inchikey = Chem.InchiToInchiKey(inchi)
            if not inchikey:
                raise Exception("Failed to obtain inchikey")
            mf = rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
            charge = Chem.GetFormalCharge(mol)
            
        except Exception as e:
            self.log_warn(f"Failed process smiles '{smiles}': {e}")
            return None
        
        chem['smiles'] = smiles
        chem['inchi'] = inchi
        chem['inchikey'] = inchikey
        chem['mf'] = mf
        chem['mw'] = mw
        chem['complexity'] = None
        chem['charge'] = charge
        chem['cmpdsynonym'] = [name]

        return chem
    

    def _strip_inchi_layers(self, inchi: str) -> str:
        if not inchi.startswith("InChI="):
            raise ValueError("Input does not appear to be a valid InChI string")

        cleaned = re.sub(r"/[pqi][^/]*", "", inchi)
        return cleaned


    def _get_mol_norm_inchi(self, mol):
        try:
            return self._strip_inchi_layers(inchi.MolToInchi(mol, options="/SNon"))
        except Exception:
            return None

    def _get_chem_norm_inchi(self, chem):
        return self._strip_inchi_layers(chem['inchi_snone'])


    def _add_new_chems(self, name_smiles):
        res_chems = list(self.chems)
        init_chems_num = len(res_chems)
        
        occupied_cids = [chem['cid'] for chem in self.cids_unmapped]
        free_cid = min(occupied_cids)-1 if occupied_cids else -1

        norm_inchi_set = set(self.norm_inchi_cid_map.keys())

        for name, smiles in self._rich_track(name_smiles, "Adding compounds"):
            chem = self.__build_basic_chem(free_cid, name, smiles)
            if chem:
                status = self._process_chem_single(chem)
                if status:
                    norm_inchi = self._get_chem_norm_inchi(chem)
                    if norm_inchi not in norm_inchi_set:
                        res_chems.append(chem)
                        norm_inchi_set.add(norm_inchi)
                        free_cid -= 1

        self._update_chems(res_chems)

        result_chems_num = len(res_chems)
        self.log(f"Processed {len(name_smiles)} entries; added {result_chems_num-init_chems_num} new compounds")
    

    def _get_fp_tanimoto(self, fp1, fp2):
        and_pop = sum((ai & bi).bit_count() for ai, bi in zip(fp1["bits"], fp2["bits"]))
        or_pop = fp1["popcount"] + fp2["popcount"] - and_pop

        return 1.0 if or_pop == 0 else and_pop / or_pop
        

    def _unmapped_complexity_criteria(self, chem):
        if chem['cid'] > 0:
            return True

        return chem['bertz_complexity'] <= self.unmapped_bertz_complexity_thr or chem['heavy_count'] <= self.unmapped_heavy_count_thr


    def _unmapped_complexity_criteria_mol(self, mol):
        try:
            return self._get_mol_bertz_complexity(mol) <= self.unmapped_bertz_complexity_thr or mol.GetNumHeavyAtoms() <= self.unmapped_heavy_count_thr
        except Exception:
            return False



if __name__ == "__main__":
    parse = ChemsProperties('data/')
    parse.plot_fp_popcount_bertz_graph()