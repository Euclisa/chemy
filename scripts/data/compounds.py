import json
import re
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, GraphDescriptors, inchi, rdMolDescriptors, Descriptors

from functools import cached_property
from types import MappingProxyType
from itertools import chain

import inspect
from collections import Counter

from rich.table import Table
from rich.rule import Rule

from . import chem_names


# Disable all RDKit warnings and info messages
RDLogger.DisableLog('rdApp.*')


class CompoundStore:

    def __init__(self, data_dir, store, logger):
        self.data_dir = data_dir
        self.store = store
        self.logger = logger

        self.chems_dir = os.path.join(self.data_dir, 'chems')
        self.chems_mapped_fn = os.path.join(self.chems_dir, 'mapped')
        self.chems_unmapped_fn = os.path.join(self.chems_dir, 'unmapped')

        os.makedirs(self.chems_mapped_fn, exist_ok=True)
        os.makedirs(self.chems_unmapped_fn, exist_ok=True)

        self.chems_edges_fn = os.path.join(self.data_dir, 'chems', 'chems_edges.jsonl')
        self.elements_fn = os.path.join(self.data_dir, 'chems', 'elements.jsonl')
        self.cids_blacklist_fn = os.path.join(self.data_dir, 'chems', 'cids_blacklist.jsonl')
        self.chem_names_blacklisted_fn = os.path.join(self.data_dir, 'unmapped_names_blacklisted.txt')
        self.chem_smiles_blacklisted_fn = os.path.join(self.data_dir, 'unmapped_smiles_blacklisted.txt')

        self.categories_fn = os.path.join(self.data_dir, 'misc', "categories.jsonl")
        self.background_cids_fn = os.path.join(self.data_dir, 'misc', 'background_cids.json')
        self._legacy_cids_filtered_synonyms_fn = os.path.join(self.data_dir, 'misc', 'cids_filtered_synonyms.jsonl')
        self.synonym_filter_decisions_fn = os.path.join(self.data_dir, 'misc', 'synonym_filter_decisions.jsonl')

        self.chems_properties_dir = os.path.join(self.data_dir, 'chems_properties')
        self.chems_properties_assets_dir = os.path.join(self.data_dir, 'assets', 'chems_properties')

        self.chems_categories_fn = os.path.join(self.chems_properties_dir, "chems_categories.jsonl")
        self.chems_wiki_fn = os.path.join(self.chems_properties_dir, 'wiki', "chems_wiki.jsonl")
        self.chems_hazards_wiki_fn = os.path.join(self.chems_properties_dir, 'wiki', "chems_hazards_wiki.jsonl")

        self.structures_dir = os.path.join(self.data_dir, 'structures')

        store.register_sorting(self.chems_mapped_fn, 'bertz_complexity')
        store.register_sorting(self.chems_unmapped_fn, 'cid')
        store.register_sorting(self.chems_categories_fn, 'cid')
        store.register_sorting(self.chems_wiki_fn, 'cid')
        store.register_sorting(self.chems_edges_fn, 'eid')
        store.register_sorting(self.elements_fn, None)
        store.register_sorting(self.categories_fn, 'code')
        store.register_sorting(self.cids_blacklist_fn, 'cid')
        store.register_sorting(self.synonym_filter_decisions_fn, 'cid')

        store.register_vault(self.chems_mapped_fn, 'chems_')
        store.register_vault(self.chems_unmapped_fn, 'chems_')

        self.complexity_thr = 700
        self.max_synonyms_thr = 150

        self.unmapped_bertz_complexity_thr = 1000
        self.unmapped_heavy_count_thr = 15

        self.unknown_name_ph = "<Unknown>"
        self.null_cid = 0

        self.CAS_PATTERN = chem_names.CAS_PATTERN

    # --- Cached properties ---

    @cached_property
    def chems(self):
        return tuple(chain(self.chems_mapped, self.chems_unmapped))

    @cached_property
    def chems_mapped(self):
        mapped = self.store.load_jsonl(self.chems_mapped_fn)
        return tuple(mapped)

    @cached_property
    def chems_unmapped(self):
        unmapped = self.store.load_jsonl(self.chems_unmapped_fn)
        return tuple(unmapped)

    @cached_property
    def cids_mapped(self):
        return frozenset(chem['cid'] for chem in self.chems_mapped)

    @cached_property
    def cids_unmapped(self):
        return frozenset(chem['cid'] for chem in self.chems_unmapped)

    @cached_property
    def cids_blacklist(self):
        return frozenset([x['cid'] for x in self.store.load_jsonl(self.cids_blacklist_fn)])

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
            norm_inchi = self.strip_inchi_layers(chem['inchi_snone'])
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
        wiki_entries = self.store.load_jsonl(self.chems_wiki_fn)
        return MappingProxyType({entry['cid']: entry['wiki'] for entry in wiki_entries})

    @cached_property
    def synonym_filter_decisions(self):
        entries = self.store.load_jsonl(self.synonym_filter_decisions_fn)
        if not entries:
            entries = self._load_legacy_synonym_filter_decisions()
        return tuple(self._normalize_synonym_filter_decision(entry) for entry in entries)

    @cached_property
    def filtered_synonyms_by_cid(self):
        filtered = dict()
        for entry in self.synonym_filter_decisions:
            if entry['decision'] != 'remove_synonym':
                continue
            filtered.setdefault(entry['cid'], set()).add(entry['norm_synonym'])
        return MappingProxyType(filtered)

    @cached_property
    def name_cids_map(self):
        _name_cids_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            norm_name = self.normalize_chem_name(chem['cmpdname'])
            _name_cids_map.setdefault(norm_name, set()).add(cid)
            for syn in chem['cmpdsynonym']:
                norm_name = self.normalize_chem_name(syn)
                _name_cids_map.setdefault(norm_name, set()).add(cid)

        return MappingProxyType({
            name: tuple(sorted(cids))
            for name, cids in _name_cids_map.items()
        })

    @cached_property
    def unique_name_cid_map(self):
        return MappingProxyType({
            name: positive_cids[0]
            for name, cids in self.name_cids_map.items()
            if len(positive_cids := tuple(cid for cid in cids if cid > 0)) == 1
        })

    @cached_property
    def ambiguous_name_cids_map(self):
        return MappingProxyType({
            name: cids
            for name, cids in self.name_cids_map.items()
            if len([cid for cid in cids if cid > 0]) > 1
        })

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
        return MappingProxyType({el['symbol']: el for el in self.store.load_jsonl(self.elements_fn)})

    # --- Synonym filter decision persistence ---

    def _load_legacy_synonym_filter_decisions(self):
        legacy_entries = self.store.load_jsonl(self._legacy_cids_filtered_synonyms_fn)
        decisions = []
        for entry in legacy_entries:
            cid = entry['cid']
            for norm_synonym in entry.get('synonyms', []):
                decisions.append({
                    'cid': cid,
                    'norm_synonym': norm_synonym,
                    'raw_synonyms': [norm_synonym],
                    'decision': 'remove_synonym',
                    'reason': 'manual_conflict',
                    'source': 'legacy_cids_filtered_synonyms',
                    'other_cids': [],
                })
        return decisions

    def _normalize_synonym_filter_decision(self, entry):
        norm_synonym = entry.get('norm_synonym')
        if norm_synonym is None and entry.get('synonym') is not None:
            norm_synonym = self.normalize_chem_name(entry['synonym'])

        return {
            'cid': int(entry['cid']),
            'norm_synonym': norm_synonym,
            'raw_synonyms': sorted(set(entry.get('raw_synonyms') or [])),
            'decision': entry.get('decision', 'remove_synonym'),
            'reason': entry.get('reason', 'manual_conflict'),
            'source': entry.get('source', 'manual_conflict_resolver'),
            'other_cids': sorted(set(entry.get('other_cids') or [])),
        }

    def write_synonym_filter_decisions(self, decisions):
        merged = dict()
        for entry in decisions:
            entry = self._normalize_synonym_filter_decision(entry)
            if entry['cid'] <= 0 or not entry['norm_synonym']:
                continue
            key = (
                entry['cid'],
                entry['norm_synonym'],
                entry['decision'],
                entry['reason'],
                entry['source'],
            )
            if key not in merged:
                merged[key] = entry
                continue
            merged[key]['raw_synonyms'] = sorted(
                set(merged[key]['raw_synonyms']) | set(entry['raw_synonyms'])
            )
            merged[key]['other_cids'] = sorted(
                set(merged[key]['other_cids']) | set(entry['other_cids'])
            )

        self.store.write_jsonl(list(merged.values()), self.synonym_filter_decisions_fn)
        self.clear_cached_property('synonym_filter_decisions')
        self.clear_cached_property('filtered_synonyms_by_cid')

    def migrate_legacy_synonym_filter_decisions(self):
        current = list(self.store.load_jsonl(self.synonym_filter_decisions_fn))
        migrated = self._load_legacy_synonym_filter_decisions()
        self.write_synonym_filter_decisions(current + migrated)
        return len(migrated)

    def audit_synonym_filtering(self, chems=None, rejected_by_reason=None):
        chems = self.chems_mapped if chems is None else tuple(chems)
        rejected_by_reason = Counter(rejected_by_reason or {})

        accepted_synonyms = sum(len(chem.get('cmpdsynonym') or []) for chem in chems)
        normalized_to_cids = dict()
        exported_synonyms = 0
        for chem in chems:
            cid = chem['cid']
            filtered_synonyms = self.filtered_synonyms_by_cid.get(cid, set())
            for syn in chem.get('cmpdsynonym') or []:
                result = chem_names.classify_name(
                    syn,
                    manual_filtered=filtered_synonyms,
                    unknown_name=self.unknown_name_ph,
                )
                if not result.accepted and result.reason:
                    rejected_by_reason[result.reason] += 1
                norm_name = self.normalize_chem_name(syn)
                normalized_to_cids.setdefault(norm_name, set()).add(cid)

        for chem in chems:
            cid = chem['cid']
            for syn in chem.get('cmpdsynonym') or []:
                norm_name = self.normalize_chem_name(syn)
                if normalized_to_cids.get(norm_name) == {cid}:
                    exported_synonyms += 1

        return {
            'compounds': len(chems),
            'accepted_synonyms': accepted_synonyms,
            'rejected_by_reason': dict(sorted(rejected_by_reason.items())),
            'unique_names': sum(1 for cids in normalized_to_cids.values() if len(cids) == 1),
            'ambiguous_names': sum(1 for cids in normalized_to_cids.values() if len(cids) > 1),
            'exportable_unique_synonyms': exported_synonyms,
            'filter_decisions': len(self.synonym_filter_decisions),
            'filter_decisions_applied': sum(len(syns) for syns in self.filtered_synonyms_by_cid.values()),
        }

    def regenerate_mapped_synonyms(self):
        migrated_decisions = self.migrate_legacy_synonym_filter_decisions()
        rejected_by_reason = Counter()
        processed_chems = []
        initial_chems_num = len(self.chems_mapped)

        for chem in self.logger.track(self.chems_mapped, "Regenerating mapped synonyms"):
            filtered_synonyms = self.filtered_synonyms_by_cid.get(chem['cid'], set())
            names_to_check = list(chem.get('cmpdsynonym') or [])
            names_to_check.extend([chem.get('cmpdname'), chem.get('iupacname')])
            for name in names_to_check:
                result = chem_names.classify_name(
                    name,
                    manual_filtered=filtered_synonyms,
                    unknown_name=self.unknown_name_ph,
                )
                if not result.accepted and result.reason:
                    rejected_by_reason[result.reason] += 1

            processed_chem = self.process_chem_single(dict(chem), force=set())
            if processed_chem:
                processed_chems.append(processed_chem)

        self.update_mapped_chems(processed_chems)
        audit = self.audit_synonym_filtering(rejected_by_reason=rejected_by_reason)
        audit['legacy_decisions_migrated'] = migrated_decisions
        audit['discarded_compounds'] = initial_chems_num - len(processed_chems)
        return audit

    # --- Cache invalidation ---

    def clear_cached_property(self, name):
        cls_attr = getattr(type(self), name, None)
        if isinstance(cls_attr, cached_property):
            if name in self.__dict__:
                del self.__dict__[name]
                return True
            return False
        raise AttributeError(f"'{name}' is not a cached_property of {type(self).__name__}")

    def _clear_all_cached(self):
        for name, attr in inspect.getmembers(type(self)):
            if isinstance(attr, cached_property) and name in self.__dict__:
                delattr(self, name)

    # --- Mapped/unmapped helpers ---

    def is_chem_mapped(self, chem):
        return chem['cid'] > 0

    def _get_mapped_chems_list(self, chems):
        return [chem for chem in chems if self.is_chem_mapped(chem)]

    def _get_unmapped_chems_list(self, chems):
        return [chem for chem in chems if not self.is_chem_mapped(chem)]

    def update_chems(self, new_chems):
        mapped_chems = self._get_mapped_chems_list(new_chems)
        self.store.write_jsonl(mapped_chems, self.chems_mapped_fn)

        unmapped_chems = self._get_unmapped_chems_list(new_chems)
        self.store.write_jsonl(unmapped_chems, self.chems_unmapped_fn)

        self._clear_all_cached()

    def update_mapped_chems(self, new_chems):
        mapped_chems = self._get_mapped_chems_list(new_chems)
        self.store.write_jsonl(mapped_chems, self.chems_mapped_fn)
        self._clear_all_cached()

    def update_unmapped_chems(self, new_chems):
        unmapped_chems = self._get_unmapped_chems_list(new_chems)
        self.store.write_jsonl(unmapped_chems, self.chems_unmapped_fn)
        self._clear_all_cached()

    def update_cids_blacklist(self, cids):
        cids = set(filter(lambda cid: cid > 0, cids))
        unique_new_cids = cids - self.cids_blacklist
        if not unique_new_cids:
            return

        with open(self.cids_blacklist_fn, 'a') as f:
            for cid in unique_new_cids:
                name = self.cid_chem_map.get(cid, {}).get('cmpdname', 'Unknown')
                f.write(json.dumps({'cid': cid, 'name': name}) + '\n')
            f.flush()

        self.clear_cached_property('cids_blacklist')

    # --- Chemical utilities ---

    def get_mol_fingerprint(self, mol):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        bitstring = fp.ToBitString()
        popcount = fp.GetNumOnBits()

        bits_per_chunk = 64
        mod = 1 << bits_per_chunk

        chunks_int = []
        for i in range(0, 1024, bits_per_chunk):
            c_int = int(bitstring[i:i+bits_per_chunk], 2)
            if c_int >= (mod >> 1):
                c_int -= mod
            chunks_int.append(c_int)

        return {'bits': chunks_int, 'popcount': popcount}

    def get_mol_bertz_complexity(self, mol):
        return GraphDescriptors.BertzCT(mol)

    def get_mol_organic_mark(self, mol):
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

        except Exception:
            self.logger.log("Error generating organic mark. Assuming True")
            return True

    def strip_inchi_layers(self, inchi_str):
        if not inchi_str.startswith("InChI="):
            raise ValueError("Input does not appear to be a valid InChI string")

        cleaned = re.sub(r"/[pqib][^/]*", "", inchi_str)
        return cleaned

    def get_mol_norm_inchi(self, mol):
        try:
            return self.strip_inchi_layers(inchi.MolToInchi(mol, options="/SNon"))
        except Exception:
            return None

    def get_chem_norm_inchi(self, chem):
        return self.strip_inchi_layers(chem['inchi_snone'])

    def get_fp_tanimoto(self, fp1, fp2):
        and_pop = sum((ai & bi).bit_count() for ai, bi in zip(fp1["bits"], fp2["bits"]))
        or_pop = fp1["popcount"] + fp2["popcount"] - and_pop

        return 1.0 if or_pop == 0 else and_pop / or_pop

    def unmapped_complexity_criteria(self, chem):
        if chem['cid'] > 0:
            return True

        return chem['bertz_complexity'] <= self.unmapped_bertz_complexity_thr or chem['heavy_count'] <= self.unmapped_heavy_count_thr

    def unmapped_complexity_criteria_mol(self, mol):
        try:
            return self.get_mol_bertz_complexity(mol) <= self.unmapped_bertz_complexity_thr or mol.GetNumHeavyAtoms() <= self.unmapped_heavy_count_thr
        except Exception:
            return False

    # --- Name processing ---

    def extract_chem_cas_numbers(self, chem):
        cas_numbers = map(lambda x: re.sub(r'^cas-', '', x.strip()), chem['cmpdsynonym'])
        cas_numbers = filter(lambda x: re.fullmatch(self.CAS_PATTERN, x), cas_numbers)
        cas_numbers = list(set(cas_numbers))
        return cas_numbers

    def __chem_name_to_ascii(self, chem_name_raw):
        return chem_names.chem_name_to_ascii(chem_name_raw)

    def clean_chem_name(self, chem_name_raw):
        return chem_names.clean_chem_name(
            chem_name_raw,
            unknown_name=self.unknown_name_ph,
        )

    def normalize_chem_name(self, chem_name_raw):
        return chem_names.normalize_chem_name(
            chem_name_raw,
            unknown_name=self.unknown_name_ph,
        )

    def good_name_criteria(self, name):
        return chem_names.classify_name(
            name,
            unknown_name=self.unknown_name_ph,
        ).accepted

    # --- Compound processing ---

    def __process_chem_synonyms(self, chem):
        synonyms = chem['cmpdsynonym']
        if not synonyms:
            return False

        if not isinstance(synonyms, list):
            synonyms = [synonyms]

        filtered_synonyms = self.filtered_synonyms_by_cid.get(chem['cid'], set())
        synonyms = [
            result.clean_name
            for result in (
                chem_names.classify_name(
                    syn,
                    manual_filtered=filtered_synonyms,
                    unknown_name=self.unknown_name_ph,
                )
                for syn in synonyms
            )
            if result.accepted
        ]
        if not synonyms:
            return False

        iupac_result = chem_names.classify_name(
            chem['iupacname'],
            manual_filtered=filtered_synonyms,
            unknown_name=self.unknown_name_ph,
        )
        if iupac_result.accepted:
            synonyms = [iupac_result.clean_name] + synonyms

        if chem['cmpdname']:
            chem['cmpdname'] = chem_names.clean_synonym(
                chem['cmpdname'],
                unknown_name=self.unknown_name_ph,
            )

        cmpdname_result = chem_names.classify_name(
            chem['cmpdname'],
            manual_filtered=filtered_synonyms,
            unknown_name=self.unknown_name_ph,
        )
        if chem['cmpdname'] is None or not cmpdname_result.accepted:
            chem['cmpdname'] = synonyms[0]
        else:
            chem['cmpdname'] = cmpdname_result.clean_name
            synonyms = [cmpdname_result.clean_name] + synonyms

        synonyms = list(dict.fromkeys(synonyms))
        chem['cmpdsynonym'] = synonyms[:self.max_synonyms_thr]

        return True

    def merge_chems(self, chem_target, chem_other):
        result_synonyms = dict()
        synonyms1 = chem_target['cmpdsynonym']
        synonyms2 = chem_other['cmpdsynonym']

        if len(synonyms1) > len(synonyms2):
            synonyms1, synonyms2 = synonyms2, synonyms1

        for i in range(len(synonyms1)):
            result_synonyms[synonyms1[i]] = 0
            result_synonyms[synonyms2[i]] = 0

        for i in range(len(synonyms1), len(synonyms2)):
            result_synonyms[synonyms2[i]] = 0

        chem_target['cmpdsynonym'] = list(dict.fromkeys(result_synonyms))
        chem_target['cas'] = list(set(chem_target['cas']) | set(chem_other['cas']))

    def process_chem_single(self, chem, force=None):
        if force is None:
            force = set()

        cid = chem['cid']

        try:
            is_mapped = self.is_chem_mapped(chem)

            if chem['charge'] != 0:
                return None

            if cid in self.cids_blacklist:
                return None

            if is_mapped and chem['complexity'] > self.complexity_thr:
                return None

            if '/i' in chem['inchi']:
                return None

            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                return None
            chem['smiles'] = Chem.MolToSmiles(mol, canonical=True)

            if 'cas' not in chem:
                chem['cas'] = self.extract_chem_cas_numbers(chem)

            if 'iupacname' not in chem:
                chem['iupacname'] = None

            if is_mapped:
                if not self.__process_chem_synonyms(chem):
                    return None
            else:
                if not self.good_name_criteria(chem['cmpdname']):
                    chem['cmpdname'] = self.unknown_name_ph

                if chem['cmpdname'] == self.unknown_name_ph:
                    chem['cmpdsynonym'] = []
                else:
                    chem['cmpdsynonym'] = [chem['cmpdname']]

            def is_hydrate_inchi(inchi_str):
                try:
                    formula_part = inchi_str.split("/")[1]
                except IndexError:
                    formula_part = inchi_str[8:]

                molecules = formula_part.split(".")

                water_pattern = re.compile(r'^(\d*)H2O$', re.IGNORECASE)
                for mol_str in molecules:
                    if water_pattern.match(mol_str):
                        return True

                return None

            if is_mapped:
                if is_hydrate_inchi(chem['inchi']) and 'hydrate' in chem['cmpdname']:
                    return None

            def need_compute_field(field):
                return field not in chem or field in force

            if need_compute_field('ECFP4_fp'):
                chem['ECFP4_fp'] = self.get_mol_fingerprint(mol)

            if need_compute_field('bertz_complexity'):
                chem['bertz_complexity'] = self.get_mol_bertz_complexity(mol)

            if need_compute_field('heavy_count'):
                chem['heavy_count'] = mol.GetNumHeavyAtoms()

            if not is_mapped and not self.unmapped_complexity_criteria(chem):
                return None

            if need_compute_field('organic'):
                chem['organic'] = self.get_mol_organic_mark(mol)

            if need_compute_field('wiki'):
                chem['wiki'] = self.cid_wiki_map.get(cid, None)

            if need_compute_field('inchi_snone'):
                chem['inchi_snone'] = inchi.MolToInchi(mol, options="/SNon")

            if need_compute_field('inchikey_snone'):
                chem['inchikey_snone'] = inchi.MolToInchiKey(mol, options="/SNon")

            if not chem['inchi_snone'] or not chem['inchikey_snone']:
                return None

            return chem

        except Exception as e:
            self.logger.log_warn(f"Exception during processing compound with CID {cid}: {e}")
            return None

    def process_chems(self, force=None):
        processed_chems = []
        for chem in self.logger.track(self.chems, "Processing compounds"):
            if (chem := self.process_chem_single(chem, force=force)):
                processed_chems.append(chem)

        return processed_chems

    def organize_chems_file(self, force=None):
        initial_chems_num = len(self.chems)

        processed_chems = self.process_chems(force=force)
        self.update_chems(processed_chems)

        self.logger.print(f"Discarded {initial_chems_num - len(processed_chems)} chems")

    def build_chem(self, cid, name, smiles):
        if not name:
            name = "<Unknown>"

        if cid >= 0:
            raise ValueError(f"Invalid CID: {cid}")

        chem = dict()
        chem['cid'] = cid
        chem['cmpdname'] = name

        try:
            mol = Chem.MolFromSmiles(smiles)
            inchi_val = Chem.MolToInchi(mol)
            if not inchi_val:
                raise Exception("Failed to obtain inchi")
            inchikey = Chem.InchiToInchiKey(inchi_val)
            if not inchikey:
                raise Exception("Failed to obtain inchikey")
            mf = rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
            charge = Chem.GetFormalCharge(mol)

        except Exception as e:
            self.logger.log_warn(f"Failed process smiles '{smiles}': {e}")
            return None

        chem['smiles'] = smiles
        chem['inchi'] = inchi_val
        chem['inchikey'] = inchikey
        chem['mf'] = mf
        chem['mw'] = mw
        chem['complexity'] = None
        chem['charge'] = charge
        chem['cmpdsynonym'] = [name]

        chem = self.process_chem_single(chem)
        if not chem:
            return None

        return chem

    def add_new_chems(self, name_smiles):
        res_chems = list(self.chems)
        init_chems_num = len(res_chems)

        occupied_cids = [chem['cid'] for chem in self.cids_unmapped]
        free_cid = min(occupied_cids)-1 if occupied_cids else -1

        norm_inchi_set = set(self.norm_inchi_cid_map.keys())

        for name, smiles in self.logger.track(name_smiles, "Adding compounds"):
            chem = self.build_chem(free_cid, name, smiles)
            if chem:
                norm_inchi = self.get_chem_norm_inchi(chem)
                if norm_inchi not in norm_inchi_set:
                    res_chems.append(chem)
                    norm_inchi_set.add(norm_inchi)
                    free_cid -= 1

        self.update_chems(res_chems)

        result_chems_num = len(res_chems)
        self.logger.log(f"Processed {len(name_smiles)} entries; added {result_chems_num-init_chems_num} new compounds")

    # --- Batch operations ---

    def generate_chems_fingerprints(self):
        for chem in self.chems:
            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")
            chem['ECFP4_fp'] = self.get_mol_fingerprint(mol)

        self.update_chems(self.chems)

    def extract_chems_cas_numbers(self):
        cnt = 0
        for chem in self.chems:
            cas_numbers = self.extract_chem_cas_numbers(chem)
            chem['cas'] = cas_numbers
            if len(cas_numbers) != 0:
                cnt += 1

        self.update_chems(self.chems)
        self.logger.print(f"Extracted CAS numbers for {cnt} chems")

    def generate_chems_bertz_complexity(self):
        for chem in self.chems:
            mol = Chem.MolFromSmiles(chem['smiles'])
            if mol is None:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")
            chem['bertz_complexity'] = self.get_mol_bertz_complexity(mol)

        self.update_chems(self.chems)

    def merge_wiki_chems(self):
        for chem in self.chems:
            chem['wiki'] = self.cid_wiki_map.get(chem['cid'], None)
        self.update_chems(self.chems)

    def generate_organic_marks_for_chems(self):
        for chem in self.chems:
            smiles = chem['smiles']
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise Exception(f"Invalid smiles for {chem['cmpdname']}")
            chem['organic'] = self.get_mol_organic_mark(mol)

        self.update_chems(self.chems)

    # --- Display helpers ---

    def display_compound_table(self, compound_i, chem, cid, extra_rows=None):
        syns = chem['cmpdsynonym']
        name = chem['cmpdname']
        inchi_val = chem['inchi']
        cas = chem['cas']
        syns_num_to_disp = 10
        top_syns_str = ', '.join(f'"{syn}"' for syn in syns[:syns_num_to_disp])

        table = Table(
            title=f"[bold cyan]Compound {compound_i}: {name}[/bold cyan] (CID: [yellow]{cid}[/yellow])",
            show_header=False,
            expand=True
        )

        if extra_rows:
            for row in extra_rows:
                table.add_row(*row)

        table.add_row(f"First {syns_num_to_disp} synonyms", f"{top_syns_str}")
        table.add_row("InChI", f"{inchi_val}")
        table.add_row("CAS", f"{cas}")
        self.logger.print(table)

    def display_compare_table(self, title, chem1, chem2, display_compound=None):
        if display_compound is None:
            display_compound = self.display_compound_table

        self.logger.print(Rule(title))

        display_compound(1, chem1, chem1['cid'])
        self.logger.print(Rule())
        display_compound(2, chem2, chem2['cid'])

    # --- Utility ---

    def count_file_lines(self, filename):
        count = 0
        with open(filename, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                count += chunk.count(b"\n")
        return count
