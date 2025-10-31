from chempy import balance_stoichiometry
import hashlib
import base64
import os

from functools import cached_property
from types import MappingProxyType

from chems_properties import ChemsProperties



class ChemsReactionProperties(ChemsProperties):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.parsed_reactions_dir = os.path.join(self.data_dir, 'reactions_parsed')
        self.reactions_details_dir = os.path.join(self.data_dir, 'reactions_details')

        self.reactions_balance_fn = os.path.join(self.parsed_reactions_dir, 'reactions_balance')

        self.unmapped_names_fn = os.path.join(self.data_dir, "unmapped_names.jsonl")
        self.chem_names_blacklisted_fn = os.path.join(self.data_dir, "unmapped_names_blacklisted.txt")
        self.chem_smiles_blacklisted_fn = os.path.join(self.data_dir, 'unmapped_smiles_blacklisted.txt')

        self._file_sorting_prefs[self.reactions_balance_fn] = 'rid'
        self._dir_vault_prefs[self.reactions_balance_fn] = 'balance'

        self._file_sorting_prefs[self.unmapped_names_fn] = ('count', True)

        self._max_balance_coeff_thr = 15

        self.sources_priority = dict()
    

    @cached_property
    def reactions_balance(self):
        balance = self._load_jsonl(self.reactions_balance_fn)
        return MappingProxyType({entry['rid']: self._entries_to_map(entry['balance'], 'cid', 'coeff') for entry in balance})
    

    def _get_all_reaction_cids(self, reaction):
        return list(entry['cid'] for entry in reaction['reagents']+reaction['products'])
    

    def _get_part_cids(self, part):
        return [entry['cid'] for entry in part]
    

    def _validate_reaction(self, reaction):
        all_cids = self._get_all_reaction_cids(reaction)

        if len(all_cids) != len(set(all_cids)):
            return False
        
        reagents_cids, products_cids = self._get_part_cids(reaction['reagents']), self._get_part_cids(reaction['products'])
        if not reagents_cids or not products_cids:
            return False
        
        if self._get_reaction_hash(reaction) != reaction['rid']:
            return False

        return True
    

    def _get_reaction_as_str(self, reaction):
        balance = self.reactions_balance.get(reaction['rid'])

        def format_components(components):
            parts = []
            for c in components:
                cid = c['cid']
                coeff = f"{balance[cid]} " if balance else ""
                parts.append(f"{coeff}{c['original_name']}")
            return " + ".join(parts)

        reagents_str = format_components(reaction['reagents'])
        products_str = format_components(reaction['products'])

        return f"{reagents_str} -> {products_str}"


    def _remove_reaction_balance_legacy_data(self, reaction):
        if 'balanced' in reaction:
            reaction.pop('balanced')

        for entry in reaction['reagents']+reaction['products']:
            if 'coeff' in entry:
                entry.pop('coeff')
        
        return reaction

    def _balance_reaction(self, reaction):
        if not self._validate_reaction(reaction):
            reaction_str = self._get_reaction_as_str(reaction)
            raise ValueError(f"Reaction is invalid: '{reaction_str}'")
        
        self._remove_reaction_balance_legacy_data(reaction)

        res_balance = {'rid': reaction['rid'], 'balance': []}

        reagents = [self.cid_mf_map[x['cid']] for x in reaction['reagents']]
        products = [self.cid_mf_map[x['cid']] for x in reaction['products']]

        try:
            reagents_coeffs, products_coeffs = balance_stoichiometry(reagents, products, underdetermined=False)
        except:
            return None

        all_coeffs = reagents_coeffs | products_coeffs
        all_coeffs = {mf: int(coeff) for mf, coeff in all_coeffs.items()}
        max_coeff = max(x for x in all_coeffs.values())
        if max_coeff > self._max_balance_coeff_thr or any(coeff < 0 for coeff in all_coeffs.values()):
            return None

        for cid in self._get_all_reaction_cids(reaction):
            mf = self.cid_mf_map[cid]
            res_balance['balance'].append({'cid': cid, 'coeff': all_coeffs[mf]})
        
        return res_balance


    def _is_react_balanced(self, reaction):
        return reaction['rid'] in self.reactions_balance

    
    def _get_reaction_hash(self, reaction):
        reagents_cids = sorted([x['cid'] for x in reaction['reagents']])
        products_cids = sorted([x['cid'] for x in reaction['products']])
        reagents_str = '(' + ','.join([str(x) for x in reagents_cids]) + ')'
        products_str = '(' + ','.join([str(x) for x in products_cids]) + ')'
        reaction_enc = (reagents_str + products_str).encode("utf-8")
        hash_bytes = hashlib.sha256(reaction_enc).digest()
        hash_b64 = base64.b64encode(hash_bytes[:16]).decode("utf-8")

        return hash_b64
    

    def _get_reaction_complexity(self, reaction):
        av_complexity = 0
        for chem in reaction['reagents']+reaction['products']:
            av_complexity += self.cid_chem_map[chem['cid']]['bertz_complexity']
        av_complexity /= len(reaction['reagents']) + len(reaction['products'])

        return av_complexity


    def _assemble_reaction(self, reaction):
        reaction['complexity'] = self._get_reaction_complexity(reaction)
        reaction['rid'] = self._get_reaction_hash(reaction)

        return reaction


    def _convert_details_to_canonic(self, details: dict):
        required_fields = ['rid', 'source']
        for field in required_fields:
            if field not in details:
                raise Exception(f"Field '{field}' must be present in reaction details entry")
        
        details.setdefault('confidence', None)
        details.setdefault('description', None)
        details.setdefault('catalysts', [])
        details.setdefault('solvents', [])
        
        if 'provenance' not in details or details['provenance'] is None:
            details['provenance'] = {'doi': None, 'patent': None}
        
        return details



if __name__ == "__main__":
    chems_parse = ChemsReactionProperties('data/')
    chems_parse.balance_parsed_reactions()