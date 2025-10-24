from chempy import balance_stoichiometry
import hashlib
import base64
import os

from chems_pubchem_parse import ChemsParsePubchem



class ChemsParseReactions(ChemsParsePubchem):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.reactions_parsed_fn = os.path.join(self.data_dir, 'reactions_parsed', "reactions_parsed.jsonl")
        self.reactions_parsed_llm_fn = os.path.join(self.data_dir, 'reactions_parsed', "reactions_parsed_llm.jsonl")
        self.reactions_parsed_ord_fn = os.path.join(self.data_dir, 'reactions_parsed', 'reactions_parsed_ord.jsonl')
        self.reactions_parsed_fixed_fn = os.path.join(self.data_dir, 'reactions_parsed', 'reactions_parsed_fixed.jsonl')

        self.reactions_details_fn = os.path.join(self.data_dir, 'reactions_details', 'reactions_details.jsonl')

        self.unmapped_names_fn = os.path.join(self.data_dir, "unmapped_names.jsonl")
        self.chem_names_blacklisted_fn = os.path.join(self.data_dir, "unmapped_names_blacklisted.txt")
        self.unmapped_smiles_fn = os.path.join(self.data_dir, 'unmapped_smiles.jsonl')
        self.chem_smiles_blacklisted_fn = os.path.join(self.data_dir, 'unmapped_smiles_blacklisted.txt')

        self._file_sorting_prefs[self.reactions_parsed_fn] = 'rid'
        self._file_sorting_prefs[self.reactions_parsed_llm_fn] = 'rid'
        self._file_sorting_prefs[self.reactions_parsed_ord_fn] = 'rid'

        self._file_sorting_prefs[self.reactions_details_fn] = 'rid'

        self._file_sorting_prefs[self.unmapped_names_fn] = ('count', True)
        self._file_sorting_prefs[self.unmapped_smiles_fn] = ('count', True)

        self.sources_priority = dict()
    

    def _get_reaction_as_str(self, reaction):
        def format_components(components):
            parts = []
            for c in components:
                coeff = f"{c['coeff']} " if c.get('coeff') is not None else ""
                parts.append(f"{coeff}{c['original_name']}")
            return " + ".join(parts)

        reagents_str = format_components(reaction['reagents'])
        products_str = format_components(reaction['products'])

        return f"{reagents_str} -> {products_str}"


    def _balance_reaction(self, reaction):
        check_rid = self._get_reaction_hash(reaction)
        if reaction['rid'] != check_rid:
            raise Exception(f"Bad RID for reaction '{reaction['rid']}'")

        reagents = [self.cid_mf_map[x['cid']] for x in reaction['reagents']]
        products = [self.cid_mf_map[x['cid']] for x in reaction['products']]

        try:
            reagents_coeffs, products_coeffs = balance_stoichiometry(reagents, products, underdetermined=False)
        except:
            reaction['balanced'] = False
            return

        all_coeffs = list(reagents_coeffs.values()) + list(products_coeffs.values())
        all_coeffs = [int(x) for x in all_coeffs]
        max_coeff = max(all_coeffs)
        if max_coeff > 15:
            reaction['balanced'] = False
            return

        for chem in reaction['reagents']:
            mf = self.cid_mf_map[chem['cid']]
            chem['coeff'] = int(reagents_coeffs[mf])
        
        for chem in reaction['products']:
            mf = self.cid_mf_map[chem['cid']]
            chem['coeff'] = int(products_coeffs[mf])
        
        reaction['balanced'] = True


    def balance_parsed_reactions(self, reactions_parsed_fn=None):
        if reactions_parsed_fn is None:
            reactions_parsed_fn = self.reactions_parsed_fn

        reactions = self._load_jsonl(reactions_parsed_fn)
        
        balanced_cnt = 0
        for react in reactions:
            self._balance_reaction(react)
            balanced_cnt += react['balanced']
        

        self._write_jsonl(reactions, reactions_parsed_fn)        
        self.log(f"Balanced {balanced_cnt} out of {len(reactions)}")

    
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
        for chem in reaction['products']:
            av_complexity += self.cid_chem_map[chem['cid']]['complexity']
        
        for chem in reaction['reagents']:
            av_complexity += self.cid_chem_map[chem['cid']]['complexity']
        av_complexity /= len(reaction['reagents']) + len(reaction['products'])

        return av_complexity


    def _assemble_reaction(self, reaction):
        reaction['complexity'] = self._get_reaction_complexity(reaction)
        reaction['rid'] = self._get_reaction_hash(reaction)

        return reaction


    def _get_parsed_reactions_participants_norm_names(self):
        reactions = self._load_jsonl(self.reactions_parsed_fn)

        norm_names = set()
        for react in reactions:
            norm_names.update(entry['norm_name'] for entry in react['reagents']+react['products'])
        
        return norm_names


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

