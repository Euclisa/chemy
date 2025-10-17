import json
import re
import shutil
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from chempy import balance_stoichiometry
import hashlib
import base64

from chems_pubchem_parse import ChemsParsePubchem

# Disable all RDKit warnings and info messages
RDLogger.DisableLog('rdApp.*')


class ChemsParseReactions(ChemsParsePubchem):

    def __init__(self, data_dir):
        super().__init__(data_dir)

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
    

    def process_raw_reactions(self, reactions_fn):
        with open(reactions_fn) as f:
            entries = [json.loads(x) for x in f.read().strip().split('\n')]
        
        for entry in entries:
            reactions = entry['reactions']
            for react in reactions:
                react = react.strip()
                reagents, products = react.split('->')
                reagents = reagents.strip().split('+')
                products = products.strip().split('+')


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
        if max_coeff > 30:
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
    

    def find_all_unicode_chars_in_raw_reactions(self):
        non_ascii = dict()
        with open(self.raw_reactions_verdict_fn) as f:
            for line in f:
                reaction = json.loads(line)['reaction']
                non_ascii_curr = [char for char in reaction if ord(char) > 127]
                for char in non_ascii_curr:
                    if char not in non_ascii:
                        non_ascii[char] = 0
                    non_ascii[char] += 1
        
        return non_ascii

    
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


    
    def _parse_reaction_str(self, reaction_str: str):
        parts = reaction_str.split('->')
        if len(parts) != 2:
            return None, set()

        reagents_str, products_str = parts
        parse_success = True
        unmapped_names = set()

        def parse_compounds(compound_str, skip_names, existing_cids=None):
            if existing_cids is None:
                existing_cids = set()

            compounds, cids = [], set()
            for name in compound_str.split('+'):
                norm = self._normalize_chem_name(name)
                if norm in skip_names:
                    continue

                clean = self._clean_chem_name(name)
                cid = self.name_cid_map.get(norm)
                if cid is None:
                    nonlocal parse_success
                    parse_success = False
                    unmapped_names.add((norm, clean))

                if cid is None or cid not in existing_cids | cids:
                    compounds.append({'norm_name': norm, 'original_name': clean, 'cid': cid})
                    cids.add(cid)
            
            compounds = list({c["cid"]: c for c in compounds}.values())
            if parse_success:
                compounds.sort(key=lambda c: c['cid'])

            return compounds, cids

        reagents, reagents_cids = parse_compounds(reagents_str, {"light", "heat", "catalyst"})
        products, products_cids = parse_compounds(products_str, {"otherproducts"}, reagents_cids)

        if products_cids & reagents_cids or not products or not reagents:
            parse_success = False

        if not parse_success:
            return None, unmapped_names

        reaction = {'reagents': reagents, 'products': products}
        reaction = self._assemble_reaction(reaction)

        return reaction, unmapped_names
