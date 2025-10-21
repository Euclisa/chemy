import json
import os

from chems_parse_reactions import ChemsParseReactions


class ChemsLLMParse(ChemsParseReactions):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.gpt_oss = "openai/gpt-oss-120b"
        self.qwen = "qwen/qwen3-235b-a22b"
        self.grok = "x-ai/grok-3-mini"
        self.gemini = "google/gemini-2.5-flash-lite"
        self.deepseek = "deepseek/deepseek-v3.2-exp"

        self.raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions.jsonl")
        self.wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "wiki_raw_reactions.jsonl")
        self.top_rare_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "top_rare_raw_reactions.jsonl")
        self.raw_reactions_verdict_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions_verdict.jsonl")
        self.products_wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', 'wiki_products_raw_reactions.jsonl')
        self.products_annot_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', 'annotated_products_raw_reactions.jsonl')

        self.reactions_thermo_llm_fn = os.path.join(self.data_dir, 'thermo', 'llm', 'reactions_thermo_llm.jsonl')

        self._file_sorting_prefs[self.reactions_thermo_llm_fn] = 'rid'

        self.sources_priority[self.gpt_oss] = 5
        self.sources_priority[self.qwen] = 4
        self.sources_priority[self.grok] = 3
    

    def _parse_reaction_scheme(self, reaction_str: str, balance=True):
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

        if balance:
            self._balance_reaction(reaction)

        return reaction, unmapped_names
    

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


    def parse_raw_llm_reactions(self, balance=True):
        unmapped_names = dict()
        parsed = []
        processed_rids = set()
        with open(self.raw_reactions_verdict_fn) as f:
            total = self._count_file_lines(self.raw_reactions_verdict_fn)
            for line in self._rich_track(f, "Parsing raw LLM reactions", total=total):
                entry = json.loads(line)
                reaction = entry['reaction']
                confidence = entry['confidence']
                if confidence < 0.4:
                    continue
                parsed_reaction, unmapped_names_curr = self._parse_reaction_scheme(reaction, balance)
                if not parsed_reaction:
                    continue
                
                processed_rids.add(parsed_reaction['rid'])
                
                parsed_reaction['confidence'] = confidence
                parsed_reaction['source'] = entry['source']

                parsed.append(parsed_reaction)
        
        self._write_jsonl(parsed, self.reactions_parsed_llm_fn)
        
        self.log(f"Successfully parsed {len(parsed)} reactions!")




if __name__ == "__main__":
    chems_parse = ChemsLLMParse('data/')
    chems_parse.parse_raw_llm_reactions()