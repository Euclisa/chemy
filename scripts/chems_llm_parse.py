import json
from rich.progress import track

from chems_llm import ChemsLLM


class ChemsLLMParse(ChemsLLM):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.sources_priority[self.gpt_oss] = 5
        self.sources_priority[self.qwen] = 4
        self.sources_priority[self.grok] = 3
    

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
            for line in self._rich_track(f, "Parsing...", total=total):
                entry = json.loads(line)
                reaction = entry['reaction']
                confidence = entry['confidence']
                if confidence < 0.4:
                    continue
                parsed_reaction, unmapped_names_curr = self._parse_reaction_str(reaction)
                for norm_name, name in unmapped_names_curr:
                    if norm_name not in unmapped_names:
                        unmapped_names[norm_name] = [0, name]
                    unmapped_names[norm_name][0] += 1
                if not parsed_reaction:
                    continue
                
                processed_rids.add(parsed_reaction['rid'])
                
                parsed_reaction['confidence'] = confidence
                parsed_reaction['source'] = entry['source']

                if balance:
                    self._balance_reaction(parsed_reaction)

                parsed.append(parsed_reaction)
        
        self._write_jsonl(parsed, self.reactions_parsed_llm_fn)
        
        self.log(f"Successfully parsed {len(parsed)} reactions!")




if __name__ == "__main__":
    chems_parse = ChemsLLMParse('data/')
    chems_parse.parse_raw_llm_reactions()