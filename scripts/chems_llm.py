import os

from chems_parse_reactions import ChemsParseReactions



class ChemsLLM(ChemsParseReactions):
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.gpt_oss = "openai/gpt-oss-120b"
        self.qwen = "qwen/qwen3-235b-a22b"
        self.grok = "x-ai/grok-3-mini"
        self.gemini = "google/gemini-2.5-flash-lite"

        self.raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions.jsonl")
        self.wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "wiki_raw_reactions.jsonl")
        self.top_rare_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', "top_rare_raw_reactions.jsonl")
        self.raw_reactions_verdict_fn = os.path.join(self.data_dir, 'raw_reactions', "raw_reactions_verdict.jsonl")
        self.products_wiki_raw_reactions_fn = os.path.join(self.data_dir, 'raw_reactions', 'wiki_products_raw_reactions.jsonl')


