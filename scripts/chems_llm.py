from chems_parse_reactions import ChemsParseReactions



class ChemsLLM(ChemsParseReactions):
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.gpt_oss = "openai/gpt-oss-120b"
        self.qwen = "qwen/qwen3-235b-a22b"
        self.grok = "x-ai/grok-3-mini"
        self.gemini = "google/gemini-2.5-flash-lite"


