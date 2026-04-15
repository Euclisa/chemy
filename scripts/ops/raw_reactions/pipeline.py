from .fetcher import RawReactionsFetcher
from .layout import RawReactionsLayout
from .presets import PRESETS, get_preset
from .validator import RawReactionsValidator


class RawReactionsPipeline:
    def __init__(self, data_dir, compounds, reactions, store, logger, llm_client, parser):
        self.layout = RawReactionsLayout(data_dir)
        self._store = store
        self._logger = logger
        self._parser = parser

        models = [parser.gpt_oss, parser.qwen]
        self._fetcher = RawReactionsFetcher(compounds, llm_client, store, logger, self.layout, models)
        self._validator = RawReactionsValidator(llm_client, store, logger, self.layout, parser, models)

    def fetch(self, preset_name: str, max_workers: int = 1):
        self._fetcher.fetch_all(get_preset(preset_name), max_workers)

    def validate(self, preset_name: str, max_workers: int = 1):
        self._validator.validate(get_preset(preset_name), max_workers)

    def run(self, preset_name: str, max_workers: int = 1):
        """Fetch raw reactions then validate them in one step."""
        self.fetch(preset_name, max_workers)
        self.validate(preset_name, max_workers)

    def run_many(self, preset_names, max_workers: int = 1):
        for name in preset_names:
            self.run(name, max_workers)

    def parse(self, preset_names=None):
        """Union verdicts from the given presets and write to reactions_parsed_llm_fn."""
        if preset_names is None:
            preset_names = list(PRESETS)

        all_entries = []
        for name in preset_names:
            all_entries.extend(self._store.load_jsonl(self.layout.verdict(name)))

        parsed = []
        processed_rids = set()
        for entry in self._logger.track(all_entries, "Parsing raw LLM reactions"):
            if entry['confidence'] < 0.4:
                continue
            parsed_reaction, _ = self._parser.parse_reaction_scheme(entry['reaction'])
            if not parsed_reaction or parsed_reaction['rid'] in processed_rids:
                continue
            processed_rids.add(parsed_reaction['rid'])
            parsed_reaction['confidence'] = entry['confidence']
            parsed_reaction['source'] = entry['source']
            parsed.append(parsed_reaction)

        self._store.write_jsonl(parsed, self._parser.reactions_parsed_llm_fn)
        self._logger.log(f"Successfully parsed {len(parsed)} reactions!")

    def find_unicode_chars(self, preset_name: str) -> dict:
        """Diagnostic: return all non-ASCII chars found in a preset's verdict file."""
        non_ascii = {}
        for entry in self._store.load_jsonl(self.layout.verdict(preset_name)):
            for char in entry['reaction']:
                if ord(char) > 127:
                    non_ascii[char] = non_ascii.get(char, 0) + 1
        return non_ascii
