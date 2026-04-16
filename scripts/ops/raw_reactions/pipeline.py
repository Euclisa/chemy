from .fetcher import RawReactionsFetcher, _is_valid_reaction_obj
from .layout import RawReactionsLayout
from .presets import PRESETS, get_preset
from .validator import RawReactionsValidator


class RawReactionsPipeline:
    def __init__(self, data_dir, compounds, reactions, store, logger, llm_client, parser):
        self.layout = RawReactionsLayout(data_dir)
        self._compounds = compounds
        self._reactions = reactions
        self._store = store
        self._logger = logger
        self._parser = parser

        models = [parser.gpt_oss, parser.qwen]
        self._fetcher = RawReactionsFetcher(compounds, llm_client, store, logger, self.layout, models)
        self._validator = RawReactionsValidator(llm_client, store, logger, self.layout, parser, models)

    def _build_existing_reactions_context(self):
        llm_reactions = self._store.load_jsonl(self._parser.reactions_parsed_llm_fn)

        by_cid = {}
        for react in llm_reactions:
            try:
                reaction_str = self._reactions.get_reaction_as_str(react)
            except (KeyError, TypeError):
                continue
            for entry in react['reagents'] + react['products']:
                by_cid.setdefault(entry['cid'], []).append(reaction_str)

        return by_cid

    def fetch(self, preset_name, max_workers=1):
        existing = self._build_existing_reactions_context()
        self._fetcher.set_existing_reactions(existing)
        self._fetcher.fetch_all(get_preset(preset_name), max_workers)

    def validate(self, preset_name, max_workers=1):
        self._validator.validate(get_preset(preset_name), max_workers)

    def run(self, preset_name, max_workers=1):
        self.fetch(preset_name, max_workers)
        self.validate(preset_name, max_workers)

    def run_many(self, preset_names, max_workers=1):
        for name in preset_names:
            self.run(name, max_workers)

    def parse(self, preset_names=None):
        if preset_names is None:
            preset_names = list(PRESETS)

        all_entries = []
        for name in preset_names:
            all_entries.extend(self._store.load_jsonl(self.layout.verdict(name)))

        parsed = []
        processed_rids = set()
        skip_stats = {
            'validator_reject': 0,
            'invalid_schema': 0,
            'unmapped': 0,
            'duplicate_rid': 0,
        }
        for entry in self._logger.track(all_entries, "Parsing structured reactions"):
            if not entry.get('valid', False):
                skip_stats['validator_reject'] += 1
                continue

            reaction_obj = entry.get('reaction')
            if not _is_valid_reaction_obj(reaction_obj):
                skip_stats['invalid_schema'] += 1
                continue
            parsed_reaction, _ = self._parser.parse_structured_reaction(reaction_obj)
            if not parsed_reaction:
                skip_stats['unmapped'] += 1
                continue
            if parsed_reaction['rid'] in processed_rids:
                skip_stats['duplicate_rid'] += 1
                continue

            processed_rids.add(parsed_reaction['rid'])
            parsed_reaction['confidence'] = entry.get('confidence', 0.0)
            parsed_reaction['source'] = entry.get('source')
            parsed.append(parsed_reaction)

        self._store.write_jsonl(parsed, self._parser.reactions_parsed_llm_fn)
        if any(skip_stats.values()):
            self._logger.log_warn(f"Parsing skipped reactions: {skip_stats}")
        self._logger.log(f"Successfully parsed {len(parsed)} reactions!")

    def find_unicode_chars(self, preset_name):
        non_ascii = {}
        for entry in self._store.load_jsonl(self.layout.verdict(preset_name)):
            reaction_obj = entry.get('reaction', {})
            if not isinstance(reaction_obj, dict):
                continue
            for compound in reaction_obj.get('reagents', []) + reaction_obj.get('products', []):
                for field in ('name', 'note'):
                    for char in compound.get(field, ''):
                        if ord(char) > 127:
                            non_ascii[char] = non_ascii.get(char, 0) + 1
        return non_ascii
