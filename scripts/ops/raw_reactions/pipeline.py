from .fetcher import RawReactionsFetcher
from .layout import RawReactionsLayout
from .models import DEFAULT_MODELS
from .presets import get_preset
from .repairer import RawReactionsRepairer
from .schema import is_valid_reaction_obj
from .services import RawReactionTaskBuilder
from .validator import ACCEPT_THR, RawReactionsValidator, is_confident


class RawReactionsPipeline:
    def __init__(self, data_dir, compounds, reactions, store, logger, llm_client, parser):
        self.layout = RawReactionsLayout(data_dir)
        self._compounds = compounds
        self._reactions = reactions
        self._store = store
        self._logger = logger
        self._parser = parser
        self._task_builder = RawReactionTaskBuilder(reactions)

        store.register_vault(self.layout.verdict(), 'verdict_')

        self.models = DEFAULT_MODELS
        self._fetcher = RawReactionsFetcher(
            compounds, llm_client, store, logger, self.layout, self.models,
        )
        self._validator = RawReactionsValidator(
            llm_client, store, logger, self.layout, parser, self.models,
        )
        self._repairer = RawReactionsRepairer(
            llm_client,
            store,
            logger,
            self.layout,
            parser,
            self._validator,
            self.models,
        )

    def _build_existing_reactions_context(self):
        llm_reactions = self._store.load_jsonl(self._parser.reactions_parsed_llm_fn)
        return self._task_builder.build_existing_reactions_context(llm_reactions)

    def fetch(self, preset_name, max_workers=1, run=None):
        existing = self._build_existing_reactions_context()
        self._fetcher.set_existing_reactions(existing)
        self._fetcher.fetch_all(get_preset(preset_name), max_workers, run=run)

    def validate(self, preset_name, max_workers=1, acceptance_threshold=ACCEPT_THR):
        self._validator.validate(
            get_preset(preset_name),
            max_workers,
            acceptance_threshold=acceptance_threshold,
        )

    def repair(self, max_workers=1, lower_bound=0.1, acceptance_threshold=ACCEPT_THR):
        self._repairer.repair(
            max_workers=max_workers,
            lower_bound=lower_bound,
            acceptance_threshold=acceptance_threshold,
        )

    def run(self, preset_name, max_workers=1, acceptance_threshold=ACCEPT_THR, run=None):
        self.fetch(preset_name, max_workers, run=run)
        self.validate(preset_name, max_workers, acceptance_threshold=acceptance_threshold)

    def parse(self, acceptance_threshold=ACCEPT_THR):
        all_entries = []
        selected_rids = set()
        accepted_original_rids = set()

        for entry in self._store.load_jsonl(self.layout.verdict()):
            rid = entry.get('rid')
            if rid is None:
                reaction_obj = entry.get('reaction')
                if is_valid_reaction_obj(reaction_obj):
                    parsed_reaction, _ = self._parser.parse_structured_reaction(reaction_obj)
                    if parsed_reaction:
                        rid = parsed_reaction['rid']
            if rid is not None:
                selected_rids.add(rid)
                if is_confident(entry, acceptance_threshold=acceptance_threshold):
                    accepted_original_rids.add(rid)
            all_entries.append(entry)

        for entry in self._store.load_jsonl(self.layout.repairs()):
            if entry.get('rid') not in selected_rids:
                continue
            if entry.get('rid') in accepted_original_rids:
                continue
            if entry.get('fixed_reaction') is None:
                continue
            if not is_confident(
                entry,
                key='validation_confidence',
                acceptance_threshold=acceptance_threshold,
            ):
                continue
            all_entries.append({
                'reaction': entry['fixed_reaction'],
                'confidence': entry.get('validation_confidence', 0.0),
                'source': (
                    f"repair:{entry.get('repair_source')}"
                    f"|validate:{entry.get('validation_source')}"
                ),
            })

        parsed = []
        processed_rids = set()
        skip_stats = {
            'confidence_reject': 0,
            'invalid_schema': 0,
            'unmapped': 0,
            'duplicate_rid': 0,
        }
        for entry in self._logger.track(all_entries, "Parsing structured reactions"):
            if not is_confident(entry, acceptance_threshold=acceptance_threshold):
                skip_stats['confidence_reject'] += 1
                continue

            reaction_obj = entry.get('reaction')
            if not is_valid_reaction_obj(reaction_obj):
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
