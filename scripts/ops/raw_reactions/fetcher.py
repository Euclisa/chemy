from .presets import POSITION_ANY, POSITION_PRODUCT, POSITION_REAGENT
from .schema import is_valid_compound, is_valid_reaction_obj, normalize_reaction_obj
from .services import (
    RawReactionGenerationService,
    RawReactionTaskBuilder,
    parse_reactions_jsonl,
)
from scripts.infra.batch_runner import run_batch
from scripts.infra.fallback import with_fallback


_is_valid_compound = is_valid_compound
_is_valid_reaction_obj = is_valid_reaction_obj
_normalize_reaction_obj = normalize_reaction_obj


class RawReactionsFetcher:
    MAX_EXISTING_CONTEXT = 20

    def __init__(self, compounds, llm_client, store, logger, layout, models):
        self.compounds = compounds
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.models = models
        self.task_builder = RawReactionTaskBuilder()
        self.generator = RawReactionGenerationService(llm_client)

    def set_existing_reactions(self, existing_by_cid):
        self.task_builder.set_existing_reactions(existing_by_cid)

    def _get_existing_context(self, cid, position):
        return self.task_builder.get_existing_context(cid, position)

    def _log_parse_stats(self, chem_name, label, stats):
        skipped = sum(
            stats.get(key, 0)
            for key in ('malformed_json', 'invalid_schema', 'null_or_empty')
        )
        if skipped:
            self.logger.log_warn(
                f"{label} skipped {skipped} lines for '{chem_name}': {stats}"
            )

    def fetch_one(self, chem, preset, model):
        """Fetch reactions for one compound using a single model.

        Returns a result dict with 'cid' and 'reactions' on success, or None
        when the model returns no usable reactions.
        """
        chem_name = chem['cmpdname']
        task = self.task_builder.build_task(chem, preset)
        result = self.generator.generate(task, model, self_review=True)
        stats = result.parse_stats
        self._log_parse_stats(chem_name, f"Fetch response from '{model}'", stats)

        if result.reviewed_parse_stats:
            self._log_parse_stats(
                chem_name,
                f"Revalidation response from '{model}'",
                result.reviewed_parse_stats,
            )

        if result.error is not None:
            self.logger.log_warn(
                f"Fetch failed for '{chem_name}' with '{model}': {result.error}"
            )
            return None

        reactions = list(result.reactions)
        if not reactions:
            return None

        self.logger.log(
            f"Got {len(reactions)} reactions for '{chem_name}' with '{model}'; "
            f"CTT: {self.llm_client.completion_tokens_total}"
        )
        return {'cid': chem['cid'], 'reactions': reactions}

    def _fetch_one_with_fallback(self, chem, preset):
        result = with_fallback(
            lambda model: self.fetch_one(chem, preset, model),
            self.models,
            logger=self.logger,
        )
        if result is not None:
            return result
        self.logger.log_warn(f"Failed to fetch reactions for '{chem['cmpdname']}'")
        return {'cid': chem['cid'], 'reactions': None}

    def fetch_all(self, preset, max_workers=1, run=None):
        raw_fn = self.layout.raw(preset.name, run=run)
        processed = {x['cid'] for x in self.store.load_jsonl(raw_fn)}
        criteria = preset.build_criteria(self.compounds)
        staged = [
            chem for chem in self.compounds.chems
            if chem['cid'] not in processed and criteria(chem)
        ]
        run_batch(
            self.llm_client,
            staged,
            lambda chem: self._fetch_one_with_fallback(chem, preset),
            raw_fn,
            self.logger,
            max_workers=max_workers,
            description=f"Fetching reactions [{preset.name}]",
        )
