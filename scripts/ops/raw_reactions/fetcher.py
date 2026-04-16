import json

from .prompts import build_fetch_prompt, build_revalidate_prompt
from .presets import POSITION_ANY, POSITION_PRODUCT, POSITION_REAGENT
from .schema import is_valid_compound, is_valid_reaction_obj, normalize_reaction_obj
from scripts.infra.batch_runner import run_batch
from scripts.infra.fallback import with_fallback


_is_valid_compound = is_valid_compound
_is_valid_reaction_obj = is_valid_reaction_obj
_normalize_reaction_obj = normalize_reaction_obj


def _count(stats, key):
    if stats is not None:
        stats[key] = stats.get(key, 0) + 1


def parse_reactions_jsonl(response, stats=None):
    reactions = []
    for line in response.strip().split('\n'):
        line = line.strip()
        if not line or line.lower() == 'null':
            _count(stats, 'null_or_empty')
            continue
        try:
            obj = json.loads(line)
        except (json.JSONDecodeError, TypeError):
            _count(stats, 'malformed_json')
            continue
        if not is_valid_reaction_obj(obj):
            _count(stats, 'invalid_schema')
            continue
        reactions.append(normalize_reaction_obj(obj))
        _count(stats, 'accepted')
    return reactions


class RawReactionsFetcher:
    MAX_EXISTING_CONTEXT = 20

    def __init__(self, compounds, llm_client, store, logger, layout, models):
        self.compounds = compounds
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.models = models
        self._existing_by_cid = {}

    def set_existing_reactions(self, existing_by_cid):
        self._existing_by_cid = existing_by_cid or {}

    def _get_existing_context(self, cid, position):
        by_position = self._existing_by_cid.get(cid, {})
        if position == POSITION_REAGENT:
            reactions = by_position.get(POSITION_REAGENT, [])
        elif position == POSITION_PRODUCT:
            reactions = by_position.get(POSITION_PRODUCT, [])
        elif position == POSITION_ANY:
            reactions = []
            seen = set()
            for side in (POSITION_REAGENT, POSITION_PRODUCT):
                for reaction in by_position.get(side, []):
                    if reaction in seen:
                        continue
                    seen.add(reaction)
                    reactions.append(reaction)
        else:
            raise ValueError(f"Unknown position: {position!r}")
        return reactions[:self.MAX_EXISTING_CONTEXT]

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
        existing = self._get_existing_context(chem['cid'], preset.position)
        fetch_prompt = build_fetch_prompt(
            chem_name, preset.position, preset.scope, existing,
        )

        response = self.llm_client.fetch_answer_str(fetch_prompt, model)
        stats = {}
        reactions = parse_reactions_jsonl(response, stats)
        self._log_parse_stats(chem_name, f"Fetch response from '{model}'", stats)

        if not reactions:
            return None

        reactions_jsonl_str = '\n'.join(json.dumps(r) for r in reactions)
        revalidate_prompt = build_revalidate_prompt(reactions_jsonl_str)
        revalidated_response = self.llm_client.fetch_answer_str(revalidate_prompt, model)
        stats = {}
        revalidated = parse_reactions_jsonl(revalidated_response, stats)
        self._log_parse_stats(chem_name, "Revalidation response", stats)

        if not revalidated:
            self.logger.log_warn(
                f"Revalidation returned no results for '{chem_name}'. Using original."
            )
            revalidated = reactions

        self.logger.log(
            f"Got {len(revalidated)} reactions for '{chem_name}' with '{model}'; "
            f"CTT: {self.llm_client.completion_tokens_total}"
        )
        return {'cid': chem['cid'], 'reactions': revalidated}

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
