import json

from .prompts import VALID_PHASES, build_fetch_prompt, build_revalidate_prompt


REACTION_KEYS = frozenset({'reagents', 'products', 'solvent', 'catalyst'})
COMPOUND_KEYS = frozenset({'name', 'phase', 'note'})


def _count(stats, key):
    if stats is not None:
        stats[key] = stats.get(key, 0) + 1


def _is_valid_compound(compound):
    if not isinstance(compound, dict):
        return False
    if not set(compound).issubset(COMPOUND_KEYS):
        return False
    if not isinstance(compound.get('name'), str) or not compound['name'].strip():
        return False
    if compound.get('phase') not in VALID_PHASES:
        return False
    if 'note' in compound and (
        not isinstance(compound['note'], str) or not compound['note'].strip()
    ):
        return False
    return True


def _is_valid_reaction_obj(obj):
    if not isinstance(obj, dict):
        return False
    if not set(obj).issubset(REACTION_KEYS):
        return False
    for key in ('reagents', 'products'):
        compounds = obj.get(key)
        if not isinstance(compounds, list) or not compounds:
            return False
        if not all(_is_valid_compound(c) for c in compounds):
            return False
    for key in ('solvent', 'catalyst'):
        if key in obj and obj[key] is not None and (
            not isinstance(obj[key], str) or not obj[key].strip()
        ):
            return False
    return True


def _normalize_reaction_obj(obj):
    for key in ('reagents', 'products'):
        for compound in obj[key]:
            compound['name'] = compound['name'].strip()
            if 'note' in compound:
                compound['note'] = compound['note'].strip()
    if not isinstance(obj.get('solvent'), str):
        obj['solvent'] = None
    else:
        obj['solvent'] = obj['solvent'].strip()
    if not isinstance(obj.get('catalyst'), str):
        obj['catalyst'] = None
    else:
        obj['catalyst'] = obj['catalyst'].strip()
    return obj


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
        if not _is_valid_reaction_obj(obj):
            _count(stats, 'invalid_schema')
            continue
        reactions.append(_normalize_reaction_obj(obj))
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

    def _get_existing_context(self, cid):
        reactions = self._existing_by_cid.get(cid, [])
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

    def fetch_one(self, chem, preset):
        chem_name = chem['cmpdname']
        existing = self._get_existing_context(chem['cid'])
        fetch_prompt = build_fetch_prompt(
            chem_name, preset.compound_role, preset.scope, existing,
        )

        reactions = None
        chosen_model = None
        for model in self.models:
            response = self.llm_client.fetch_answer_str(fetch_prompt, model)
            stats = {}
            reactions = parse_reactions_jsonl(response, stats)
            self._log_parse_stats(chem_name, f"Fetch response from '{model}'", stats)
            if reactions:
                chosen_model = model
                break

        if not reactions:
            self.logger.log_warn(f"Failed to fetch reactions for '{chem_name}'")
            return {'cid': chem['cid'], 'reactions': None}

        reactions_jsonl_str = '\n'.join(json.dumps(r) for r in reactions)
        revalidate_prompt = build_revalidate_prompt(reactions_jsonl_str)
        revalidated_response = self.llm_client.fetch_answer_str(
            revalidate_prompt, chosen_model,
        )
        stats = {}
        revalidated = parse_reactions_jsonl(revalidated_response, stats)
        self._log_parse_stats(chem_name, "Revalidation response", stats)

        if not revalidated:
            self.logger.log_warn(
                f"Revalidation returned no results for '{chem_name}'. Using original."
            )
            revalidated = reactions

        self.logger.log(
            f"Got {len(revalidated)} reactions for '{chem_name}' with '{chosen_model}'; "
            f"CTT: {self.llm_client.completion_tokens_total}"
        )
        return {'cid': chem['cid'], 'reactions': revalidated}

    def fetch_all(self, preset, max_workers=1):
        raw_fn = self.layout.raw(preset.name)
        processed = {x['cid'] for x in self.store.load_jsonl(raw_fn)}
        criteria = preset.build_criteria(self.compounds)
        staged = [
            chem for chem in self.compounds.chems
            if chem['cid'] not in processed and criteria(chem)
        ]
        self.llm_client.submit_entries_to_llm(
            raw_fn,
            staged,
            max_workers,
            lambda chem: self.fetch_one(chem, preset),
            self.logger,
            description=f"Fetching reactions [{preset.name}]",
        )
