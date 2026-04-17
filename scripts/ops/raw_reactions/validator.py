from .prompts import VALIDATE_BATCH_INSTRUCT, format_reactions_for_validation
from .schema import is_valid_reaction_obj
from .services import RawReactionValidationService, parse_indexed_verdicts
from scripts.infra.batch_runner import run_batch
from scripts.infra.fallback import with_fallback


ACCEPT_THR = 0.4


def is_confident(entry, key='confidence', acceptance_threshold=ACCEPT_THR):
    confidence = entry.get(key, 0.0)
    return confidence is not None and confidence >= acceptance_threshold


def get_reaction_id(parser, reaction_obj):
    if not is_valid_reaction_obj(reaction_obj):
        return None, 'invalid_schema'
    parsed, _ = parser.parse_structured_reaction(reaction_obj)
    if not parsed:
        return None, 'unmapped'
    return parsed['rid'], None


class RawReactionsValidator:
    def __init__(self, llm_client, store, logger, layout, parser, models):
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.parser = parser
        self.models = models
        self.service = RawReactionValidationService(llm_client, parser)

    def _str_verdict_to_bool(self, verdict):
        verdict = verdict.lower()
        return 'invalid' not in verdict and 'valid' in verdict

    def _extract_verdicts(self, response):
        verdicts = parse_indexed_verdicts(response, len(response.split('\n')))
        if verdicts is not None:
            return verdicts
        return [self._str_verdict_to_bool(v) for v in response.split('\n')]

    def _get_verdicts_safe(self, batch, model, max_retries=3):
        for _ in range(max_retries):
            prompt = f"{VALIDATE_BATCH_INSTRUCT}\n{format_reactions_for_validation(batch)}"
            response = self.llm_client.fetch_answer_str(prompt, model)
            verdicts = self._extract_verdicts(response)
            if len(verdicts) == len(batch):
                return verdicts
        return None

    def _get_reaction_id(self, reaction_obj):
        return get_reaction_id(self.parser, reaction_obj)

    def _build_result(self, reaction_entry, positives, rounds_done, model):
        confidence = positives / rounds_done if rounds_done > 0 else 0.0
        result = reaction_entry.copy()
        result['confidence'] = confidence
        result['rounds'] = rounds_done
        result['positives'] = positives
        result['source'] = model
        return result

    def _as_verdict_result(self, result):
        result = result.copy()
        result.pop('rid', None)
        return result

    def validate_batch(
        self,
        reactions,
        model,
        max_rounds=9,
        acceptance_threshold=ACCEPT_THR,
        early_stop=True,
    ):
        """Validate a batch of reactions with a single model.

        Returns a list of result dicts on success, or None if the model
        produces persistent format errors (signal to caller to try another model).
        """
        results = self.service.validate_entries(
            reactions,
            model=model,
            max_rounds=max_rounds,
            acceptance_threshold=acceptance_threshold,
            early_stop=early_stop,
        )
        if results is None:
            self.logger.log_warn(
                f"Falling to another model due to format errors ('{model}')"
            )
            return None
        for result in results:
            self.logger.log(
                f"Validated reaction; "
                f"confidence: {result['confidence']:.2f} "
                f"({result['positives']}/{result['rounds']}); "
                f"CTT: {self.llm_client.completion_tokens_total}"
            )
        return [self._as_verdict_result(result) for result in results]

    def _log_skip_stats(self, preset_name, stats):
        if any(stats.values()):
            self.logger.log_warn(
                f"Validation skipped reactions [{preset_name}]: {dict(stats)}"
            )

    def validate(self, preset, max_workers=1, acceptance_threshold=ACCEPT_THR):
        verdict_fn = self.layout.verdict()

        skip_stats = {
            'invalid_existing_verdict_schema': 0,
            'unmapped_existing_verdict': 0,
            'invalid_raw_schema': 0,
            'unmapped_raw': 0,
            'duplicate_rid': 0,
            'empty_raw_entry': 0,
        }

        processed_rids = set()
        for entry in self.store.load_jsonl(verdict_fn):
            rid, reason = self._get_reaction_id(entry.get('reaction'))
            if rid is not None:
                processed_rids.add(rid)
            elif reason == 'invalid_schema':
                skip_stats['invalid_existing_verdict_schema'] += 1
            elif reason == 'unmapped':
                skip_stats['unmapped_existing_verdict'] += 1

        reactions = []
        for raw_fn in self.layout.raw_all(preset.name):
            for entry in self.store.load_jsonl(raw_fn):
                raw_reactions = entry.get('reactions') or []
                if not raw_reactions:
                    skip_stats['empty_raw_entry'] += 1
                for reaction_obj in raw_reactions:
                    rid, reason = self._get_reaction_id(reaction_obj)
                    if rid is None:
                        if reason == 'invalid_schema':
                            skip_stats['invalid_raw_schema'] += 1
                        elif reason == 'unmapped':
                            skip_stats['unmapped_raw'] += 1
                        continue
                    if rid in processed_rids:
                        skip_stats['duplicate_rid'] += 1
                        continue
                    reactions.append({'cid': entry['cid'], 'reaction': reaction_obj})
                    processed_rids.add(rid)

        self._log_skip_stats(preset.name, skip_stats)

        def _batch_with_fallback(batch, max_rounds, acceptance_threshold):
            return with_fallback(
                lambda model: self.validate_batch(
                    batch, model, max_rounds, acceptance_threshold,
                ),
                self.models,
                logger=self.logger,
            )

        run_batch(
            self.llm_client,
            reactions,
            _batch_with_fallback,
            verdict_fn,
            self.logger,
            max_workers=max_workers,
            routine_args=[9, acceptance_threshold],
            batch_size=5,
            description=f"Validating reactions [{preset.name}]",
            store=self.store,
        )
