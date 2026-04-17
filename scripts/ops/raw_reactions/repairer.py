from .schema import is_valid_reaction_obj
from .prompts import REPAIR_BATCH_INSTRUCT, format_repair_candidates
from .validator import ACCEPT_THR, is_confident
from .services import parse_repair_response
from scripts.infra.batch_runner import run_batch
from scripts.infra.fallback import with_fallback


class RawReactionsRepairer:
    def __init__(self, llm_client, store, logger, layout, parser, validator, models):
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.parser = parser
        self.validator = validator
        self.models = models

    def _get_reaction_id(self, reaction_obj):
        if not is_valid_reaction_obj(reaction_obj):
            return None, None, 'invalid_schema'
        parsed, _ = self.parser.parse_structured_reaction(reaction_obj)
        if not parsed:
            return None, None, 'unmapped'
        return parsed['rid'], parsed, None

    def _base_result(self, group):
        return {
            'rid': group['rid'],
            'candidate_count': group['candidate_count'],
            'original_confidence': group['original_confidence'],
        }

    def _failure_result(self, group, reason, repair_source=None):
        result = self._base_result(group)
        result['fixed_reaction'] = None
        result['reason'] = reason
        if repair_source is not None:
            result['repair_source'] = repair_source
        return result

    def _success_result(self, group, fixed_reaction, fixed_rid, repair_source, validation):
        result = self._base_result(group)
        result.update({
            'fixed_reaction': fixed_reaction,
            'fixed_rid': fixed_rid,
            'repair_source': repair_source,
            'validation_source': validation.get('source'),
            'validation_confidence': validation.get('confidence', 0.0),
            'rounds': validation.get('rounds', 0),
            'positives': validation.get('positives', 0),
        })
        return result

    def _choose_representative(self, entries):
        return max(
            entries,
            key=lambda entry: (
                entry.get('confidence', 0.0),
                entry.get('positives', 0),
                entry.get('rounds', 0),
            ),
        )

    def _build_candidate_groups(self, lower_bound, acceptance_threshold):
        verdict_fn = self.layout.verdict()
        repaired_rids = {
            entry.get('rid')
            for entry in self.store.load_jsonl(self.layout.repairs())
            if entry.get('rid') is not None
        }
        groups = {}
        skip_stats = {
            'invalid_schema': 0,
            'unmapped': 0,
            'already_repaired': 0,
        }

        for entry in self.store.load_jsonl(verdict_fn):
            rid, _, reason = self._get_reaction_id(entry.get('reaction'))
            if rid is None:
                skip_stats[reason] += 1
                continue
            if rid in repaired_rids:
                skip_stats['already_repaired'] += 1
                continue

            group = groups.setdefault(rid, {
                'rid': rid,
                'entries': [],
                'accepted': False,
                'borderline': False,
            })
            confidence = entry.get('confidence', 0.0)
            group['accepted'] = group['accepted'] or confidence >= acceptance_threshold
            group['borderline'] = group['borderline'] or (
                lower_bound <= confidence < acceptance_threshold
            )
            group['entries'].append(entry)

        candidates = []
        for group in groups.values():
            if group['accepted'] or not group['borderline']:
                continue
            representative = self._choose_representative(group['entries'])
            candidates.append({
                'rid': group['rid'],
                'reaction': representative['reaction'],
                'cid': representative.get('cid'),
                'candidate_count': len(group['entries']),
                'original_confidence': representative.get('confidence', 0.0),
            })

        if any(skip_stats.values()):
            self.logger.log_warn(f"Repair skipped reactions: {skip_stats}")
        return candidates

    def _parse_repair_response(self, response, expected_count):
        return parse_repair_response(response, expected_count)

    def _get_fix(self, candidates, model, max_retries=2):
        """Attempt to get repair suggestions from a single model.

        Returns (repairs, model) on success, or (None, None) on persistent
        format errors (signal to caller to try another model).
        """
        prompt = f"{REPAIR_BATCH_INSTRUCT}\n\n{format_repair_candidates(candidates)}"
        for _ in range(max_retries):
            response = self.llm_client.fetch_answer_str(prompt, model)
            repairs = self._parse_repair_response(response, len(candidates))
            if repairs is not None:
                return repairs, model
        return None, None

    def repair_batch(
        self,
        candidates,
        repair_model,
        validation_model,
        acceptance_threshold=ACCEPT_THR,
    ):
        """Run one repair pass with explicit repair and validation models.

        repair_model:      single model used to generate fixes.
        validation_model:  single model used to validate generated fixes.

        Returns a list of result dicts, or None if repair parsing failed
        or validation parsing failed (signal to caller to try another model).
        """
        repairs, repair_source = self._get_fix(candidates, repair_model)
        if repairs is None:
            return None

        results = []
        validation_entries = []
        validation_context = []
        for group, fixed_reaction in zip(candidates, repairs):
            if fixed_reaction is None:
                results.append(self._failure_result(group, 'null_repair', repair_source))
                continue
            if not is_valid_reaction_obj(fixed_reaction):
                results.append(self._failure_result(group, 'invalid_schema', repair_source))
                continue
            fixed_rid, _, reason = self._get_reaction_id(fixed_reaction)
            if fixed_rid is None:
                results.append(self._failure_result(group, reason, repair_source))
                continue
            validation_entries.append({'cid': group.get('cid'), 'reaction': fixed_reaction})
            validation_context.append((group, fixed_reaction, fixed_rid))

        if validation_entries:
            validations = self.validator.validate_batch(
                validation_entries,
                model=validation_model,
                max_rounds=9,
                acceptance_threshold=acceptance_threshold,
            )

            if validations is None:
                return None

            for validation, (group, fixed_reaction, fixed_rid) in zip(
                validations,
                validation_context,
            ):
                if is_confident(validation, acceptance_threshold=acceptance_threshold):
                    results.append(
                        self._success_result(
                            group,
                            fixed_reaction,
                            fixed_rid,
                            repair_source,
                            validation,
                        )
                    )
                else:
                    results.append(
                        self._failure_result(group, 'validation_rejected', repair_source)
                    )

        return results

    def repair(
        self,
        max_workers=1,
        lower_bound=0.1,
        acceptance_threshold=ACCEPT_THR,
        repair_models=None,
        validation_models=None,
    ):
        repair_models = repair_models or self.models
        validation_models = validation_models or self.validator.models

        candidates = self._build_candidate_groups(
            lower_bound=lower_bound,
            acceptance_threshold=acceptance_threshold,
        )

        def _batch_with_fallback(batch, acceptance_threshold):
            for repair_model in repair_models:
                results = with_fallback(
                    lambda validation_model: self.repair_batch(
                        batch,
                        repair_model=repair_model,
                        validation_model=validation_model,
                        acceptance_threshold=acceptance_threshold,
                    ),
                    validation_models,
                    logger=self.logger,
                )
                if results is not None:
                    return results
            return []

        run_batch(
            self.llm_client,
            candidates,
            _batch_with_fallback,
            self.layout.repairs(),
            self.logger,
            max_workers=max_workers,
            routine_args=[acceptance_threshold],
            batch_size=5,
            description="Repairing reactions",
        )
