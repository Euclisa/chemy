from .prompts import VALIDATE_BATCH_INSTRUCT, format_reactions_for_validation
from .fetcher import _is_valid_reaction_obj


ACCEPT_THR = 0.5


class RawReactionsValidator:
    def __init__(self, llm_client, store, logger, layout, parser, models):
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.parser = parser
        self.models = models

    def _str_verdict_to_bool(self, verdict):
        verdict = verdict.lower()
        return 'invalid' not in verdict and 'valid' in verdict

    def _extract_verdicts(self, response):
        return [self._str_verdict_to_bool(v) for v in response.split('\n')]

    def _get_verdicts_safe(self, batch, model, max_retries=3):
        for _ in range(max_retries):
            prompt = f"{VALIDATE_BATCH_INSTRUCT}\n{format_reactions_for_validation(batch)}"
            response = self.llm_client.fetch_answer_str(prompt, model)
            verdicts = self._extract_verdicts(response)
            if len(verdicts) == len(batch):
                return verdicts
        return None

    def _build_result(self, reaction, positives, rounds_done, model):
        confidence = positives / rounds_done if rounds_done > 0 else 0.0
        result = reaction.copy()
        result['confidence'] = confidence
        result['valid'] = confidence >= ACCEPT_THR
        result['rounds'] = rounds_done
        result['positives'] = positives
        result['source'] = model
        return result

    def _validate_batch(self, reactions, max_rounds):
        for model in self.models:
            active = list(range(len(reactions)))
            positives = [0] * len(reactions)
            rounds_done = [0] * len(reactions)
            results = []
            bad = False

            for round_i in range(max_rounds):
                if not active:
                    break

                batch = [reactions[i] for i in active]
                verdicts = self._get_verdicts_safe(batch, model)

                if verdicts is None:
                    self.logger.log_warn(
                        f"Falling to another model due to format errors ('{model}')"
                    )
                    bad = True
                    break

                remaining = max_rounds - round_i - 1
                next_active = []
                for j, idx in enumerate(active):
                    rounds_done[idx] += 1
                    positives[idx] += int(verdicts[j])

                    best_possible = (positives[idx] + remaining) / max_rounds
                    worst_possible = positives[idx] / max_rounds

                    if best_possible < ACCEPT_THR or worst_possible >= ACCEPT_THR:
                        result = self._build_result(
                            reactions[idx], positives[idx], rounds_done[idx], model,
                        )
                        results.append(result)
                        self.logger.log(
                            f"Validated reaction; "
                            f"confidence: {result['confidence']:.2f} "
                            f"({positives[idx]}/{rounds_done[idx]}); "
                            f"CTT: {self.llm_client.completion_tokens_total}"
                        )
                    else:
                        next_active.append(idx)

                active = next_active

            if bad:
                continue

            for idx in active:
                result = self._build_result(
                    reactions[idx], positives[idx], rounds_done[idx], model,
                )
                results.append(result)

            return results

        return None

    def _get_reaction_id(self, reaction_obj):
        if not _is_valid_reaction_obj(reaction_obj):
            return None, 'invalid_schema'
        parsed, _ = self.parser.parse_structured_reaction(reaction_obj)
        if not parsed:
            return None, 'unmapped'
        return parsed['rid'], None

    def _log_skip_stats(self, preset_name, stats):
        if any(stats.values()):
            self.logger.log_warn(
                f"Validation skipped reactions [{preset_name}]: {dict(stats)}"
            )

    def validate(self, preset, max_workers=1):
        raw_fn = self.layout.raw(preset.name)
        verdict_fn = self.layout.verdict(preset.name)

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
        self.llm_client.submit_entries_to_llm(
            verdict_fn,
            reactions,
            max_workers,
            self._validate_batch,
            self.logger,
            routine_args=[9],
            batch_size=5,
            description=f"Validating reactions [{preset.name}]",
        )
