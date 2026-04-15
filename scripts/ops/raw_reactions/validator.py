from .prompts import VALIDATE_BATCH_INSTRUCT, format_reactions_for_validation


class RawReactionsValidator:
    def __init__(self, llm_client, store, logger, layout, parser, models):
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.parser = parser  # ReactionLLMParser — provides parse_reaction_scheme
        self.models = models  # ordered fallback list

    def _str_verdict_to_bool(self, verdict):
        verdict = verdict.lower()
        return 'invalid' not in verdict and 'valid' in verdict

    def _extract_verdicts(self, response):
        return [self._str_verdict_to_bool(v) for v in response.split('\n')]

    def _validate_batch(self, reactions, valid_cnt):
        for model in self.models:
            valid_i = 0
            mistakes_cnt = 0
            mistakes_thr = 3
            confidences = [0.0] * len(reactions)
            confidence_thr = 0.5
            bad = False
            results = []
            staged = list(reactions)

            while valid_i < valid_cnt and staged:
                prompt = f"{VALIDATE_BATCH_INSTRUCT}\n{format_reactions_for_validation(staged)}"
                response = self.llm_client.fetch_answer_str(prompt, model)
                verdicts = self._extract_verdicts(response)

                if len(verdicts) != len(staged):
                    if mistakes_cnt == mistakes_thr:
                        self.logger.log_warn(
                            f"Falling to another model due to mistakes ('{model}')"
                        )
                        bad = True
                        break
                    mistakes_cnt += 1
                    continue

                finished_indices = set()
                remaining_tries = valid_cnt - valid_i - 1
                for i, verdict in enumerate(verdicts):
                    confidences[i] += verdict / valid_cnt
                    max_confidence = remaining_tries / valid_cnt + confidences[i]
                    est_confidence = (confidences[i] + max_confidence) / 2
                    if confidences[i] >= confidence_thr or max_confidence <= confidence_thr:
                        react = staged[i].copy()
                        react['valid'] = confidences[i] > confidence_thr
                        react['confidence'] = est_confidence
                        react['source'] = model
                        results.append(react)
                        finished_indices.add(i)
                        self.logger.log(
                            f"Processed reaction '{staged[i]['reaction']}'; "
                            f"confidence: {est_confidence:.2f}; "
                            f"CTT: {self.llm_client.completion_tokens_total}"
                        )

                staged = [staged[i] for i in range(len(staged)) if i not in finished_indices]
                confidences = [confidences[i] for i in range(len(confidences)) if i not in finished_indices]
                valid_i += 1

            if not bad:
                return results

        return None

    def _get_reaction_id(self, react):
        react_parsed, _ = self.parser.parse_reaction_scheme(react['reaction'])
        if not react_parsed:
            self.logger.log_warn(f"Verdict file contains non-parsable reaction: '{react}'")
            return None
        return react_parsed['rid']

    def validate(self, preset, max_workers: int = 1):
        raw_fn = self.layout.raw(preset.name)
        verdict_fn = self.layout.verdict(preset.name)

        processed_rids = {
            rid for entry in self.store.load_jsonl(verdict_fn)
            if (rid := self._get_reaction_id(entry)) is not None
        }

        reactions = []
        for entry in self.store.load_jsonl(raw_fn):
            for react in (entry['reactions'] or []):
                react_parsed, _ = self.parser.parse_reaction_scheme(react)
                if not react_parsed or react_parsed['rid'] in processed_rids:
                    continue
                reactions.append({'cid': entry['cid'], 'reaction': react})
                processed_rids.add(react_parsed['rid'])

        self.llm_client.submit_entries_to_llm(
            verdict_fn,
            reactions,
            max_workers,
            self._validate_batch,
            self.logger,
            routine_args=[9],
            batch_size=10,
            description=f"Validating reactions [{preset.name}]",
        )
