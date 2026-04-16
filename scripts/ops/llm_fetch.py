import re

from scripts.infra.batch_runner import run_batch


class LLMFetchOps:
    def __init__(self, compounds, reactions, store, logger, properties, llm_client, reaction_parser):
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.properties = properties
        self.llm_client = llm_client
        self.reaction_parser = reaction_parser

        self.gpt_oss = reaction_parser.gpt_oss
        self.qwen = reaction_parser.qwen
        self.grok = reaction_parser.grok
        self.gemini = reaction_parser.gemini
        self.deepseek = reaction_parser.deepseek
        self.LLM_HAZARD_CATEGORIES = reaction_parser.LLM_HAZARD_CATEGORIES

    def restrict_ctt(self, max_ctt):
        return self.llm_client.restrict_ctt(max_ctt)

    def fetch_llm_answer_str(self, message, model, reasoning_effort="medium", max_completion_tokens=10000):
        return self.llm_client.fetch_answer_str(
            message,
            model,
            reasoning_effort=reasoning_effort,
            max_completion_tokens=max_completion_tokens,
        )

    def _get_processed_entries(self, out_fn, key):
        entries = self.store.load_jsonl(out_fn)
        if callable(key):
            return {k for x in entries if (k := key(x)) is not None}
        return {x[key] for x in entries}

    def _get_reactions_from_response(self, response):
        return [line for line in response.split('\n') if '->' in line]

    def _submit_entries_to_llm(
        self,
        out_fn,
        entries,
        max_workers,
        routine,
        routine_args=None,
        batch_size=None,
        description="Generation",
    ):
        if routine_args is None:
            routine_args = []
        run_batch(
            self.llm_client,
            entries,
            routine,
            out_fn,
            self.logger,
            max_workers=max_workers,
            routine_args=routine_args,
            batch_size=batch_size,
            description=description,
        )

    def _str_verdict_to_bool(self, verdict):
        verdict = verdict.lower()
        return 'invalid' not in verdict and 'valid' in verdict

    def _extract_verdicts_from_response(self, response):
        return [self._str_verdict_to_bool(verdict) for verdict in response.split('\n')]

    def _get_chem_description(self, chem, valid_cnt):
        instruct = (
            "Write an engaging and informative plain text description of the compound these synonyms refer to. "
            "You may decide freely what aspects to include — composition, properties, uses, history, "
            "relevance, or anything else meaningful — as long as the result is interesting to read and "
            "based on reliable chemical knowledge. "
            "If nothing meaningful can be said about it, respond exactly with the word None.\n"
            "Result should be in plain text format. Markdown or other formatting types are strictly forbidden."
            "Guidelines:\n"
            "- Everything you write must be true or based on well-accepted chemical knowledge.\n"
            "- You may include lesser-known but reliable information, but do not guess or invent facts.\n"
            "- The text should feel lively, colorful, and pleasant to read — avoid dull or purely technical writing.\n"
            "- You may write at any length, as long as the text remains engaging and free of boring or redundant details.\n"
            "- The output should be a smooth, uninterrupted narrative without explicit topic breaks or labeled sections.\n"
        )
        instruct_validate = (
            "Decide if the following text is factually correct and does not contain misleading or speculative "
            "information. If the description is fully correct or at least safely accurate within general chemical "
            "knowledge, return 'Valid'. If it contains anything blatantly wrong, incorrect, or unverified, "
            "return 'Invalid'. Do not write anything except 'Valid' or 'Invalid'."
        )
        instruct_fix = (
            "Review the provided chemical compound description for errors or inconsistencies. Identify any mistakes "
            "and return a corrected, accurate version of the description, preserving its original narrative tone. "
            "Return only the corrected description, without any explanations or extra text."
        )

        max_syns = 3
        chem_name = chem['cmpdname']
        synonyms = ', '.join([f'"{x}"' for x in chem['cmpdsynonym'][:max_syns]])
        models_schedule = [self.gpt_oss, self.gpt_oss, self.qwen, self.qwen]
        av_confidence = 0

        def validate_description(descr, confidence_thr, model):
            confidence = 0
            for valid_i in range(valid_cnt):
                verdict = self.fetch_llm_answer_str(f"{instruct_validate}\n\n{descr}", model)
                if self._str_verdict_to_bool(verdict):
                    confidence += 1 / valid_cnt

                remaining_tries = valid_cnt - valid_i - 1
                max_confidence = remaining_tries / valid_cnt + confidence
                if max_confidence < confidence_thr or confidence >= confidence_thr:
                    confidence += remaining_tries / 2 / valid_cnt
                    break

            return confidence

        for try_i, model in enumerate(models_schedule):
            description = self.fetch_llm_answer_str(f"{synonyms}\n\n{instruct}", model)
            if len(description) < 20:
                self.logger.log_warn(
                    f"Failed to generate description for {chem_name} on try {try_i+1}/{len(models_schedule)} "
                    f"('{model}')"
                )
                continue

            confidence_thr = 0.49
            confidence = validate_description(description, confidence_thr, model)
            if confidence >= confidence_thr:
                self.logger.log(
                    f"('{model}') Generated description for {chem_name} of length {len(description)}. "
                    f"confidence: {confidence}; CTT: {self.llm_client.completion_tokens_total}"
                )
                return {'cid': chem['cid'], 'description': description, 'confidence': confidence, 'source': model}

            description = self.fetch_llm_answer_str(f"{instruct_fix}\n\n{description}", model)
            confidence = validate_description(description, confidence_thr, model)
            if confidence >= confidence_thr:
                self.logger.log(
                    f"('{model}') Generated description for {chem_name} of length {len(description)} after "
                    f"fixing. confidence: {confidence}; CTT: {self.llm_client.completion_tokens_total}"
                )
                return {'cid': chem['cid'], 'description': description, 'confidence': confidence, 'source': model}

            self.logger.log_warn(
                f"Low validation confidence for {chem_name} on try {try_i+1}/{len(models_schedule)} "
                f"('{model}'): {confidence}"
            )
            av_confidence += confidence / len(models_schedule)

        self.logger.log_warn(
            f"Failed to generate description for {chem_name} due to low validation confidence: {av_confidence}"
        )
        return None

    def get_chems_descriptions(self, max_workers=1):
        chem_reactions_occurence = self.properties.get_chems_reactions_occurence()
        processed = self._get_processed_entries(self.reaction_parser.chems_descriptions_fn, 'cid')
        staged_chems = sorted(
            [chem for chem in self.compounds.chems if chem['cid'] not in processed],
            key=lambda chem: chem_reactions_occurence.get(chem['cid'], 0),
            reverse=True,
        )
        self._submit_entries_to_llm(
            self.reaction_parser.chems_descriptions_fn,
            staged_chems,
            max_workers,
            self._get_chem_description,
            routine_args=[6],
            description="Compounds description generation",
        )

    def _get_reactions_description(self, reactions, valid_cnt):
        instruct = (
            "Please generate detailed short plain text descriptions for each of the chemical reactions "
            "provided below. For each reaction, include the general type of reaction, its purpose if "
            "applicable, key transformations, typical reaction conditions, and common solvents or catalysts "
            "if needed. Prefix each description with the corresponding reaction number. Provide only the "
            "descriptions, do not add any extra text or commentary and do not restate or repeat the reactions "
            "themselves. Provide only information that you are confident is correct. Do not use formatting, "
            "but plain text only. If you cannot provide reliable information for a specific reaction, simply "
            "write 'None' as the description."
        )
        instruct_validate = (
            "You will be given pairs of chemical reaction schemes and their descriptions. Each reaction scheme "
            "shows reagents and products (unbalanced). Your task is to validate whether each description "
            "accurately describes the corresponding reaction. Return 'Valid' if the description is chemically "
            "accurate and contains no significant errors. Return 'Invalid' if the description contains factual "
            "errors, incorrect mechanisms, impossible conditions, or misidentified reaction types. Minor "
            "imprecision is acceptable if the core chemistry is correct. Respond with only 'Valid' or 'Invalid' "
            "for each pair, one per line, in the same order as presented."
        )

        models_schedule = [self.gpt_oss, self.qwen]
        for try_i, model in enumerate(models_schedule):
            reactions_formatted_str = '\n'.join([f"{i}. {react['reaction']}" for i, react in enumerate(reactions)])
            messages = [
                {"role": "system", "content": ""},
                {"role": "user", "content": f"{instruct}\n\n{reactions_formatted_str}"},
            ]

            response = self.llm_client._fetch_answer(messages, model)
            descriptions = re.split(r'\n(?=\d+[.)]\s)', response)
            descriptions = [d.strip() for d in descriptions if d.strip()]
            descriptions = [re.sub(r'^\d+[.)]\s+', '', d) for d in descriptions]
            if len(descriptions) != len(reactions):
                self.logger.log_warn(
                    f"Failed to get descriptions ({len(descriptions)} != {len(reactions)}) ('{model}') "
                    f"({try_i+1}/{len(models_schedule)})"
                )
                continue

            indices = [i for i, desc in enumerate(descriptions) if len(desc) > 20]
            descriptions = [descriptions[i] for i in indices]
            reactions = [reactions[i] for i in indices]

            valid_i = 0
            mistakes_cnt = 0
            mistakes_thr = 2
            confidences = [0 for _ in range(len(descriptions))]
            confidence_thr = 0.39
            results = []
            while valid_i < valid_cnt and len(descriptions) > 0:
                formatted = []
                for i in range(len(descriptions)):
                    react = reactions[i]['reaction']
                    formatted.append(f'{i+1}. Scheme: "{react}"\nDescription: "{descriptions[i]}"')
                formatted = '\n\n'.join(formatted)

                response = self.fetch_llm_answer_str(f"{instruct_validate}\n\n{formatted}", model)
                verdicts = self._extract_verdicts_from_response(response)
                if len(verdicts) != len(descriptions):
                    if mistakes_cnt == mistakes_thr:
                        self.logger.log("Returning early due to mistakes at validation phase")
                        return results
                    mistakes_cnt += 1
                    continue

                finished_indices = set()
                remaining_tries = valid_cnt - valid_i - 1
                for i in range(len(verdicts)):
                    confidences[i] += verdicts[i] / valid_cnt
                    max_confidence = remaining_tries / valid_cnt + confidences[i]
                    est_confidence = (confidences[i] + max_confidence) / 2
                    if confidences[i] >= confidence_thr:
                        reaction_str = reactions[i]['reaction']
                        self.logger.log(
                            f"Generated description for '{reaction_str}'; confidence: {est_confidence}; "
                            f"CTT: {self.llm_client.completion_tokens_total}"
                        )

                        entry = {
                            'rid': reactions[i]['rid'],
                            'description': descriptions[i],
                            'confidence': est_confidence,
                            'source': model,
                        }
                        results.append(self.reactions.convert_details_to_canonic(entry))
                        finished_indices.add(i)
                    elif max_confidence < confidence_thr:
                        finished_indices.add(i)

                reactions = [reactions[i] for i in range(len(reactions)) if i not in finished_indices]
                descriptions = [
                    descriptions[i] for i in range(len(descriptions)) if i not in finished_indices
                ]
                confidences = [
                    confidences[i] for i in range(len(confidences)) if i not in finished_indices
                ]
                valid_i += 1

            return results

        self.logger.log_warn("Failed to generate reactions descriptions")
        return None

    def get_reactions_descriptions(self, max_workers=1):
        reactions = self.store.load_jsonl(self.reaction_parser.reactions_parsed_llm_fn)
        chem_reactions_occurence = self.properties.get_chems_reactions_occurence(reactions)

        reactions_power = {}
        for react in reactions:
            all_cids = self.reactions.get_all_reaction_cids(react)
            reactions_power[react['rid']] = sum(chem_reactions_occurence[x] for x in all_cids) / len(all_cids)

        processed = self._get_processed_entries(self.reaction_parser.reactions_details_llm_fn, 'rid')
        reactions = [react for react in reactions if react['rid'] not in processed]
        reactions.sort(key=lambda react: reactions_power[react['rid']], reverse=True)
        reactions = [
            {'reaction': self.reactions.get_reaction_as_str(react), 'rid': react['rid']}
            for react in reactions
        ]

        self._submit_entries_to_llm(
            self.reaction_parser.reactions_details_llm_fn,
            reactions,
            max_workers,
            self._get_reactions_description,
            routine_args=[6],
            batch_size=6,
            description="Reactions descriptions generation",
        )
        self.logger.log(f"Submitted {len(reactions)} reactions for descriptions generation")

    def _fix_unbalanced_reactions_batch(self, reactions):
        model = self.gpt_oss
        instruct = (
            "You will be given a list of unbalanced chemical reaction schemes. Each reaction is generally "
            "valid but contains an error that prevents it from being balanced. Your task is to correct each "
            "reaction so that it can be balanced. Write the corrected reactions, one per line. Do not specify "
            "balance coefficients. If you need to add any compounds to the reaction, use their full chemical "
            "names, not formulas. If a compound is already correctly placed, leave it as it is. Do not write "
            "anything other than the reaction schemes. If a particular reaction cannot be fixed for any reason, "
            "write strictly 'None' on that line."
        )

        reactions_formatted_str = '\n'.join(
            [f"{i+1}. {self.reactions.get_reaction_as_str(react)}" for i, react in enumerate(reactions)]
        )

        def attempt_fix():
            for _ in range(2):
                response = self.fetch_llm_answer_str(f"{instruct}\n\n{reactions_formatted_str}", model)
                response = response.strip().split('\n')
                if len(response) != len(reactions):
                    msg = f": {response[0][:100]}..." if len(response) == 1 else ""
                    self.logger.log_warn(f"Mistake {len(response)} != {len(reactions)}{msg}")
                    continue

                result = []
                for i, react_str in enumerate(response):
                    parsed, _ = self.reaction_parser.parse_reaction_scheme(react_str)
                    if not parsed or not self.reactions.balance_reaction(parsed):
                        continue
                    result.append({'old_rid': reactions[i]['rid'], 'fix': parsed})
                return result
            return []

        def get_intersection(*results):
            if not results:
                return []

            def key(entry):
                return (entry['old_rid'], entry['fix']['rid'])

            intersection_keys = {key(entry) for entry in results[0]}
            for res in results[1:]:
                intersection_keys &= {key(entry) for entry in res}
                if not intersection_keys:
                    return []

            return [entry for entry in results[0] if key(entry) in intersection_keys]

        return get_intersection(*[attempt_fix() for _ in range(2)])

    def fix_unbalanced_reactions(self, max_workers=1):
        processed = self._get_processed_entries(self.reaction_parser.reactions_parsed_fixed_fn, 'old_rid')
        reactions = sorted(list(self.properties.parsed_reactions), key=lambda reaction: reaction['complexity'])
        reactions_staged = [
            react
            for react in reactions
            if not self.reactions.is_react_balanced(react) and react['rid'] not in processed
        ]

        self._submit_entries_to_llm(
            self.reaction_parser.reactions_parsed_fixed_fn,
            reactions_staged,
            max_workers,
            self._fix_unbalanced_reactions_batch,
            batch_size=5,
            description="Unbalanced reactions fix",
        )

    def _get_chems_hazards_categories_batch(self, chems):
        if not chems:
            return None

        categories_str = ', '.join(self.LLM_HAZARD_CATEGORIES)
        instruct = (
            "You will be given a list of chemical compounds. For each compound, assign one or more categories "
            f"from the following list using your common knowledge: {categories_str}.\n"
            "Write the categories exactly as they appear in the list, separated by commas. Each compound's "
            "categories should appear on a separate line in the same order as the input compounds. If a "
            "compound does not fit any category, write 'None' on its corresponding line. Do not include any "
            "additional text or explanation."
        )

        chems_formatted_str = '\n'.join(chem['cmpdname'] for chem in chems)
        model = self.gpt_oss
        runs_num = 6
        occurrence_thr_ratio = 0.49
        mistakes_thr = 2
        mistakes_cnt = 0
        run_i = 0
        categories_count = [dict() for _ in chems]

        while run_i < runs_num:
            response = self.fetch_llm_answer_str(f"{instruct}\n\n{chems_formatted_str}", model)
            response = response.strip().split('\n')
            if len(response) == len(chems):
                valid_run = True
                for i, entry in enumerate(response):
                    cats = [c.lower() for c in re.sub(r'\s+', '', entry).split(',')]
                    if any(c not in self.LLM_HAZARD_CATEGORIES and c != 'none' for c in cats):
                        valid_run = False
                        break

                    for cat in cats:
                        if cat != 'none':
                            categories_count[i][cat] = categories_count[i].get(cat, 0) + 1

                if valid_run:
                    run_i += 1
                    continue

            mistakes_cnt += 1
            if mistakes_cnt >= mistakes_thr:
                self.logger.log_warn("Failed to generate categories after multiple attempts")
                return None

        occurrence_thr = round(runs_num * occurrence_thr_ratio)
        results = []
        for i, chem in enumerate(chems):
            final_cats = [
                cat for cat, count in categories_count[i].items() if count >= occurrence_thr
            ]
            if final_cats:
                self.logger.log(
                    f'Got categories for "{chem["cmpdname"]}": {", ".join(final_cats)}; '
                    f'CTT: {self.llm_client.completion_tokens_total}'
                )
                results.append({'cid': chem['cid'], 'categories': final_cats, 'source': model})

        return results

    def get_chems_hazards_categories(self, max_workers=1):
        processed = self._get_processed_entries(self.reaction_parser.chems_hazard_categories_llm_fn, 'cid')
        chems_staged = [chem for chem in self.compounds.chems if chem['cid'] not in processed]
        self._submit_entries_to_llm(
            self.reaction_parser.chems_hazard_categories_llm_fn,
            chems_staged,
            max_workers,
            self._get_chems_hazards_categories_batch,
            batch_size=10,
            description="Hazard categories generation",
        )

    def _get_chems_nfpa_ratings_batch(self, chems):
        if not chems:
            return None

        instruct = (
            "You will be given a list of chemical compounds. For each compound, assign NFPA ratings based on "
            "your common knowledge. Provide the result in the format: <health>,<flammability>,<stability>, "
            "with three comma-separated integers from 0 to 4. Each compound's ratings should appear on a "
            "separate line, in the same order as the input compounds. Do not include any additional text or explanation."
        )

        chems_formatted_str = '\n'.join(chem['cmpdname'] for chem in chems)
        models_schedule = [self.grok]
        runs_num = 5
        mistakes_thr = 2
        mistakes_cnt = 0
        run_i = 0
        model_i = 0
        ratings = [{'health': 0.0, 'flammability': 0.0, 'instability': 0.0} for _ in chems]

        while run_i < runs_num:
            model = models_schedule[model_i]
            response = self.fetch_llm_answer_str(f"{instruct}\n\n{chems_formatted_str}", model)
            response = response.strip().split('\n')
            if len(response) == len(chems):
                valid_run = True
                for i, entry in enumerate(response):
                    try:
                        entry_values = list(map(float, re.sub(r'\s+', '', entry).split(',')))
                        if len(entry_values) != 3 or any(rate < 0 or rate > 4 for rate in entry_values):
                            valid_run = False
                            break
                        h, f, s = entry_values
                        ratings[i]['health'] += h
                        ratings[i]['flammability'] += f
                        ratings[i]['instability'] += s
                    except Exception:
                        valid_run = False
                        break

                if valid_run:
                    run_i += 1
                    continue

            mistakes_cnt += 1
            if mistakes_cnt >= mistakes_thr:
                model_i += 1
                if model_i == len(models_schedule):
                    self.logger.log_warn("Failed to generate NFPA ratings after multiple attempts")
                    return None
                self.logger.log(f"Falling to the next model: '{models_schedule[model_i]}'")
                mistakes_cnt = 0

        results = []
        for i, chem in enumerate(chems):
            nfpa = {
                'health': ratings[i]['health'] / runs_num,
                'flammability': ratings[i]['flammability'] / runs_num,
                'instability': ratings[i]['instability'] / runs_num,
            }
            results.append({'cid': chem['cid'], 'nfpa': nfpa, 'source': models_schedule[model_i]})
            self.logger.log(
                f'Got NFPA for "{chem["cmpdname"]}": {nfpa}; '
                f'CTT: {self.llm_client.completion_tokens_total}'
            )

        return results

    def get_chems_nfpa_ratings(self, max_workers=1):
        processed = self._get_processed_entries(self.reaction_parser.chems_nfpa_llm_fn, 'cid')
        chems_staged = [chem for chem in self.compounds.chems if chem['cid'] not in processed]
        self._submit_entries_to_llm(
            self.reaction_parser.chems_nfpa_llm_fn,
            chems_staged,
            max_workers,
            self._get_chems_nfpa_ratings_batch,
            batch_size=10,
            description="NFPA ratings generation",
        )
