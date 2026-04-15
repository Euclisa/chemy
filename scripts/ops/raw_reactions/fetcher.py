from .prompts import build_fetch_prompt, build_revalidate_prompt


class RawReactionsFetcher:
    def __init__(self, compounds, llm_client, store, logger, layout, models):
        self.compounds = compounds
        self.llm_client = llm_client
        self.store = store
        self.logger = logger
        self.layout = layout
        self.models = models  # ordered fallback list

    def _get_reactions_from_response(self, response):
        return [line for line in response.split('\n') if '->' in line]

    def fetch_one(self, chem: dict, preset) -> dict:
        chem_name = chem['cmpdname']
        fetch_prompt = build_fetch_prompt(chem_name, preset.compound_role, preset.scope)
        revalidate_prompt = build_revalidate_prompt()

        reactions = None
        chosen_model = None
        for model in self.models:
            response = self.llm_client.fetch_answer_str(fetch_prompt, model)
            reactions = self._get_reactions_from_response(response)
            if reactions:
                chosen_model = model
                break

        if not reactions:
            self.logger.log_warn(f"Failed to fetch reactions for '{chem_name}'")
            return {'cid': chem['cid'], 'reactions': None}

        messages = [
            {"role": "system", "content": ""},
            {"role": "user", "content": fetch_prompt},
            {"role": "assistant", "content": '\n'.join(reactions)},
            {"role": "user", "content": revalidate_prompt},
        ]
        response = self.llm_client._fetch_answer(messages, chosen_model)
        reactions_revalid = self._get_reactions_from_response(response)
        if not reactions_revalid:
            self.logger.log_warn(
                f"Failed to fetch revalidated reactions for '{chem_name}'. Assuming all are valid..."
            )
            reactions_revalid = reactions

        self.logger.log(
            f"Got {len(reactions_revalid)} reactions for '{chem_name}' with '{chosen_model}'; "
            f"CTT: {self.llm_client.completion_tokens_total}"
        )
        return {'cid': chem['cid'], 'reactions': reactions_revalid}

    def fetch_all(self, preset, max_workers: int = 1):
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
