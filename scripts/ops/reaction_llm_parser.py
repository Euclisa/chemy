import os


class ReactionLLMParser:
    def __init__(self, data_dir, compounds, reactions, store, logger):
        self.data_dir = data_dir
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger

        self.gpt_oss = "openai/gpt-oss-120b"
        self.qwen = "qwen/qwen3-235b-a22b"
        self.grok = "x-ai/grok-4-fast"
        self.gemini = "google/gemini-2.5-flash-lite"
        self.deepseek = "deepseek/deepseek-v3.2-exp"

        self.reactions_parsed_fixed_fn = os.path.join(
            self.reactions.parsed_reactions_dir,
            'reactions_parsed_fixed.jsonl',
        )
        self.reactions_parsed_llm_fn = os.path.join(self.reactions.parsed_reactions_dir, 'llm')
        os.makedirs(self.reactions_parsed_llm_fn, exist_ok=True)

        self.reactions_details_llm_fn = os.path.join(self.reactions.reactions_details_dir, 'llm')
        os.makedirs(self.reactions_details_llm_fn, exist_ok=True)

        self.reactions_descriptions_fn = os.path.join(
            self.data_dir,
            'reactions_details',
            'reactions_descriptions.jsonl',
        )

        self.chems_properties_llm_dir = os.path.join(self.compounds.chems_properties_dir, 'llm')
        os.makedirs(self.chems_properties_llm_dir, exist_ok=True)
        self.chems_descriptions_fn = os.path.join(self.chems_properties_llm_dir, 'chems_descriptions.jsonl')
        self.chems_hazard_categories_llm_fn = os.path.join(
            self.chems_properties_llm_dir,
            'chems_hazard_categories_llm.jsonl',
        )
        self.chems_nfpa_llm_fn = os.path.join(
            self.chems_properties_llm_dir,
            'chems_nfpa_llm.jsonl',
        )

        store.register_sorting(self.reactions_parsed_llm_fn, 'rid')
        store.register_sorting(self.reactions_details_llm_fn, 'rid')
        store.register_vault(self.reactions_parsed_llm_fn, 'reactions_')
        store.register_vault(self.reactions_details_llm_fn, 'reactions_details_')

        self.LLM_HAZARD_CATEGORIES = {
            'explosive',
            'acute_toxic',
            'flammable',
            'oxidizer',
            'corrosive',
            'serious_health_hazard',
            'environment_hazard',
        }

    def parse_reaction_scheme(self, reaction_str: str):
        parts = reaction_str.split('->')
        if len(parts) != 2:
            return None, set()

        reagents_str, products_str = parts
        parse_success = True
        unmapped_names = set()

        def parse_compounds(compound_str, skip_names, existing_cids=None):
            nonlocal parse_success

            existing_cids = set() if existing_cids is None else existing_cids
            compounds = []
            cids = set()
            for name in compound_str.split('+'):
                norm = self.compounds.normalize_chem_name(name)
                if norm in skip_names:
                    continue

                clean = self.compounds.clean_chem_name(name)
                cid = self.compounds.name_cid_map.get(norm)
                if cid is None:
                    parse_success = False
                    unmapped_names.add((norm, clean))

                if cid is None or cid not in existing_cids | cids:
                    compounds.append({'norm_name': norm, 'original_name': clean, 'cid': cid})
                    cids.add(cid)

            compounds = list({c["cid"]: c for c in compounds}.values())
            if parse_success:
                compounds.sort(key=lambda c: c['cid'])

            return compounds, cids

        reagents, reagents_cids = parse_compounds(reagents_str, {"light", "heat", "catalyst"})
        products, products_cids = parse_compounds(products_str, {"otherproducts"}, reagents_cids)

        if products_cids & reagents_cids or not products or not reagents:
            parse_success = False

        if not parse_success:
            return None, unmapped_names

        reaction = self.reactions.assemble_reaction({'reagents': reagents, 'products': products})
        return reaction, unmapped_names

