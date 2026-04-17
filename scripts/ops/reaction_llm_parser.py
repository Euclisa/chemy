import os

from .raw_reactions.models import DEEPSEEK, GEMINI, GPT_OSS, GROK, QWEN


class ReactionLLMParser:
    def __init__(self, data_dir, compounds, reactions, store, logger):
        self.data_dir = data_dir
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger

        self.gpt_oss = GPT_OSS
        self.qwen = QWEN
        self.grok = GROK
        self.gemini = GEMINI
        self.deepseek = DEEPSEEK

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

    def parse_structured_reaction(self, reaction_dict):
        if not isinstance(reaction_dict, dict):
            return None, set()

        parse_success = True
        unmapped_names = set()

        def parse_compounds(compounds, skip_names, existing_cids=None):
            nonlocal parse_success
            existing_cids = set() if existing_cids is None else existing_cids
            result = []
            cids = set()

            for compound in compounds:
                name = compound.get('name', '')
                phase = compound.get('phase')
                note = compound.get('note')
                norm = self.compounds.normalize_chem_name(name)
                if norm in skip_names:
                    continue

                cid = self.compounds.unique_name_cid_map.get(norm)
                if cid is None:
                    parse_success = False
                    unmapped_names.add((norm, name))

                if cid is None or cid not in existing_cids | cids:
                    entry = {'norm_name': norm, 'original_name': name, 'cid': cid}
                    if phase:
                        entry['phase'] = phase
                    if note:
                        entry['note'] = note
                    result.append(entry)
                    cids.add(cid)

            result = list({c['cid']: c for c in result}.values())
            if parse_success:
                result.sort(key=lambda c: c['cid'])

            return result, cids

        reagents, reagents_cids = parse_compounds(
            reaction_dict.get('reagents', []),
            {"light", "heat", "catalyst"},
        )
        products, products_cids = parse_compounds(
            reaction_dict.get('products', []),
            {"otherproducts"},
            reagents_cids,
        )

        if products_cids & reagents_cids or not products or not reagents:
            parse_success = False

        if not parse_success:
            return None, unmapped_names

        reaction = self.reactions.assemble_reaction({
            'reagents': reagents,
            'products': products,
        })

        if reaction_dict.get('solvent'):
            reaction['solvent'] = reaction_dict['solvent']
        if reaction_dict.get('catalyst'):
            reaction['catalyst'] = reaction_dict['catalyst']

        return reaction, unmapped_names
