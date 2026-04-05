import os
import numpy as np

from functools import cached_property


class PropertiesCompiler:
    LLM_HAZARD_CATEGORIES = {
        'explosive',
        'acute_toxic',
        'flammable',
        'oxidizer',
        'corrosive',
        'serious_health_hazard',
        'environment_hazard',
    }

    def __init__(self, data_dir, compounds, reactions, store, logger, crc, solubility, hazards):
        self.data_dir = data_dir
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.crc = crc
        self.solubility = solubility
        self.hazards = hazards

        self.reactions_parsed_dir = os.path.join(self.data_dir, 'reactions_parsed')
        self.reactions_details_dir = os.path.join(self.data_dir, 'reactions_details')

        self.reactions_parsed_ord_fn = os.path.join(self.reactions_parsed_dir, 'ord')
        self.reactions_parsed_llm_fn = os.path.join(self.reactions_parsed_dir, 'llm')
        self.reactions_details_ord_fn = os.path.join(self.reactions_details_dir, 'ord')
        self.reactions_details_llm_fn = os.path.join(self.reactions_details_dir, 'llm')

        self.chems_properties_compiled_fn = os.path.join(
            self.compounds.chems_properties_dir,
            'chems_properties.jsonl',
        )
        self.chems_curiosity_fn = os.path.join(self.compounds.chems_properties_dir, 'chems_curiosity.jsonl')

        store.register_sorting(self.chems_properties_compiled_fn, 'cid')
        store.register_sorting(self.chems_curiosity_fn, ('curiosity', True))

        self.property_max_synonyms = 7
        self.property_max_cas = 5
        self.sources_priority = {
            'ord': 10,
            'openai/gpt-oss-120b': 5,
            'qwen/qwen3-235b-a22b': 4,
            'x-ai/grok-4-fast': 3,
        }

    def clear_cached_property(self, name):
        cls_attr = getattr(type(self), name, None)
        if isinstance(cls_attr, cached_property) and name in self.__dict__:
            del self.__dict__[name]

    def _clear_parsed_caches(self):
        self.clear_cached_property('parsed_reactions')
        self.clear_cached_property('parsed_reactions_balanced')
        self.clear_cached_property('reactions_details')

    @cached_property
    def parsed_reactions(self):
        ord_reactions = self.store.load_jsonl(self.reactions_parsed_ord_fn)
        llm_reactions = self.store.load_jsonl(self.reactions_parsed_llm_fn)
        parsed_reactions = {react['rid']: react for react in llm_reactions + ord_reactions}
        return tuple(parsed_reactions.values())

    @cached_property
    def parsed_reactions_balanced(self):
        return tuple(react for react in self.parsed_reactions if react['rid'] in self.reactions.reactions_balance)

    @cached_property
    def reactions_details(self):
        ord_details = self.store.load_jsonl(self.reactions_details_ord_fn)
        llm_details = self.store.load_jsonl(self.reactions_details_llm_fn)
        parsed_reactions = {entry['rid']: entry for entry in llm_details + ord_details}
        return tuple(parsed_reactions.values())

    def write_parsed_reactions(self, reactions):
        ord_reactions = []
        llm_reactions = []
        for react in reactions:
            if react['source'] == 'ord':
                ord_reactions.append(react)
            else:
                llm_reactions.append(react)

        self.store.write_jsonl(ord_reactions, self.reactions_parsed_ord_fn)
        self.store.write_jsonl(llm_reactions, self.reactions_parsed_llm_fn)
        self._clear_parsed_caches()

    def merge_parsed_files(self, out_fn, *parsed_reactions_files):
        rid_reaction = dict()
        for fn in parsed_reactions_files:
            reactions = self.store.load_jsonl(fn)
            for react in reactions:
                rid = react['rid']
                if rid not in rid_reaction:
                    rid_reaction[rid] = react
                    continue

                old_react = rid_reaction[rid]
                old_source = old_react.get('source')
                new_source = react.get('source')
                new_source_priority = self.sources_priority.get(new_source, -1)
                old_source_priority = self.sources_priority.get(old_source, -1)
                if new_source_priority > old_source_priority:
                    rid_reaction[rid] = react

        self.store.write_jsonl(list(rid_reaction.values()), out_fn, backup=False)

    def generate_edges(self):
        edge_reaction_id_map = dict()
        for react in self.logger.track(self.parsed_reactions, "Generating edges..."):
            react_id = react['rid']
            for reagent in react['reagents']:
                for product in react['products']:
                    edge = (reagent['cid'], product['cid'])
                    edge_reaction_id_map.setdefault(edge, []).append(react_id)

        def get_eid(edge):
            source = edge[0]
            target = edge[1]
            return ((source + target) * target) % 2**64

        edges = []
        for edge, reaction_ids in edge_reaction_id_map.items():
            edges.append(
                {
                    'eid': get_eid(edge),
                    'source': edge[0],
                    'target': edge[1],
                    'reactions': reaction_ids,
                }
            )

        self.store.write_jsonl(edges, self.compounds.chems_edges_fn, backup=False)
        self.logger.log(f"Generated {len(edge_reaction_id_map)} edges")

    def discard_orphaned_unmapped_chems(self):
        connected_cids = set()
        for reaction in self.parsed_reactions:
            connected_cids.update(self.reactions.get_all_reaction_cids(reaction))

        for entry in self.reactions_details:
            for cmpd in entry['solvents'] + entry['catalysts']:
                connected_cids.add(cmpd['cid'])

        res_chems = [chem for chem in self.compounds.chems_unmapped if chem['cid'] in connected_cids]
        self.compounds.update_unmapped_chems(res_chems)

    def balance_parsed_reactions(self, reactions_parsed_fn=None):
        if reactions_parsed_fn is None:
            reactions = self.store.load_jsonl(self.reactions_parsed_llm_fn)
            reactions += self.store.load_jsonl(self.reactions_parsed_ord_fn)
        else:
            reactions = self.store.load_jsonl(reactions_parsed_fn)

        balances = self.store.load_jsonl(self.reactions.reactions_balance_fn)
        balanced_rids = set(entry['rid'] for entry in balances)

        balanced_cnt = 0
        for react in self.logger.track(reactions, "Balancing parsed reactions"):
            if react['rid'] in balanced_rids:
                continue

            bal = self.reactions.balance_reaction(react)
            if bal:
                balances.append(bal)
                balanced_cnt += 1

        self.store.write_jsonl(balances, self.reactions.reactions_balance_fn)
        self.reactions.clear_cached_property('reactions_balance')
        self.clear_cached_property('parsed_reactions_balanced')
        self.logger.log(
            f"Processed {len(reactions)}; balanced {balanced_cnt}; total balanced: {len(balances)}"
        )

    def get_chems_reactions_occurence(self, reactions=None):
        reactions = self.parsed_reactions if reactions is None else reactions

        chem_reactions_occurence = {chem['cid']: 0 for chem in self.compounds.chems}
        for react in reactions:
            for cid in self.reactions.get_all_reaction_cids(react):
                chem_reactions_occurence[cid] += 1

        return chem_reactions_occurence

    def compile_chems_properties(self):
        solubility = self.store.load_jsonl(self.solubility.chems_solubility_fn)
        crc_inorganic = self.store.load_jsonl(self.crc.crc_inorganic_constants_fn)
        crc_organic = self.store.load_jsonl(self.crc.crc_organic_constants_fn)
        crc_flammability = self.store.load_jsonl(self.crc.crc_flammability_fn)

        result = dict()

        def add_property(cid, prop, value):
            if cid in result:
                result[cid].append({'property': prop, 'value': f"{value}"})

        for chem in self.logger.track(self.compounds.chems, 'Adding basic properties', transient=True):
            cid = chem['cid']
            result[cid] = []

            add_property(cid, 'PubChem CID', cid)
            add_property(cid, 'Synonyms', ', '.join(chem['cmpdsynonym'][:self.property_max_synonyms]))
            add_property(cid, 'Molecular formula', chem['mf'])
            add_property(cid, 'Molecular weight', f"{chem['mw']} g/mol")
            if chem['wiki']:
                add_property(cid, 'Wikipedia', chem['wiki'])
            if chem['iupacname']:
                add_property(cid, 'IUPAC name', chem['iupacname'])
            if chem['cas']:
                add_property(cid, 'CAS number', ', '.join(chem['cas'][:self.property_max_cas]))
            add_property(cid, 'SMILES', chem['smiles'])
            add_property(cid, 'InChI', chem['inchi'])
            add_property(cid, 'InChI-key', chem['inchikey'])
            add_property(cid, 'Pubchem complexity', chem['complexity'])
            add_property(cid, 'Bertz complexity', chem['bertz_complexity'])

        for entry in self.logger.track(solubility, 'Adding solubility', transient=True):
            value = '\n'.join(f"{item['solvent_name']}: {item['value']}" for item in entry['solubility'])
            add_property(entry['cid'], 'Solubility', value)

        def build_crc_str(entry):
            crc_str = str(entry['value'])
            if entry.get('approx'):
                crc_str = '~ ' + crc_str
            if entry.get('decomposes', False):
                crc_str += ' (decomposes)'
            if entry.get('sublimes', False):
                crc_str += ' (sublimes)'
            return crc_str

        processed_crc_cids = set()
        for entry in self.logger.track(
            crc_inorganic + crc_organic,
            'Adding CRC constants',
            transient=True,
        ):
            cid = entry['cid']
            if cid in processed_crc_cids:
                continue

            if entry['physical_form']:
                add_property(cid, 'Appearence', entry['physical_form'])
            if entry['mp'] and entry['mp']['value'] is not None:
                add_property(cid, 'Melting point', build_crc_str(entry['mp']))
            if entry['bp'] and entry['bp']['value'] is not None:
                add_property(cid, 'Boiling point', build_crc_str(entry['bp']))
            if entry['density']:
                add_property(cid, 'Density', f"{entry['density']} g/cm^3")
            if entry.get('refractive_index') is not None:
                add_property(cid, 'Refractive index', entry['refractive_index'])

            processed_crc_cids.add(cid)

        for entry in self.logger.track(
            crc_flammability,
            'Adding CRC flammability data',
            transient=True,
        ):
            cid = entry['cid']
            if entry['flash_point']:
                fp = build_crc_str(entry['flash_point'])
                add_property(cid, 'Flash point', f"{fp} °C")
            if entry['flash_limits']:
                add_property(cid, 'Flash limits', entry['flash_limits'])
            if entry['ignition_temp']:
                it = build_crc_str(entry['ignition_temp'])
                add_property(cid, 'Ignition temperature', f"{it} °C")

        return [{'cid': cid, 'properties': props} for cid, props in result.items()]

    def generate_curiosity_index(self):
        hazards = self.store.load_jsonl(self.hazards.chems_hazards_fn)

        cid_to_curiosity = dict()
        for chem in self.compounds.chems:
            cid_to_curiosity[chem['cid']] = chem['complexity'] / self.compounds.complexity_thr

        for entry in hazards:
            cid = entry['cid']

            nfpa_coeff = 0
            if entry['nfpa'] is not None:
                nfpa_coeff = sum(entry['nfpa']['value'].values()) / 12

            pictograms_coeff = 0
            if entry['pictograms'] is not None:
                pictograms_coeff = len(entry['pictograms']['value']) / len(self.LLM_HAZARD_CATEGORIES)

            if cid in cid_to_curiosity:
                cid_to_curiosity[cid] += nfpa_coeff + pictograms_coeff

        chems_occurence = self.get_chems_reactions_occurence()
        occur_values = np.array(list(chems_occurence.values()), dtype=float)
        sorted_vals = np.sort(occur_values)
        n = len(sorted_vals)

        def percentile_rank(value):
            return np.searchsorted(sorted_vals, value, side='right') / n

        for cid, occurence in chems_occurence.items():
            if cid in cid_to_curiosity:
                cid_to_curiosity[cid] += 3 * percentile_rank(occurence)

        max_curiosity = max(cid_to_curiosity.values())
        curiosity_entries = [
            {'cid': cid, 'curiosity': curiosity / max_curiosity}
            for cid, curiosity in cid_to_curiosity.items()
        ]
        self.store.write_jsonl(curiosity_entries, self.chems_curiosity_fn)

    def get_parsed_reactions_participants_norm_names(self):
        norm_names = set()
        for react in self.parsed_reactions:
            norm_names.update(entry['norm_name'] for entry in react['reagents'] + react['products'])
        return norm_names
