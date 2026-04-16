import json
import os
import re

from random import sample

from rdkit import Chem


class OrdParser:
    def __init__(self, data_dir, compounds, reactions, store, logger, ord_data='submodules/ord-data/'):
        self.data_dir = data_dir
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger

        self.ord_data = ord_data
        self.ord_source = "ord"

        self.reactions_parsed_ord_fn = os.path.join(self.reactions.parsed_reactions_dir, 'ord')
        os.makedirs(self.reactions_parsed_ord_fn, exist_ok=True)

        self.reactions_details_ord_fn = os.path.join(self.reactions.reactions_details_dir, 'ord')
        os.makedirs(self.reactions_details_ord_fn, exist_ok=True)

        self.smiles_cid_map_ord_fn = os.path.join(self.data_dir, 'misc', 'smiles_cid_map_ord.jsonl')
        self.unmapped_smiles_fn = os.path.join(self.data_dir, 'unmapped_smiles.jsonl')

        store.register_sorting(self.reactions_parsed_ord_fn, 'rid')
        store.register_sorting(self.reactions_details_ord_fn, 'rid')
        store.register_sorting(self.unmapped_smiles_fn, ('count', True))
        store.register_sorting(self.smiles_cid_map_ord_fn, 'cid')

        store.register_vault(self.reactions_parsed_ord_fn, 'reactions_')
        store.register_vault(self.reactions_details_ord_fn, 'reactions_details_')

    def _load_ord_helpers(self):
        from ord_schema.message_helpers import load_message
        from ord_schema.proto import dataset_pb2
        from google.protobuf.json_format import MessageToJson

        return load_message, dataset_pb2, MessageToJson

    def extract_file(self, in_fn, out_fn):
        load_message, dataset_pb2, MessageToJson = self._load_ord_helpers()

        dataset = load_message(
            os.path.join(self.ord_data, in_fn),
            dataset_pb2.Dataset,
        )

        rxn = dataset.reactions[0]
        rxn_json = json.loads(
            MessageToJson(
                message=rxn,
                including_default_value_fields=False,
                preserving_proto_field_name=True,
                indent=2,
                sort_keys=False,
                use_integers_for_enums=False,
                descriptor_pool=None,
                float_precision=None,
                ensure_ascii=True,
            )
        )

        with open(out_fn, 'w') as f:
            f.write(json.dumps(rxn_json, indent=2))

    def extract_ord_reactions(self, out_fn):
        load_message, dataset_pb2, MessageToJson = self._load_ord_helpers()

        reactions_written = 0
        overall = 0
        for dirpath, _, filenames in os.walk(self.ord_data):
            for fn in filenames:
                file_path = os.path.join(dirpath, fn)
                if not file_path.endswith(".pb.gz"):
                    continue

                try:
                    dataset = load_message(file_path, dataset_pb2.Dataset)
                    for reaction in dataset.reactions:
                        bad = False
                        reaction_json = json.loads(
                            MessageToJson(
                                message=reaction,
                                including_default_value_fields=False,
                                preserving_proto_field_name=True,
                                indent=2,
                                sort_keys=False,
                                use_integers_for_enums=False,
                                descriptor_pool=None,
                                float_precision=None,
                                ensure_ascii=True,
                            )
                        )
                        inputs_dict = reaction_json["inputs"]
                        inputs = []
                        for key in list(inputs_dict.keys()):
                            for input_entry in inputs_dict[key]["components"]:
                                if input_entry['reaction_role'] == "REACTANT":
                                    identifiers = [
                                        item for item in input_entry['identifiers'] if item['type'] == "SMILES"
                                    ]
                                    if not identifiers:
                                        bad = True
                                        break
                                    inputs.append(identifiers[0]['value'])

                        if bad:
                            continue

                        outcomes_list = reaction_json["outcomes"]
                        if len(outcomes_list) > 1:
                            print(file_path)
                            continue

                        outcomes = []
                        for prod in outcomes_list[0]["products"]:
                            prod_smiles = [item for item in prod['identifiers'] if item['type'] == "SMILES"]
                            if not prod_smiles:
                                bad = True
                                break
                            outcomes.append(prod_smiles[0]['value'])

                        if bad:
                            continue

                        chems = inputs + outcomes
                        criteria = True
                        for chem in chems:
                            mol = Chem.MolFromSmiles(chem)
                            if not mol:
                                bad = True
                                break

                            if not self.compounds.unmapped_complexity_criteria_mol(mol):
                                criteria = False
                                break

                        if bad:
                            continue

                        overall += 1
                        if criteria:
                            with open(out_fn, 'a') as f:
                                f.write(json.dumps(reaction_json) + '\n')
                            reactions_written += 1
                            if reactions_written % 1000 == 0:
                                self.logger.log(f"Written: {reactions_written}; All: {overall}")

                except Exception as exc:
                    self.logger.log_err(f"Exception occured: {exc}")

    def split_ord_file(self, ord_file, out_dir, prefix="ord", lines_per_file=20000):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        with open(ord_file) as f_in:
            i = 1
            while True:
                f_index = i // lines_per_file
                out_fn = os.path.join(out_dir, f"{prefix}_{f_index}.jsonl")
                with open(out_fn, 'w') as f_out:
                    print(f"Writing in '{out_fn}'")
                    while (i % lines_per_file) != 0:
                        line = f_in.readline()
                        if not line:
                            return

                        f_out.write(line)
                        i += 1
                    i += 1

    def sample_ord(self, ord_file, samples_num, out_fn):
        with open(ord_file) as f:
            reactions = [json.loads(x) for x in f.read().strip().split('\n')]

        res = sample(reactions, samples_num)
        with open(out_fn, 'w') as f:
            f.write(json.dumps(res, indent=1))

    def clean_ord_reactions(self, ord_files_dir, out_fn):
        fixed_cnt = 0
        with open(out_fn, 'w') as f_out:
            for dirpath, _, filenames in os.walk(ord_files_dir):
                for fn in filenames:
                    with open(os.path.join(dirpath, fn)) as f:
                        for line in f:
                            bad = False

                            reaction = json.loads(line)
                            inputs_dict = reaction["inputs"]
                            reactants = []
                            solvents = []
                            catalysts = []
                            for key in list(inputs_dict.keys()):
                                for input_entry in inputs_dict[key]["components"]:
                                    input_smiles = [
                                        item for item in input_entry['identifiers'] if item['type'] == "SMILES"
                                    ]
                                    if not input_smiles:
                                        bad = True
                                        break

                                    input_name = [
                                        item for item in input_entry['identifiers'] if item['type'] == "NAME"
                                    ]
                                    entry = {'smiles': input_smiles[0]['value']}
                                    entry['name'] = input_name[0]['value'] if input_name else None
                                    if input_entry['reaction_role'] == "REACTANT":
                                        reactants.append(entry)
                                    elif input_entry['reaction_role'] == "SOLVENT":
                                        solvents.append(entry)
                                    elif input_entry['reaction_role'] == "CATALYST":
                                        catalysts.append(entry)

                            if bad:
                                continue

                            products_raw = reaction["outcomes"][0]["products"]
                            products = []
                            for prod in products_raw:
                                prod_smiles = [
                                    item for item in prod['identifiers'] if item['type'] == "SMILES"
                                ][0]['value']
                                prod_name = [item for item in prod['identifiers'] if item['type'] == "NAME"]

                                entry = {'smiles': prod_smiles}
                                entry['name'] = prod_name[0]['value'] if prod_name else None
                                products.append(entry)

                            def fix_chems(chems):
                                nonlocal fixed_cnt
                                fixes = {
                                    "[Li][CH2]CCC": "[Li+].CCC[CH2-]",
                                    "[H][H]": "[HH]",
                                    "[Al+3].[H-].[H-].[H-].[H-].[Li+]": "[Li+].[AlH4-]",
                                    "[Cu][I]": "[Cu+].[I-]",
                                    "[Al+3].[Cl-].[Cl-].[Cl-]": "[Al](Cl)(Cl)Cl",
                                    "O.[Cl-].[Na+]": "[O-]Cl.[Na+]",
                                    "[Mn+4].[O-2].[O-2]": "O=[Mn]=O",
                                    "CC(C)[CH2][Al+][CH2]C(C)C.[H-]": "[H-].CC(C)C[Al+]CC(C)C",
                                    "N#[C][Cu]": "[C-]#N.[Cu+]",
                                    "N#[C][Na]": "[C-]#N.[Na+]",
                                    "[C].[Pd]": "C.[Pd]",
                                }

                                fixed = []
                                negative = []
                                positive = []
                                for chem in chems:
                                    smiles = chem['smiles']
                                    mol = Chem.MolFromSmiles(smiles)
                                    if not mol:
                                        return None
                                    smiles = Chem.MolToSmiles(mol)
                                    if smiles in fixes:
                                        smiles = fixes[smiles]
                                    chem['smiles'] = smiles
                                    charge = Chem.GetFormalCharge(mol)
                                    if charge == 0:
                                        fixed.append(chem)
                                    elif charge > 0:
                                        positive.append(chem)
                                    else:
                                        negative.append(chem)

                                if len(negative) == 1 and len(positive) == 1:
                                    smiles = f"{positive[0]['smiles']}.{negative[0]['smiles']}"
                                    mol = Chem.MolFromSmiles(smiles)
                                    if not mol or Chem.GetFormalCharge(mol) != 0:
                                        return None
                                    smiles = Chem.MolToSmiles(mol)
                                    if smiles in fixes:
                                        smiles = fixes[smiles]

                                    fixed_cnt += 1
                                    fixed.append({'smiles': smiles})

                                return fixed

                            reactants = fix_chems(reactants)
                            if not reactants:
                                continue
                            products = fix_chems(products)
                            if not products:
                                continue

                            description = None
                            if "notes" in reaction and "procedure_details" in reaction["notes"]:
                                description = reaction["notes"]['procedure_details']

                            provenance = dict()
                            if "provenance" in reaction:
                                provenance_field = reaction['provenance']
                                provenance['doi'] = provenance_field['doi'] if 'doi' in provenance_field else None
                                provenance['patent'] = (
                                    provenance_field['patent'] if 'patent' in provenance_field else None
                                )

                            cleaned_entry = {
                                "reagents": reactants,
                                "products": products,
                                "solvents": solvents if solvents else None,
                                "catalysts": catalysts if catalysts else None,
                                "description": description,
                                "provenance": provenance if provenance else None,
                            }

                            f_out.write(json.dumps(cleaned_entry) + '\n')

        self.logger.log(f"Total fixed: {fixed_cnt}")

    def parse_raw_ord_reactions(self, ord_reactions_fn):
        smiles_cid_ord_map = self.store.load_jsonl_map(self.smiles_cid_map_ord_fn, 'smiles', 'cid')
        smiles_cid_ord_map = {
            smiles: cid
            for smiles, cid in smiles_cid_ord_map.items()
            if cid != self.compounds.null_cid and cid in self.compounds.cid_chem_map
        }

        unmapped_smiles = dict()
        smiles_name_map = dict()
        mapped_count = 0
        overall_count = 0
        processed_rids = set()
        parsed_reactions_list = []
        parsed_details_list = []
        with open(ord_reactions_fn) as f_in:
            total = self.compounds.count_file_lines(ord_reactions_fn)
            for line in self.logger.track(f_in, "Parsing ORD reactions", total=total):
                overall_count += 1
                reaction = json.loads(line)
                parse_success = True

                def process_substance(substance, target_list):
                    nonlocal parse_success

                    def process_smiles(smiles):
                        parts = smiles.split('.')
                        if len(parts) < 2:
                            return smiles

                        filtered = []
                        for part in parts:
                            mol = Chem.MolFromSmiles(part)
                            if mol is None:
                                filtered.append(part)
                                continue
                            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                            if formula not in ('H2O', 'HOH', 'O'):
                                filtered.append(part)
                        return '.'.join(filtered)

                    smiles = process_smiles(substance['smiles'])
                    mol = Chem.MolFromSmiles(smiles)
                    if not mol:
                        parse_success = False
                        return None

                    smiles = Chem.MolToSmiles(mol, canonical=True)
                    norm_inchi = self.compounds.get_mol_norm_inchi(mol)
                    if not norm_inchi:
                        parse_success = False
                        return None

                    def process_name(name):
                        if not name:
                            return None

                        name = re.sub(
                            r'\s*(mono|di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?hydrate$',
                            '',
                            name,
                            flags=re.IGNORECASE,
                        )
                        return self.compounds.normalize_chem_name(name.strip(), is_clean=True)

                    original_name = substance.get('name')
                    if smiles in smiles_cid_ord_map:
                        cid = smiles_cid_ord_map[smiles]
                    elif norm_inchi in self.compounds.norm_inchi_cid_map:
                        cid = self.compounds.norm_inchi_cid_map[norm_inchi]
                    elif smiles in self.compounds.smiles_cid_map:
                        cid = self.compounds.smiles_cid_map[smiles]
                    else:
                        unmapped_smiles[smiles] = unmapped_smiles.setdefault(smiles, 0) + 1
                        if (
                            smiles_name_map.get(smiles) is None
                            and self.compounds.good_name_criteria(original_name)
                        ):
                            smiles_name_map[smiles] = original_name

                        parse_success = False
                        return None

                    chem = self.compounds.cid_chem_map[cid]
                    name = chem['cmpdname']
                    norm_name = self.compounds.normalize_chem_name(name)
                    target_list.append({'norm_name': norm_name, 'original_name': name, 'cid': cid})

                if reaction["description"] is None or len(reaction["description"]) < 150:
                    continue

                provenance = reaction['provenance']
                if provenance is None or provenance['patent'] is None or provenance['doi'] is None:
                    continue

                parsed_reaction = {'reagents': [], 'products': []}
                for reagent in reaction['reagents']:
                    process_substance(reagent, parsed_reaction['reagents'])
                for product in reaction['products']:
                    process_substance(product, parsed_reaction['products'])

                parsed_details = {'solvents': [], 'catalysts': []}
                if reaction['solvents']:
                    for solvent in reaction['solvents']:
                        process_substance(solvent, parsed_details['solvents'])
                if reaction['catalysts']:
                    for catalyst in reaction['catalysts']:
                        process_substance(catalyst, parsed_details['catalysts'])

                parsed_details['provenance'] = reaction['provenance']
                parsed_details['description'] = reaction['description']
                if parsed_details['provenance'] is None:
                    parsed_details['provenance'] = {'doi': None, 'patent': None}

                def deduplicate_part(entry, part):
                    entry[part] = list({chem['cid']: chem for chem in (entry.get(part) or [])}.values())

                if parse_success:
                    deduplicate_part(parsed_details, 'solvents')
                    deduplicate_part(parsed_details, 'catalysts')
                    deduplicate_part(parsed_reaction, 'reagents')
                    deduplicate_part(parsed_reaction, 'products')

                    parsed_reaction = self.reactions.assemble_reaction(parsed_reaction)
                    if parsed_reaction['rid'] in processed_rids:
                        continue

                    processed_rids.add(parsed_reaction['rid'])
                    parsed_reaction["source"] = self.ord_source
                    parsed_reaction['confidence'] = None

                    if not self.reactions.validate_reaction(parsed_reaction):
                        continue

                    parsed_details["rid"] = parsed_reaction['rid']
                    parsed_details["source"] = self.ord_source
                    parsed_details['confidence'] = None

                    parsed_reactions_list.append(parsed_reaction)
                    parsed_details_list.append(parsed_details)

                    mapped_count += 1
                    if mapped_count % 10000 == 0:
                        self.logger.log(f"Mapped {mapped_count} reactions out of {overall_count}")

        self.store.write_jsonl(parsed_reactions_list, self.reactions_parsed_ord_fn)
        self.store.write_jsonl(parsed_details_list, self.reactions_details_ord_fn)

        unmapped_smiles_list = [
            {'smiles': smiles, 'name': smiles_name_map.get(smiles), 'count': unmapped_smiles[smiles]}
            for smiles in unmapped_smiles
        ]
        self.store.write_jsonl(unmapped_smiles_list, self.unmapped_smiles_fn)

        self.logger.log(f"Unmapped smiles: {len(unmapped_smiles)}")
        self.logger.log(f"Mapped {mapped_count} reactions out of {overall_count}")

    def add_unmapped_ord_chems(self, count_thr=2):
        unmapped = self.store.load_jsonl(self.unmapped_smiles_fn)
        name_smiles = []
        for entry in unmapped:
            if entry['count'] < count_thr:
                continue
            name_smiles.append((entry['name'], entry['smiles']))

        self.compounds.add_new_chems(name_smiles)

    def populate_unmapped_smiles_cid_map(self):
        smiles_cid_ord_map = self.store.load_jsonl_map(self.smiles_cid_map_ord_fn, 'smiles', 'cid')
        unmapped_entries = self.store.load_jsonl(self.unmapped_smiles_fn)

        try:
            def process_entry(entry, skip_prompt=False):
                smiles = entry['smiles']
                if smiles in smiles_cid_ord_map:
                    return True

                name = entry['name']
                if not name:
                    return True

                norm_name = self.compounds.normalize_chem_name(name, is_clean=True)
                if norm_name in self.compounds.unique_name_cid_map:
                    unmapped_chem = self.compounds.build_chem(-1, name, smiles)
                    if not unmapped_chem:
                        return True

                    cid = self.compounds.unique_name_cid_map[norm_name]
                    mapped_chem = self.compounds.cid_chem_map[cid]

                    mapped_chem_inchi = self.compounds.get_chem_norm_inchi(mapped_chem)
                    unmapped_chem_inchi = self.compounds.get_chem_norm_inchi(unmapped_chem)

                    while True:
                        if mapped_chem_inchi != unmapped_chem_inchi or cid < 0:
                            if skip_prompt:
                                return True

                            self.compounds.display_compare_table(
                                "[bold yellow]Confirm merge?[/bold yellow]",
                                unmapped_chem,
                                mapped_chem,
                            )
                            decision = input("Merge? [yn]: ").lower().strip()
                        else:
                            decision = 'y'

                        if decision in {'y', 'yes'}:
                            smiles_cid_ord_map[smiles] = cid
                            self.logger.print(f"Merging '{name}' with CID {cid}")
                        elif decision in {'n', 'no'}:
                            smiles_cid_ord_map[smiles] = self.compounds.null_cid
                        elif decision in {'exit', 'stop'}:
                            return False
                        else:
                            self.logger.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                            continue
                        break

                return True

            for entry in unmapped_entries:
                process_entry(entry, skip_prompt=True)

            for entry in unmapped_entries:
                proceed = process_entry(entry, skip_prompt=False)
                if not proceed:
                    break

        finally:
            self.store.write_jsonl_map(smiles_cid_ord_map, 'smiles', 'cid', self.smiles_cid_map_ord_fn)
