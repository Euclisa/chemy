import os


class SqlExporter:
    def __init__(self, compounds, reactions, store, logger, properties, misc, hazards, db_name=None):
        self.compounds = compounds
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.properties = properties
        self.misc = misc
        self.hazards = hazards

        self.db_name = db_name if db_name is not None else os.getenv("CHEMS_DB_NAME")

    def populate_db(self):
        import psycopg2
        from psycopg2.extras import execute_values

        with psycopg2.connect(database=self.db_name) as conn, self.logger.progress() as progress:
            with conn.cursor() as cur:
                task = progress.add_task("Populating db", total=20)

                def execute_values_advance(cursor, sql, data):
                    execute_values(cursor, sql, data)
                    progress.advance(task)

                all_cids = {chem['cid'] for chem in self.compounds.chems}

                sql = (
                    "INSERT INTO compounds (cid, name, mf, mw, charge, smiles, inchi, inchikey, "
                    "complexity, bertz_complexity, organic) VALUES %s ON CONFLICT (cid) DO NOTHING "
                )
                execute_values_advance(
                    cur,
                    sql,
                    [
                        (
                            chem['cid'],
                            chem['cmpdname'],
                            chem['mf'],
                            chem['mw'],
                            chem['charge'],
                            chem['smiles'],
                            chem['inchi'],
                            chem['inchikey'],
                            chem['complexity'],
                            chem['bertz_complexity'],
                            chem['organic'],
                        )
                        for chem in self.compounds.chems
                    ],
                )

                sql = (
                    "INSERT INTO compound_synonyms (cid, synonym) VALUES %s "
                    "ON CONFLICT (cid, synonym) DO NOTHING"
                )
                data = [
                    (chem['cid'], syn)
                    for chem in self.compounds.chems_mapped
                    for syn in chem['cmpdsynonym']
                    if syn
                ]
                execute_values_advance(cur, sql, data)

                sql = "INSERT INTO compound_fingerprints (cid, bits, popcount) VALUES %s "
                data = [
                    (chem['cid'], chem['ECFP4_fp']['bits'], chem['ECFP4_fp']['popcount'])
                    for chem in self.compounds.chems
                ]
                execute_values_advance(cur, sql, data)

                sql = "INSERT INTO compound_cas (cid, cas) VALUES %s"
                data = [
                    (chem['cid'], cas)
                    for chem in self.compounds.chems
                    for cas in chem['cas']
                ]
                execute_values_advance(cur, sql, data)

                sql = "INSERT INTO compound_wiki (cid, wiki) VALUES %s "
                data = [
                    (entry['cid'], entry['wiki'])
                    for entry in self.compounds.chems
                    if entry['wiki'] is not None
                ]
                execute_values_advance(cur, sql, data)

                sql_nfpa = (
                    "INSERT INTO compound_nfpa (cid, health, flammability, instability) VALUES %s"
                )
                sql_pictograms = (
                    "INSERT INTO compound_hazard_pictograms (cid, pictogram) VALUES %s"
                )

                hazards = self.store.load_jsonl(self.hazards.chems_hazards_fn)
                nfpa_data = []
                pictograms_data = []
                for entry in hazards:
                    cid = entry['cid']
                    if cid not in all_cids:
                        continue

                    if entry['nfpa']:
                        nfpa = entry['nfpa']['value']
                        nfpa_data.append(
                            (cid, nfpa.get('health'), nfpa.get('flammability'), nfpa.get('instability'))
                        )

                    if entry['pictograms']:
                        for pic in entry['pictograms']['value']:
                            pictograms_data.append((cid, pic))

                execute_values_advance(cur, sql_nfpa, nfpa_data)
                execute_values_advance(cur, sql_pictograms, pictograms_data)

                sql = "INSERT INTO categories (code, name, domain) VALUES %s"
                categories_data = self.store.load_jsonl(self.compounds.categories_fn)
                categories_data = [(c['code'], c['name'], c['domain']) for c in categories_data]
                execute_values_advance(cur, sql, categories_data)

                sql = "INSERT INTO compound_categories (cid, category_code) VALUES %s"
                chems_categories_data = self.store.load_jsonl(self.compounds.chems_categories_fn)
                chems_categories_data = [
                    (entry['cid'], category)
                    for entry in chems_categories_data
                    for category in entry['categories']
                    if entry['cid'] in all_cids
                ]
                execute_values_advance(cur, sql, chems_categories_data)

                sql = (
                    "INSERT INTO compound_properties (cid, property_name, property_value, rank) VALUES %s"
                )
                chems_properties_data = self.properties.compile_chems_properties()
                chems_properties_data = [
                    (entry['cid'], prop['property'], prop['value'], rank)
                    for entry in chems_properties_data
                    for rank, prop in enumerate(entry['properties'])
                ]
                execute_values_advance(cur, sql, chems_properties_data)

                sql = "INSERT INTO compound_descriptions (cid, description) VALUES %s"
                chems_descriptions_fn = os.path.join(
                    self.compounds.chems_properties_dir,
                    'llm',
                    'chems_descriptions.jsonl',
                )
                chems_descriptions_data = self.store.load_jsonl(chems_descriptions_fn)
                chems_descriptions_data = [
                    (entry['cid'], entry['description'])
                    for entry in chems_descriptions_data
                    if entry['cid'] in all_cids
                ]
                execute_values_advance(cur, sql, chems_descriptions_data)

                sql = "INSERT INTO compound_commonness_sorting (cid, rank) VALUES %s"
                commonness_sorting_data = self.misc.get_commonness_chems_sorting()
                commonness_sorting_data = [
                    (cid, rank) for rank, cid in enumerate(commonness_sorting_data)
                ]
                execute_values_advance(cur, sql, commonness_sorting_data)

                sql = "INSERT INTO compound_complexity_sorting (cid, rank) VALUES %s"
                complexity_sorting_data = self.misc.get_complexity_chems_sorting()
                complexity_sorting_data = [
                    (cid, rank) for rank, cid in enumerate(complexity_sorting_data)
                ]
                execute_values_advance(cur, sql, complexity_sorting_data)

                sql = "INSERT INTO compound_curiosity_sorting (cid, rank) VALUES %s"
                curiosity_sorting_data = self.misc.get_curiosity_chems_sorting()
                curiosity_sorting_data = [
                    (cid, rank) for rank, cid in enumerate(curiosity_sorting_data)
                ]
                execute_values_advance(cur, sql, curiosity_sorting_data)

                rid_details_map = {entry['rid']: entry for entry in self.properties.reactions_details}

                sql = (
                    "INSERT INTO reactions (rid, complexity, source, balanced, confidence) VALUES %s"
                )
                sql_reactants = "INSERT INTO reaction_reactants (cid, coeff, rid) VALUES %s"
                sql_products = "INSERT INTO reaction_products (cid, coeff, rid) VALUES %s"
                sql_solvents = "INSERT INTO reaction_solvents (cid, rid) VALUES %s"
                sql_catalysts = "INSERT INTO reaction_catalysts (cid, rid) VALUES %s"
                sql_details = (
                    "INSERT INTO reaction_details (rid, doi, patent, description, source, confidence) "
                    "VALUES %s"
                )

                reactions = []
                reaction_reactants = []
                reaction_products = []
                reaction_solvents = []
                reaction_catalysts = []
                reaction_details = []
                for react in self.properties.parsed_reactions:
                    rid = react['rid']
                    balance = self.reactions.reactions_balance.get(rid)

                    reactions.append(
                        (
                            rid,
                            react['complexity'],
                            react['source'],
                            self.reactions.is_react_balanced(react),
                            react['confidence'],
                        )
                    )

                    for entry in react['reagents']:
                        coeff = None if not balance else balance[entry['cid']]
                        reaction_reactants.append((entry['cid'], coeff, rid))

                    for entry in react['products']:
                        coeff = None if not balance else balance[entry['cid']]
                        reaction_products.append((entry['cid'], coeff, rid))

                    if rid in rid_details_map:
                        for entry in rid_details_map[rid]['solvents']:
                            reaction_solvents.append((entry['cid'], rid))

                        for entry in rid_details_map[rid]['catalysts']:
                            reaction_catalysts.append((entry['cid'], rid))

                        reaction_details.append(
                            (
                                rid,
                                rid_details_map[rid]['provenance']['doi'],
                                rid_details_map[rid]['provenance']['patent'],
                                rid_details_map[rid]['description'],
                                rid_details_map[rid]['source'],
                                rid_details_map[rid]['confidence'],
                            )
                        )

                execute_values_advance(cur, sql, reactions)
                execute_values_advance(cur, sql_reactants, reaction_reactants)
                execute_values_advance(cur, sql_products, reaction_products)
                execute_values_advance(cur, sql_solvents, reaction_solvents)
                execute_values_advance(cur, sql_catalysts, reaction_catalysts)
                execute_values_advance(cur, sql_details, reaction_details)

                conn.commit()
