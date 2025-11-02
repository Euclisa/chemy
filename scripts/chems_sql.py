import psycopg2
from psycopg2.extras import execute_values
import json
import os

from chems_misc import ChemsMisc


class ChemsSql(ChemsMisc):
    def __init__(self, data_dir, db_name=None):
        super().__init__(data_dir)

        if db_name is None:
            self.db_name = os.getenv("CHEMS_DB_NAME")
        else:
            self.db_name = db_name


    def populate_db(self):

        with self._rich_progress() as progress:
            with psycopg2.connect(database=self.db_name) as conn:
                with conn.cursor() as cur:
                    task = progress.add_task("Populating db", total=20)

                    def execute_values_advance(cursor, sql, data):
                        execute_values(cursor, sql, data)
                        progress.advance(task)

                    all_cids = set([chem['cid'] for chem in self.chems])

                    sql = \
                    "INSERT INTO compounds (cid, name, mf, mw, charge, smiles, inchi, inchikey, complexity, bertz_complexity, organic) " \
                    "VALUES %s " \
                    "ON CONFLICT (cid) DO NOTHING "

                    
                    execute_values_advance(cur, sql, [(chem['cid'],
                                                chem['cmpdname'],
                                                chem['mf'],
                                                chem['mw'],
                                                chem['charge'],
                                                chem['smiles'],
                                                chem['inchi'],
                                                chem['inchikey'],
                                                chem['complexity'],
                                                chem['bertz_complexity'],
                                                chem['organic']) for chem in self.chems])

                    sql = \
                    "INSERT INTO compound_synonyms (cid, synonym) " \
                    "VALUES %s " \
                    "ON CONFLICT (cid, synonym) DO NOTHING"
                    data = [(chem['cid'], syn) for chem in self.chems for syn in chem['cmpdsynonym'] if syn]
                    execute_values_advance(cur, sql, data)


                    sql = \
                    "INSERT INTO compound_fingerprints (cid, ECFP4_fp, popcount) " \
                    "VALUES %s "
                    data = [(chem['cid'], chem['ECFP4_fp']['bits'], chem['ECFP4_fp']['popcount']) for chem in self.chems]
                    execute_values_advance(cur, sql, data)



                    sql = \
                    "INSERT INTO compound_cas (cid, cas) " \
                    "VALUES %s"
                    data = [(chem['cid'], cas) for chem in self.chems for cas in chem['cas']]
                    execute_values_advance(cur, sql, data)


                    sql = \
                    "INSERT INTO compound_wiki (cid, wiki) " \
                    "VALUES %s "
                    data = [(x['cid'], x['wiki']) for x in self.chems if x['wiki'] is not None]
                    execute_values_advance(cur, sql, data)


                    sql_nfpa = \
                    "INSERT INTO compound_nfpa (cid, health, flammability, instability) " \
                    "VALUES %s"
                    sql_pictograms = \
                    "INSERT INTO compound_hazard_pictograms (cid, pictogram) " \
                    "VALUES %s"

                    hazards = self._load_jsonl(self.chems_hazards_fn)
                    nfpa_data = []
                    pictograms_data = []
                    for entry in hazards:
                        cid = entry['cid']
                        if cid not in all_cids:
                            continue

                        if entry['nfpa']:
                            nfpa = entry['nfpa']['value']
                            nfpa_data.append((cid, nfpa.get('health'), nfpa.get('flammability'), nfpa.get('instability')))
                        
                        if entry['pictograms']:
                            for pic in entry['pictograms']['value']:
                                pictograms_data.append((cid, pic))

                    execute_values_advance(cur, sql_nfpa, nfpa_data)
                    execute_values_advance(cur, sql_pictograms, pictograms_data)


                    sql = \
                    "INSERT INTO categories (code, name, domain) " \
                    "VALUES %s"
                    categories_data = self._load_jsonl(self.categories_fn)
                    categories_data = [(c['code'], c['name'], c['domain']) for c in categories_data]
                    execute_values_advance(cur, sql, categories_data)


                    sql = \
                    "INSERT INTO compound_categories (cid, category_code) " \
                    "VALUES %s"
                    chems_categories_data = self._load_jsonl(self.chems_categories_fn)
                    chems_categories_data = [(entry['cid'], cat) for entry in chems_categories_data for cat in entry['categories'] if entry['cid'] in all_cids]
                    execute_values_advance(cur, sql, chems_categories_data)


                    sql = \
                    "INSERT INTO compound_properties (cid, property_name, property_value) " \
                    "VALUES %s"
                    chems_properties_data = self._compile_chems_properties()
                    chems_properties_data = [(entry['cid'], prop['property'], prop['value']) for entry in chems_properties_data for prop in entry['properties']]
                    execute_values_advance(cur, sql, chems_properties_data)

                    
                    sql = \
                    "INSERT INTO compound_descriptions (cid, description) " \
                    "VALUES %s"
                    chems_descriptions_data = self._load_jsonl(self.chems_descriptions_fn)
                    chems_descriptions_data = [(entry['cid'], entry['description']) for entry in chems_descriptions_data if entry['cid'] in all_cids]
                    execute_values_advance(cur, sql, chems_descriptions_data)


                    sql = \
                    "INSERT INTO compound_commonness_sorting (cid, rank) " \
                    "VALUES %s"
                    commonness_sorting_data = self._get_commonnes_chems_sorting()
                    commonness_sorting_data = [(cid, rank) for rank, cid in enumerate(commonness_sorting_data)]
                    execute_values_advance(cur, sql, commonness_sorting_data)

                    sql = \
                    "INSERT INTO compound_complexity_sorting (cid, rank) " \
                    "VALUES %s"
                    complexity_sorting_data = self._get_complexity_chems_sorting()
                    complexity_sorting_data = [(cid, rank) for rank, cid in enumerate(complexity_sorting_data)]
                    execute_values_advance(cur, sql, complexity_sorting_data)

                    sql = \
                    "INSERT INTO compound_curiosity_sorting (cid, rank) " \
                    "VALUES %s"
                    curiosity_sorting_data = self._get_curiosity_chems_sorting()
                    curiosity_sorting_data = [(cid, rank) for rank, cid in enumerate(curiosity_sorting_data)]
                    execute_values_advance(cur, sql, curiosity_sorting_data)


                    
                    rid_details_map = dict()
                    for entry in self.reactions_details:
                        rid_details_map[entry['rid']] = entry
                    
                    
                    sql = \
                    "INSERT INTO reactions (rid, complexity, source, balanced, confidence) " \
                    "VALUES %s"
                    sql_reactants = \
                    "INSERT INTO reaction_reactants (cid, coeff, rid) " \
                    "VALUES %s"
                    sql_products = \
                    "INSERT INTO reaction_products (cid, coeff, rid) " \
                    "VALUES %s"
                    sql_solvents = \
                    "INSERT INTO reaction_solvents (cid, rid) " \
                    "VALUES %s"
                    sql_catalysts = \
                    "INSERT INTO reaction_catalysts (cid, rid) " \
                    "VALUES %s"
                    sql_details = \
                    "INSERT INTO reaction_details (rid, doi, patent, description, source, confidence) " \
                    "VALUES %s"

                    reactions = []
                    reaction_reactants = []
                    reaction_products = []
                    reaction_solvents = []
                    reaction_catalysts = []
                    reaction_details = []
                    for react in self.parsed_reactions:
                        rid = react['rid']
                        balance = self.reactions_balance.get(rid)

                        reactions.append((rid, react['complexity'], react['source'], self._is_react_balanced(react), react['confidence']))

                        for x in react['reagents']:
                            cid = x['cid']
                            coeff = None if not balance else balance[cid]
                            reaction_reactants.append((x['cid'], coeff, rid))
                        
                        for x in react['products']:
                            cid = x['cid']
                            coeff = None if not balance else balance[cid]
                            reaction_products.append((x['cid'], coeff, rid))
                        
                        if rid in rid_details_map:
                            for x in rid_details_map[rid]['solvents']:
                                reaction_solvents.append((x['cid'], rid))
                            
                            for x in rid_details_map[rid]['catalysts']:
                                reaction_catalysts.append((x['cid'], rid))

                            reaction_details.append(
                                (
                                rid,
                                rid_details_map[rid]['provenance']['doi'],
                                rid_details_map[rid]['provenance']['patent'],
                                rid_details_map[rid]['description'],
                                rid_details_map[rid]['source'],
                                rid_details_map[rid]['confidence'])
                                )
                        
                    execute_values_advance(cur, sql, reactions)
                    execute_values_advance(cur, sql_reactants, reaction_reactants)
                    execute_values_advance(cur, sql_products, reaction_products)
                    execute_values_advance(cur, sql_solvents, reaction_solvents)
                    execute_values_advance(cur, sql_catalysts, reaction_catalysts)
                    execute_values_advance(cur, sql_details, reaction_details)
    
                    conn.commit()



if __name__ == "__main__":
    sql = ChemsSql('data', 'chemistry')
    sql.populate_db()