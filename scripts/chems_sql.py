import psycopg2
from psycopg2.extras import execute_values
import json
import os

from chems_db import ChemsDB


class ChemsSql(ChemsDB):
    def __init__(self, data_dir, db_name=None):
        super().__init__(data_dir)

        if db_name is None:
            self.db_name = os.getenv("CHEMS_DB_NAME")


    def populate_db(self):
        with open(self.chems_fn) as f:
            chems = [json.loads(x) for x in f.read().strip().split('\n')]
        
        all_cids = set([chem['cid'] for chem in chems])
        
        conn = psycopg2.connect(database=self.db_name)
        cur = conn.cursor()

        sql = \
        "INSERT INTO compounds (cid, name, mf, mw, charge, smiles, inchi, inchikey, complexity, bertz_complexity, organic) " \
        "VALUES %s " \
        "ON CONFLICT (cid) DO NOTHING "

        
        execute_values(cur, sql, [(chem['cid'],
                                    chem['cmpdname'],
                                    chem['mf'],
                                    chem['mw'],
                                    chem['charge'],
                                    chem['smiles'],
                                    chem['inchi'],
                                    chem['inchikey'],
                                    chem['complexity'],
                                    chem['bertz_complexity'],
                                    chem['organic']) for chem in chems])

        sql = \
        "INSERT INTO compound_synonyms (cid, synonym) " \
        "VALUES %s " \
        "ON CONFLICT (cid, synonym) DO NOTHING"
        data = [(chem['cid'], syn) for chem in chems for syn in chem['cmpdsynonym'] if syn]
        execute_values(cur, sql, data)


        sql = \
        "INSERT INTO compound_fingerprints (cid, ECFP4_fp, popcount) " \
        "VALUES %s "
        data = [(chem['cid'], chem['ECFP4_fp']['bits'], chem['ECFP4_fp']['popcount']) for chem in chems]
        execute_values(cur, sql, data)


        sql = \
        "INSERT INTO compound_cas (cid, cas) " \
        "VALUES %s"
        data = [(chem['cid'], cas) for chem in chems for cas in chem['cas']]
        execute_values(cur, sql, data)


        sql = \
        "INSERT INTO compound_wiki (cid, wiki) " \
        "VALUES %s "
        data = [(x['cid'], x['wiki']) for x in chems if x['wiki'] is not None]
        execute_values(cur, sql, data)

        hazards = self._load_jsonl(self.hazards_chems_fn)
        sql = \
        "INSERT INTO compound_nfpa (cid, health, flammability, instability) " \
        "VALUES %s"

        nfpa_data = [(entry['cid'], entry['nfpa'].get('healthHazard'), entry['nfpa'].get('fireHazard'), entry['nfpa'].get('instability')) for entry in hazards if entry['cid'] in all_cids]
        execute_values(cur, sql, nfpa_data)


        sql = \
        "INSERT INTO compound_hazard_statements (cid, statement) " \
        "VALUES %s"
        statements_data = [(entry['cid'], statement) for entry in hazards for statement in entry['statements'] if entry['cid'] in all_cids]
        execute_values(cur, sql, statements_data)

        sql = \
        "INSERT INTO compound_hazard_pictograms (cid, pictogram) " \
        "VALUES %s"
        pictograms_data = [(entry['cid'], pic) for entry in hazards for pic in entry['pictograms'] if entry['cid'] in all_cids]
        execute_values(cur, sql, pictograms_data)

        categories = self._load_jsonl(self.categories_fn)
        sql = \
        "INSERT INTO categories (code, name, domain) " \
        "VALUES %s"
        categories_data = [(c['code'], c['name'], c['domain']) for c in categories]
        execute_values(cur, sql, categories_data)



        chems_categories = self._load_jsonl(self.chems_categories_fn)
        sql = \
        "INSERT INTO compound_categories (cid, category_code) " \
        "VALUES %s"
        chems_categories_data = [(entry['cid'], cat) for entry in chems_categories for cat in entry['categories'] if entry['cid'] in all_cids]
        execute_values(cur, sql, chems_categories_data)

        with open(self.chems_descriptions_fn) as f:
            chems_descriptions_data = []
            for x in f.read().strip().split('\n'):
                entry = json.loads(x)
                if entry['cid'] in all_cids:
                    chems_descriptions_data.append((entry['cid'], entry['description']))
        
        sql = \
        "INSERT INTO compound_descriptions (cid, description) " \
        "VALUES %s"
        execute_values(cur, sql, chems_descriptions_data)

        with open(self.background_cids_fn) as f:
            background_cids = [(cid,) for cid in json.loads(f.read())]
        
        sql = \
        "INSERT INTO background_compounds (cid) " \
        "VALUES %s"
        execute_values(cur, sql, background_cids)


        reactions = self._load_jsonl(self.reactions_parsed_fn)
        details = self._load_jsonl(self.reactions_details_fn)
        
        rid_local_id_map = dict()
        for i, react in enumerate(reactions):
            rid_local_id_map[react['rid']] = i+1
        
        rid_details_map = dict()
        for entry in details:
            rid_details_map[entry['rid']] = entry
        
        
        sql = \
        "INSERT INTO reactions (rid, complexity, source, balanced, confidence) " \
        "VALUES %s"
        sql_reactants = \
        "INSERT INTO reaction_reactants (cid, rid) " \
        "VALUES %s"
        sql_products = \
        "INSERT INTO reaction_products (cid, rid) " \
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
        execute_values(cur, sql, [(x['rid'], x['complexity'], x['source'], x['balanced'], x['confidence']) for x in reactions])
        execute_values(cur, sql_reactants, [(x['cid'], react['rid']) for react in reactions for x in react['reagents']])
        execute_values(cur, sql_products, [(x['cid'], react['rid']) for react in reactions for x in react['products']])
        execute_values(cur, sql_solvents, [(x['cid'], rid) for rid in rid_details_map for x in rid_details_map[rid]['solvents']  if rid in rid_local_id_map])
        execute_values(cur, sql_catalysts, [(x['cid'], rid) for rid in rid_details_map for x in rid_details_map[rid]['catalysts'] if rid in rid_local_id_map])
        execute_values(cur, sql_details, [(rid, rid_details_map[rid]['provenance']['doi'], rid_details_map[rid]['provenance']['patent'], rid_details_map[rid]['description'], rid_details_map[rid]['source'], rid_details_map[rid]['confidence']) for rid in rid_details_map if rid in rid_local_id_map])

        
        conn.commit()

        cur.close()
        conn.close()