from chems_llm_parse import ChemsLLMParse
from chems_ord_parse import ChemsOrdParse
from chems_sql import ChemsSql


class ChemsMain(ChemsLLMParse, ChemsOrdParse, ChemsSql):
    
    def __init__(self, data_dir):
        super().__init__(data_dir)
    
    def merge_parsed_files(self, out_fn, *parsed_reactions_files):
        rid_reaction = dict()
        total_reactions = 0
        for fn in parsed_reactions_files:
            reactions = self._load_jsonl(fn)
            
            total_reactions += len(reactions)
            
            for react in reactions:
                rid = react['rid']
                if rid not in rid_reaction:
                    rid_reaction[rid] = react
                else:
                    old_react = rid_reaction[rid]
                    old_source = old_react.get('source')
                    new_source = react.get('source')
                    new_source_priority = self.sources_priority.get(new_source, -1)
                    old_source_priority = self.sources_priority.get(old_source, -1)
                    if new_source_priority > old_source_priority:
                        rid_reaction[rid] = react

        reactions_res = list(rid_reaction.values())
        self._write_jsonl(reactions_res, out_fn, backup=False)
    

    def merge_reactions(self):
        self.merge_parsed_files(self.reactions_parsed_fn, self.reactions_parsed_llm_fn, self.reactions_parsed_ord_fn)
    
    def merge_details(self):
        self.merge_parsed_files(self.reactions_details_fn, self.reactions_details_llm_fn, self.reactions_details_ord_fn)
    

    def generate_edges(self):
        self.merge_reactions()

        reactions = self._load_jsonl(self.reactions_parsed_fn)
        
        edge_reaction_id_map = dict()
        for react in self._rich_track(reactions, "Generating edges..."):
            react_id = react['rid']
            for r in react['reagents']:
                r_cid = r['cid']
                for p in react['products']:
                    p_cid = p['cid']
                    edge = (r_cid, p_cid)
                    if edge not in edge_reaction_id_map:
                        edge_reaction_id_map[edge] = []
                    edge_reaction_id_map[edge].append(react_id)
        
        def get_eid(edge):
            source = edge[0]
            target = edge[1]
            return ((source+target)*target) % 2**64
        
        edges = []
        for edge in edge_reaction_id_map:
            entry = {'eid': get_eid(edge), 'source': edge[0], 'target': edge[1], 'reactions': edge_reaction_id_map[edge]}
            edges.append(entry)

        self._write_jsonl(edges, self.chems_edges_fn, backup=False)
        self.log(f"Generated {len(edge_reaction_id_map)} edges")
    

    def test(self):
        raw_verdict = self._load_jsonl(self.raw_reactions_verdict_fn)
        rid_entry = dict()
        cnt = 0
        for entry in raw_verdict:
            react = entry['reaction']
            react_parsed, _ = self._parse_reaction_str(react)
            if not react_parsed:
                continue

            rid = react_parsed['rid']
            if rid in rid_entry and rid_entry[rid]['confidence'] > entry['confidence']:
                continue
            
            entry['rid'] = react_parsed['rid']
            rid_entry[rid] = entry
            if entry['confidence'] >= 0.4:
                cnt += 1
        
        self.print(cnt)

        self._write_jsonl(list(rid_entry.values()), self.raw_reactions_verdict_fn)

        




if __name__ == "__main__":
    chems = ChemsMain('data/')
    chems.test()