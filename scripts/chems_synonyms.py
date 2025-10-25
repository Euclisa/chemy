from rich.table import Table
from rich.rule import Rule

from chems_reaction_properties import ChemsReactionProperties


class ChemsSynonyms(ChemsReactionProperties):

    def __init__(self, data_dir):
        super().__init__(data_dir)


    def __update_cids_filtered_synonyms(self, new_entries):
        for cid, synonyms in new_entries.items():
            if not synonyms:
                continue
            self.cids_filtered_synonyms.setdefault(cid, set()).update(synonyms)
        
        res_entries = [{'cid': cid, 'synonyms': sorted(list(syns))} for cid, syns in self.cids_filtered_synonyms.items()]
        self._write_jsonl(res_entries, self.cids_filtered_synonyms_fn)


    def resolve_conflicting_synonyms(self, only_relevant=True):

        print("Resolve options:")
        print("  s0   - discard conflicted synonym from both compounds")
        print("  s1   - retain conflicted synonym for 1st compound")
        print("  s2   - retain conflicted synonym for 2nd compound")
        print("  c1   - retain 1st compound")
        print("  c2   - retain 2nd compound")
        print("  exit - save all previous decisions and abort")
        print()

        def display_conflict(conflict_i, conflict_norm_name, cid1, cid2):
            
            def get_conflict_inds(chem):
                return [
                    i for i, x in enumerate(chem['cmpdsynonym'])
                    if self._normalize_chem_name(x, is_clean=True) == conflict_norm_name
                ]
            
            def get_conflict_synonyms_str(syns, conflict_inds):
                return '; '.join(f"'{syns[i]}' ({i}/{len(syns)})" for i in conflict_inds)
            
            def display_compound(compound_i, chem, cid):
                syns = chem['cmpdsynonym']
                conflict_inds = get_conflict_inds(chem)
                name = chem['cmpdname']
                inchi = chem['inchi']
                cas = chem['cas']
                conflict_synonyms_str = get_conflict_synonyms_str(syns, conflict_inds)
                syns_num_to_disp = 10
                top_syns_str = ', '.join(f'"{syn}"' for syn in syns[:syns_num_to_disp])

                table = Table(
                    title=f"[bold cyan]Compound {compound_i}: {name}[/bold cyan] (CID: [yellow]{cid}[/yellow])",
                    show_header=False,
                    expand=True
                )
                table.add_row("Conflict synonyms", f"[magenta]{conflict_synonyms_str}[/magenta]")
                table.add_row(f"First {syns_num_to_disp} synonyms", f"{top_syns_str}")
                table.add_row("InChI", f"{inchi}")
                table.add_row("CAS", f"{cas}")
                self.print(table)

            self.print(Rule(f"[bold yellow]CONFLICT {conflict_i}: '{conflict_norm_name}'[/bold yellow]"))

            chem1, chem2 = self.cid_chem_map[cid1], self.cid_chem_map[cid2]
            display_compound(1, chem1, cid1)
            self.print(Rule())
            display_compound(2, chem2, cid2)

        conflict_map = dict()
        cids_to_delete = set()
        cids_syns_to_del = dict()
        try:
            for chem in self.chems:
                cid = chem['cid']
                cids_syns_to_del[cid] = set()
                norm_syns = set(self._normalize_chem_name(syn, is_clean=True) for syn in chem['cmpdsynonym'])
                for norm_name in norm_syns:
                    conflict_map.setdefault(norm_name, []).append(cid)
            
            stop = False
            conflict_map = {name: cids for name, cids in conflict_map.items() if len(cids) > 1}
            if only_relevant:
                relevant_names = self._get_parsed_reactions_participants_norm_names()
                conflict_map = {name: cids for name, cids in conflict_map.items() if name in relevant_names}

            self.print()
            self.print(f"{len(conflict_map)} conflicting names pending resolution")
            self.print()
            conflict_i = 0
            for norm_name in conflict_map:
                conflict_cids = [cid for cid in conflict_map[norm_name] if cid not in cids_to_delete]
                if len(conflict_cids) <= 1:
                    continue

                while len(conflict_cids) >= 2:
                    cid1, cid2 = conflict_cids[0], conflict_cids[1]
                    conflict_i += 1

                    display_conflict(conflict_i, norm_name, cid1, cid2)

                    decision = input("* Decision: ").strip()
                    if decision == 's0':
                        cids_syns_to_del[cid1].add(norm_name)
                        cids_syns_to_del[cid2].add(norm_name)
                        conflict_cids = conflict_cids[2:]
                        self.print(f"* Removed synonym from both compounds")
                    elif decision == 's1':
                        cids_syns_to_del[cid2].add(norm_name)
                        conflict_cids.pop(1)
                        self.print(f"* Removed synonym from CID {cid2}")
                    elif decision == 's2':
                        cids_syns_to_del[cid1].add(norm_name)
                        conflict_cids.pop(0)
                        self.print(f"* Removed synonym from CID {cid1}")
                    elif decision == 'c1':
                        cids_to_delete.add(cid2)
                        conflict_cids.pop(1)
                        self.print(f"* Removed compound with CID {cid2}")
                    elif decision == 'c2':
                        cids_to_delete.add(cid1)
                        conflict_cids.pop(0)
                        self.print(f"* Removed compound with CID {cid1}")
                    elif decision == 'exit':
                        stop = True
                        break
                    else:
                        self.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                        continue
                    
                    print()
                
                if stop:
                    break

        finally:
            self._update_cids_blacklist(cids_to_delete)
            self.__update_cids_filtered_synonyms(cids_syns_to_del)

            resolved_chems = []
            for chem in self.chems:
                cid = chem['cid']
                if cid not in cids_to_delete:
                    if cid in cids_syns_to_del and cids_syns_to_del[cid]:
                        syns_to_del = cids_syns_to_del[cid]
                        chem['cmpdsynonym'] = list(filter(lambda x: self._normalize_chem_name(x, is_clean=True) not in syns_to_del, chem['cmpdsynonym']))
                        if not chem['cmpdsynonym']:
                            continue
                        if self._normalize_chem_name(chem['cmpdname'], is_clean=True) in syns_to_del:
                            chem['cmpdname'] = chem['cmpdsynonym'][0]
                    resolved_chems.append(chem)

            self._update_chems(resolved_chems)


if __name__ == "__main__":
    synonyms = ChemsSynonyms('data/')
    synonyms.resolve_conflicting_synonyms()