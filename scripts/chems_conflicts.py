import re

from rich.table import Table
from rich.rule import Rule

from chems_reaction_properties import ChemsReactionProperties


class ChemsConflicts(ChemsReactionProperties):

    def __init__(self, data_dir):
        super().__init__(data_dir)


    def __update_cids_filtered_synonyms(self, new_entries):
        cids_filtered_synonyms = dict(self.cids_filtered_synonyms)
        for cid, synonyms in new_entries.items():
            if cid < 0 or not synonyms:
                continue
            cids_filtered_synonyms.setdefault(cid, set()).update(synonyms)
        
        res_entries = [{'cid': cid, 'synonyms': sorted(list(syns))} for cid, syns in cids_filtered_synonyms.items()]

        self._write_jsonl(res_entries, self.cids_filtered_synonyms_fn)
        self._clear_cached_property('cids_filtered_synonyms')

    

    def _display_conflict_table_cid(self, conflict_i, conflict_str, cid1, cid2, display_compound=None):
        chem1, chem2 = self.cid_chem_map[cid1], self.cid_chem_map[cid2]
        self._display_conflict_table(conflict_i, conflict_str, chem1, chem2, display_compound=display_compound)
    

    def _display_conflict_table(self, conflict_i, conflict_str, chem1, chem2, display_compound=None):
        self._display_compare_table(f"[bold yellow]CONFLICT {conflict_i}: '{conflict_str}'[/bold yellow]", chem1, chem2, display_compound=display_compound)


    def resolve_conflicting_synonyms(self, only_relevant=True):

        self.print("Resolve options:")
        self.print("  s0   - discard conflicted synonym from both compounds")
        self.print("  s1   - retain conflicted synonym for 1st compound")
        self.print("  s2   - retain conflicted synonym for 2nd compound")
        self.print("  c0   - discard both compounds")
        self.print("  c1   - retain 1st compound")
        self.print("  c2   - retain 2nd compound")
        self.print("  save - save all previous decisions and abort")
        self.print("  exit - abort without saving")
        self.print()

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
                conflict_synonyms_str = get_conflict_synonyms_str(syns, conflict_inds)

                self._display_compound_table(compound_i, chem, cid, [("Conflict synonyms", f"[magenta]{conflict_synonyms_str}[/magenta]")])
            
            self._display_conflict_table_cid(conflict_i, conflict_norm_name, cid1, cid2, display_compound=display_compound)

        conflict_map = dict()
        cids_to_delete = set()
        cids_syns_to_del = dict()
        stop = False
        save = True
        try:
            for chem in self.chems:
                cid = chem['cid']
                cids_syns_to_del[cid] = set()
                norm_syns = set(self._normalize_chem_name(syn, is_clean=True) for syn in chem['cmpdsynonym'])
                for norm_name in norm_syns:
                    conflict_map.setdefault(norm_name, []).append(cid)
            
            conflict_map = {name: cids for name, cids in conflict_map.items() if len(cids) > 1}
            if only_relevant:
                relevant_names = self._get_parsed_reactions_participants_norm_names()
                conflict_map = {name: cids for name, cids in conflict_map.items() if name in relevant_names}

            self.print()
            self.print(f"{len(conflict_map)} conflicting names pending resolution")
            self.print()
            conflict_i = 0
            for norm_name in conflict_map:
                conflict_cids = sorted([cid for cid in conflict_map[norm_name] if cid not in cids_to_delete])
                if len(conflict_cids) <= 1:
                    continue

                while len(conflict_cids) >= 2:
                    cid1, cid2 = conflict_cids[0], conflict_cids[1]

                    conflict_i += 1

                    if cid1*cid2 > 0:
                        if cid1 > 0:
                            display_conflict(conflict_i, norm_name, cid1, cid2)
                            decision = input("* Decision: ").strip()
                        else:
                            decision = 's0'
                    else:
                        decision = 'c2'

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
                    elif decision == 'c0':
                        cids_to_delete.add(cid1)
                        cids_to_delete.add(cid2)
                        conflict_cids = conflict_cids[2:]
                        self.print(f"* Removed both compounds")
                    elif decision == 'c1':
                        cids_to_delete.add(cid2)
                        conflict_cids.pop(1)
                        self.print(f"* Removed compound with CID {cid2}")
                    elif decision == 'c2':
                        cids_to_delete.add(cid1)
                        conflict_cids.pop(0)
                        self.print(f"* Removed compound with CID {cid1}")
                    elif decision == 'save':
                        stop = True
                        break
                    elif decision == 'exit':
                        stop = True
                        save = False
                        break
                    else:
                        self.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                        continue
                
                if stop:
                    break

        finally:
            if not save:
                return

            self._update_cids_blacklist(cids_to_delete)
            self.__update_cids_filtered_synonyms(cids_syns_to_del)

            resolved_chems = []
            for chem in self.chems:
                cid = chem['cid']
                if cid not in cids_to_delete:
                    if cid in cids_syns_to_del and cids_syns_to_del[cid]:
                        syns_to_del = cids_syns_to_del[cid]
                        chem['cmpdsynonym'] = list(filter(lambda x: self._normalize_chem_name(x, is_clean=True) not in syns_to_del, chem['cmpdsynonym']))
                        if not chem['cmpdsynonym'] and self._is_chem_mapped(chem):
                            self.log(f"Removed compound {cid} (no synonyms left)")
                            continue

                        if self._normalize_chem_name(chem['cmpdname'], is_clean=True) in syns_to_del:
                            if self._is_chem_mapped(chem):
                                chem['cmpdname'] = chem['cmpdsynonym'][0]
                            else:
                                chem['cmpdname'] = self.unknown_name_ph
                                chem['cmpdsynonym'] = []
                    resolved_chems.append(chem)

            self._update_chems(resolved_chems)


    def resolve_conflicting_inchi(self):

        self.print("Resolve options:")
        self.print("  m1   - retain first compound with merge")
        self.print("  m2   - retain second compound with merge")
        self.print("  c1   - retain first compound without merge")
        self.print("  c2   - retain second compound without merge")
        self.print("  save - save all previous decisions and abort")
        self.print("  exit - abort without saving")
        self.print()

        conflict_map = dict()
        for chem in self.chems:
            cid = chem['cid']
            inchi_norm = self._get_chem_norm_inchi(chem)    
            conflict_map.setdefault(inchi_norm, []).append(cid)
        
        conflict_map = {inchi: sorted(cids) for inchi, cids in conflict_map.items() if len(cids) > 1}

        self.print(f"{len(conflict_map)} conflicting inchi pending resolution")
        self.print()

        cids_to_delete = set()
        stop = False
        save = True
        try:
            conflict_i = 0
            for conflict_inchi, conflict_cids in conflict_map.items():
                while len(conflict_cids) >= 2:
                    cid1, cid2 = conflict_cids[0], conflict_cids[1]

                    chem1, chem2 = self.cid_chem_map[cid1], self.cid_chem_map[cid2]
                    conflict_i += 1
                    
                    # if one of cids < 0 then delete the one (it must be cid1)
                    if cid1*cid2 > 0:
                        if cid1 > 0:
                            if len(chem1['inchi']) > len(chem2['inchi']):
                                decision = 'm2'
                            else:
                                decision = 'm1'
                        else:
                            self._display_conflict_table_cid(conflict_i, conflict_inchi, cid1, cid2)
                            decision = input("* Decision: ").strip()
                    else:
                        decision = 'c2'

                    if decision == 'm1':
                        self._merge_chems(chem1, chem2)
                        conflict_cids.pop(1)
                        cids_to_delete.add(cid2)
                        self.log(f"Discarded compound with CID {cid2}")
                    elif decision == 'm2':
                        self._merge_chems(chem2, chem1)
                        conflict_cids.pop(0)
                        cids_to_delete.add(cid1)
                        self.log(f"Discarded compound with CID {cid1}")
                    elif decision == 'c1':
                        conflict_cids.pop(1)
                        cids_to_delete.add(cid2)
                        self.log(f"Discarded compound with CID {cid2}")
                    elif decision == 'c2':
                        conflict_cids.pop(0)
                        cids_to_delete.add(cid1)
                        self.log(f"Discarded compound with CID {cid1}")
                    elif decision == 'save':
                        stop = True
                        break
                    elif decision == 'exit':
                        stop = True
                        save = False
                        break
                    else:
                        self.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                        continue
                    
                if stop:
                    break
        
        finally:
            if not save:
                return

            updated_chems = [chem for chem in self.chems if chem['cid'] not in cids_to_delete]

            self._update_cids_blacklist(cids_to_delete)
            self._update_chems(updated_chems)




    


if __name__ == "__main__":
    synonyms = ChemsConflicts('data/')
    synonyms.resolve_conflicting_inchi()