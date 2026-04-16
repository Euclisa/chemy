class ConflictResolver:
    def __init__(self, compounds, logger, store, properties):
        self.compounds = compounds
        self.logger = logger
        self.store = store
        self.properties = properties

    def _build_synonym_filter_decision(self, cid, norm_synonym, conflict_cids):
        chem = self.compounds.cid_chem_map[cid]
        raw_synonyms = [
            syn for syn in chem['cmpdsynonym']
            if self.compounds.normalize_chem_name(syn, is_clean=True) == norm_synonym
        ]
        return {
            'cid': cid,
            'norm_synonym': norm_synonym,
            'raw_synonyms': raw_synonyms,
            'decision': 'remove_synonym',
            'reason': 'manual_conflict',
            'source': 'manual_conflict_resolver',
            'other_cids': [other_cid for other_cid in conflict_cids if other_cid != cid],
        }

    def find_synonym_conflicts(self, only_relevant=True):
        conflict_map = {
            name: cids
            for name, cids in self.compounds.name_cids_map.items()
            if len([cid for cid in cids if cid > 0]) > 1
        }
        if only_relevant:
            relevant_names = self.properties.get_parsed_reactions_participants_norm_names()
            conflict_map = {
                name: cids for name, cids in conflict_map.items() if name in relevant_names
            }

        return conflict_map

    def apply_synonym_decisions(self, decisions, cids_to_delete=None):
        cids_to_delete = set() if cids_to_delete is None else set(cids_to_delete)
        self.compounds.update_cids_blacklist(cids_to_delete)

        decisions = [
            self.compounds._normalize_synonym_filter_decision(decision)
            for decision in decisions
            if decision['cid'] > 0
        ]
        self.compounds.write_synonym_filter_decisions(
            list(self.compounds.synonym_filter_decisions) + decisions
        )

        cids_syns_to_del = dict()
        for decision in decisions:
            if decision['decision'] != 'remove_synonym':
                continue
            cids_syns_to_del.setdefault(decision['cid'], set()).add(decision['norm_synonym'])

        resolved_chems = []
        for chem in self.compounds.chems:
            cid = chem['cid']
            if cid in cids_to_delete:
                continue

            if cid in cids_syns_to_del and cids_syns_to_del[cid]:
                syns_to_del = cids_syns_to_del[cid]
                chem['cmpdsynonym'] = list(
                    filter(
                        lambda value: self.compounds.normalize_chem_name(
                            value,
                            is_clean=True,
                        ) not in syns_to_del,
                        chem['cmpdsynonym'],
                    )
                )
                if not chem['cmpdsynonym'] and self.compounds.is_chem_mapped(chem):
                    self.logger.log(f"Removed compound {cid} (no synonyms left)")
                    continue

                if self.compounds.normalize_chem_name(
                    chem['cmpdname'],
                    is_clean=True,
                ) in syns_to_del:
                    if self.compounds.is_chem_mapped(chem):
                        chem['cmpdname'] = chem['cmpdsynonym'][0]
                    else:
                        chem['cmpdname'] = self.compounds.unknown_name_ph
                        chem['cmpdsynonym'] = []

            resolved_chems.append(chem)

        self.compounds.update_chems(resolved_chems)

    def _display_conflict_table_cid(self, conflict_i, conflict_str, cid1, cid2, display_compound=None):
        chem1 = self.compounds.cid_chem_map[cid1]
        chem2 = self.compounds.cid_chem_map[cid2]
        self._display_conflict_table(
            conflict_i,
            conflict_str,
            chem1,
            chem2,
            display_compound=display_compound,
        )

    def _display_conflict_table(self, conflict_i, conflict_str, chem1, chem2, display_compound=None):
        if display_compound is None:
            display_compound = self.compounds.display_compound_table

        self.compounds.display_compare_table(
            f"[bold yellow]CONFLICT {conflict_i}: '{conflict_str}'[/bold yellow]",
            chem1,
            chem2,
            display_compound=display_compound,
        )

    def resolve_conflicting_synonyms(self, only_relevant=True):
        self.logger.print("Resolve options:")
        self.logger.print("  s0   - discard conflicted synonym from both compounds")
        self.logger.print("  s1   - retain conflicted synonym for 1st compound")
        self.logger.print("  s2   - retain conflicted synonym for 2nd compound")
        self.logger.print("  c0   - discard both compounds")
        self.logger.print("  c1   - retain 1st compound")
        self.logger.print("  c2   - retain 2nd compound")
        self.logger.print("  save - save all previous decisions and abort")
        self.logger.print("  exit - abort without saving")
        self.logger.print()

        def display_conflict(conflict_i, conflict_norm_name, cid1, cid2):
            def get_conflict_inds(chem):
                return [
                    i for i, value in enumerate(chem['cmpdsynonym'])
                    if self.compounds.normalize_chem_name(value, is_clean=True) == conflict_norm_name
                ]

            def get_conflict_synonyms_str(syns, conflict_inds):
                return '; '.join(f"'{syns[i]}' ({i}/{len(syns)})" for i in conflict_inds)

            def display_compound(compound_i, chem, cid):
                syns = chem['cmpdsynonym']
                conflict_inds = get_conflict_inds(chem)
                conflict_synonyms_str = get_conflict_synonyms_str(syns, conflict_inds)
                self.compounds.display_compound_table(
                    compound_i,
                    chem,
                    cid,
                    [("Conflict synonyms", f"[magenta]{conflict_synonyms_str}[/magenta]")],
                )

            self._display_conflict_table_cid(
                conflict_i,
                conflict_norm_name,
                cid1,
                cid2,
                display_compound=display_compound,
            )

        cids_to_delete = set()
        decisions = []
        stop = False
        save = True
        try:
            conflict_map = self.find_synonym_conflicts(only_relevant=only_relevant)

            self.logger.print()
            self.logger.print(f"{len(conflict_map)} conflicting names pending resolution")
            self.logger.print()
            conflict_i = 0
            for norm_name in conflict_map:
                conflict_cids = sorted(
                    [cid for cid in conflict_map[norm_name] if cid not in cids_to_delete]
                )
                if len(conflict_cids) <= 1:
                    continue

                while len(conflict_cids) >= 2:
                    cid1, cid2 = conflict_cids[0], conflict_cids[1]
                    conflict_i += 1

                    if cid1 * cid2 > 0:
                        if cid1 > 0:
                            display_conflict(conflict_i, norm_name, cid1, cid2)
                            decision = input("* Decision: ").strip()
                        else:
                            decision = 's0'
                    else:
                        decision = 'c2'

                    if decision == 's0':
                        decisions.append(
                            self._build_synonym_filter_decision(cid1, norm_name, conflict_cids)
                        )
                        decisions.append(
                            self._build_synonym_filter_decision(cid2, norm_name, conflict_cids)
                        )
                        conflict_cids = conflict_cids[2:]
                        self.logger.print("* Removed synonym from both compounds")
                    elif decision == 's1':
                        decisions.append(
                            self._build_synonym_filter_decision(cid2, norm_name, conflict_cids)
                        )
                        conflict_cids.pop(1)
                        self.logger.print(f"* Removed synonym from CID {cid2}")
                    elif decision == 's2':
                        decisions.append(
                            self._build_synonym_filter_decision(cid1, norm_name, conflict_cids)
                        )
                        conflict_cids.pop(0)
                        self.logger.print(f"* Removed synonym from CID {cid1}")
                    elif decision == 'c0':
                        cids_to_delete.add(cid1)
                        cids_to_delete.add(cid2)
                        conflict_cids = conflict_cids[2:]
                        self.logger.print("* Removed both compounds")
                    elif decision == 'c1':
                        cids_to_delete.add(cid2)
                        conflict_cids.pop(1)
                        self.logger.print(f"* Removed compound with CID {cid2}")
                    elif decision == 'c2':
                        cids_to_delete.add(cid1)
                        conflict_cids.pop(0)
                        self.logger.print(f"* Removed compound with CID {cid1}")
                    elif decision == 'save':
                        stop = True
                        break
                    elif decision == 'exit':
                        stop = True
                        save = False
                        break
                    else:
                        self.logger.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                        continue

                if stop:
                    break

        finally:
            if not save:
                return

            self.apply_synonym_decisions(decisions, cids_to_delete=cids_to_delete)

    def resolve_conflicting_inchi(self):
        self.logger.print("Resolve options:")
        self.logger.print("  m1   - retain first compound with merge")
        self.logger.print("  m2   - retain second compound with merge")
        self.logger.print("  c1   - retain first compound without merge")
        self.logger.print("  c2   - retain second compound without merge")
        self.logger.print("  save - save all previous decisions and abort")
        self.logger.print("  exit - abort without saving")
        self.logger.print()

        conflict_map = dict()
        for chem in self.compounds.chems:
            cid = chem['cid']
            inchi_norm = self.compounds.get_chem_norm_inchi(chem)
            conflict_map.setdefault(inchi_norm, []).append(cid)

        conflict_map = {
            inchi: sorted(cids) for inchi, cids in conflict_map.items() if len(cids) > 1
        }

        self.logger.print(f"{len(conflict_map)} conflicting inchi pending resolution")
        self.logger.print()

        cids_to_delete = set()
        stop = False
        save = True
        try:
            conflict_i = 0
            for conflict_inchi, conflict_cids in conflict_map.items():
                while len(conflict_cids) >= 2:
                    cid1, cid2 = conflict_cids[0], conflict_cids[1]
                    chem1 = self.compounds.cid_chem_map[cid1]
                    chem2 = self.compounds.cid_chem_map[cid2]
                    conflict_i += 1

                    if cid1 * cid2 > 0:
                        if cid1 > 0:
                            decision = 'm2' if len(chem1['inchi']) > len(chem2['inchi']) else 'm1'
                        else:
                            self._display_conflict_table_cid(conflict_i, conflict_inchi, cid1, cid2)
                            decision = input("* Decision: ").strip()
                    else:
                        decision = 'c2'

                    if decision == 'm1':
                        self.compounds.merge_chems(chem1, chem2)
                        conflict_cids.pop(1)
                        cids_to_delete.add(cid2)
                        self.logger.log(f"Discarded compound with CID {cid2}")
                    elif decision == 'm2':
                        self.compounds.merge_chems(chem2, chem1)
                        conflict_cids.pop(0)
                        cids_to_delete.add(cid1)
                        self.logger.log(f"Discarded compound with CID {cid1}")
                    elif decision == 'c1':
                        conflict_cids.pop(1)
                        cids_to_delete.add(cid2)
                        self.logger.log(f"Discarded compound with CID {cid2}")
                    elif decision == 'c2':
                        conflict_cids.pop(0)
                        cids_to_delete.add(cid1)
                        self.logger.log(f"Discarded compound with CID {cid1}")
                    elif decision == 'save':
                        stop = True
                        break
                    elif decision == 'exit':
                        stop = True
                        save = False
                        break
                    else:
                        self.logger.print(f"[red]!!! Invalid decision '{decision}' !!![/red]")
                        continue

                if stop:
                    break

        finally:
            if not save:
                return

            updated_chems = [chem for chem in self.compounds.chems if chem['cid'] not in cids_to_delete]
            self.compounds.update_cids_blacklist(cids_to_delete)
            self.compounds.update_chems(updated_chems)
