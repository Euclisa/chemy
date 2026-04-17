import json
import os
from collections import Counter, defaultdict

from .layout import RawReactionsLayout
from .models import DEFAULT_MODELS
from .services import batched
from scripts.data import chem_names
from scripts.infra.parallel_runner import iter_parallel


DEFAULT_ROUNDS = 5
DEFAULT_MIN_VOTES = 3
DEFAULT_BATCH_SIZE = 7


def build_synonym_prompt(names):
    payload = "\n".join(json.dumps({"name": name}) for name in names)
    return (
        "For each input chemical name, list common equivalent names for the same "
        "compound only. Use names that could identify the compound in a reaction.\n"
        "Return exactly one JSON object per input line, in the same order:\n"
        '{"name":"<input name>","synonyms":["<name1>","<name2>"]}\n'
        "Rules:\n"
        "- If unsure, use an empty synonyms list.\n"
        "- Return only JSONL, no extra text.\n\n"
        f"{payload}"
    )


def parse_synonym_jsonl(response, input_norms, normalize_name):
    results = defaultdict(list)
    input_norms = set(input_norms)

    for raw_line in (response or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
        except (json.JSONDecodeError, TypeError):
            continue
        if not isinstance(obj, dict):
            continue

        norm_name = normalize_name(obj.get("name", ""))
        if norm_name not in input_norms:
            continue

        synonyms = obj.get("synonyms") or []
        if not isinstance(synonyms, list):
            continue
        results[norm_name].extend(syn for syn in synonyms if isinstance(syn, str))

    return results


class RawReactionSynonymEnricher:
    def __init__(self, data_dir, compounds, store, logger, llm_client, parser):
        self.data_dir = data_dir
        self.compounds = compounds
        self.store = store
        self.logger = logger
        self.llm_client = llm_client
        self.parser = parser
        self.layout = RawReactionsLayout(data_dir)
        self.models = DEFAULT_MODELS

        self.output_dir = os.path.join(data_dir, "raw_reactions", "synonym_enrichment")
        self.default_audit_jsonl_fn = os.path.join(self.output_dir, "candidates.jsonl")
        self.default_audit_md_fn = os.path.join(self.output_dir, "candidates.md")

    def _list_all_raw_files(self):
        base = os.path.join(self.data_dir, "raw_reactions")
        if not os.path.isdir(base):
            return []

        raw_files = []
        for dirname in sorted(os.listdir(base)):
            dirpath = os.path.join(base, dirname)
            if not os.path.isdir(dirpath):
                continue
            for fn in sorted(os.listdir(dirpath)):
                if fn.startswith("raw_") and fn.endswith(".jsonl"):
                    raw_files.append(os.path.join(dirpath, fn))
        return raw_files

    def _raw_files(self, preset_name=None):
        if preset_name:
            return self.layout.raw_all(preset_name)
        return self._list_all_raw_files()

    def collect_unmapped_names(self, preset_name=None, limit=None):
        counts = Counter()
        variants = defaultdict(Counter)
        examples = defaultdict(list)

        for raw_fn in self._raw_files(preset_name):
            for entry in self.store.load_jsonl(raw_fn):
                for reaction_obj in entry.get("reactions") or []:
                    _, unmapped = self.parser.parse_structured_reaction(reaction_obj)
                    for norm_name, raw_name in unmapped:
                        clean_name = self.compounds.clean_chem_name(raw_name)
                        counts[norm_name] += 1
                        variants[norm_name][clean_name] += 1
                        if len(examples[norm_name]) < 3:
                            examples[norm_name].append(reaction_obj)

        entries = []
        for norm_name, count in counts.items():
            raw_variants = [
                {"name": name, "count": variant_count}
                for name, variant_count in variants[norm_name].most_common()
            ]
            entries.append({
                "name": raw_variants[0]["name"],
                "norm_name": norm_name,
                "count": count,
                "raw_variants": raw_variants,
                "examples": examples[norm_name],
            })

        entries.sort(key=lambda entry: (-entry["count"], entry["norm_name"]))
        if limit is not None:
            entries = entries[:limit]
        return entries

    def _response_to_vote_maps(self, response, batch):
        input_norms = [entry["norm_name"] for entry in batch]
        parsed = parse_synonym_jsonl(
            response,
            input_norms,
            self.compounds.normalize_chem_name,
        )
        votes = defaultdict(Counter)
        raw_synonyms = defaultdict(lambda: defaultdict(Counter))

        for norm_name, synonyms in parsed.items():
            seen_in_round = set()
            for synonym in synonyms:
                result = chem_names.classify_name(
                    synonym,
                    unknown_name=self.compounds.unknown_name_ph,
                )
                if not result.accepted:
                    continue
                if result.norm_name in seen_in_round:
                    raw_synonyms[norm_name][result.norm_name][result.clean_name] += 1
                    continue
                seen_in_round.add(result.norm_name)
                votes[norm_name][result.norm_name] += 1
                raw_synonyms[norm_name][result.norm_name][result.clean_name] += 1

        return votes, raw_synonyms

    def _resolve_candidate(self, norm_synonym):
        cids = self.compounds.name_cids_map.get(norm_synonym, ())
        positive_cids = tuple(cid for cid in cids if cid > 0)
        if len(positive_cids) == 1:
            cid = positive_cids[0]
            return {
                "decision": "matched",
                "matched_cid": cid,
                "matched_name": self.compounds.cid_chem_map[cid]["cmpdname"],
            }
        if len(positive_cids) > 1:
            return {
                "decision": "ambiguous_db_match",
                "matched_cids": list(positive_cids),
            }
        return {"decision": "no_db_match"}

    def _build_audit_entry(self, entry, votes, raw_synonyms, min_votes):
        candidate_entries = []
        for norm_synonym, vote_count in sorted(
            votes.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            if vote_count < min_votes:
                continue
            raw_counts = raw_synonyms[norm_synonym]
            raw_forms = [
                {"synonym": synonym, "count": count}
                for synonym, count in raw_counts.most_common()
            ]
            candidate = {
                "synonym": raw_forms[0]["synonym"],
                "norm_synonym": norm_synonym,
                "votes": vote_count,
                "raw_synonyms": raw_forms,
            }
            candidate.update(self._resolve_candidate(norm_synonym))
            candidate_entries.append(candidate)

        matched_cids = sorted({
            candidate["matched_cid"]
            for candidate in candidate_entries
            if candidate.get("decision") == "matched"
        })
        add_result = chem_names.classify_name(
            entry["name"],
            unknown_name=self.compounds.unknown_name_ph,
        )

        audit_entry = {
            **entry,
            "candidate_synonyms": candidate_entries,
            "approved": False,
            "target_cid": None,
            "target_name": None,
            "add_synonym": None,
            "status": "no_unique_match",
        }

        if len(matched_cids) > 1:
            audit_entry["status"] = "conflicting_candidate_matches"
        elif len(matched_cids) == 1:
            target_cid = matched_cids[0]
            if not add_result.accepted:
                audit_entry["status"] = f"bad_original_name:{add_result.reason}"
            else:
                guard = self._addition_guard(target_cid, add_result.clean_name)
                if guard is None:
                    audit_entry.update({
                        "approved": True,
                        "target_cid": target_cid,
                        "target_name": self.compounds.cid_chem_map[target_cid]["cmpdname"],
                        "add_synonym": add_result.clean_name,
                        "status": "ready",
                    })
                else:
                    audit_entry["status"] = guard

        return audit_entry

    def _mine_batch(self, batch, model, rounds, min_votes):
        names = [entry["name"] for entry in batch]
        prompt = build_synonym_prompt(names)
        votes_by_name = defaultdict(Counter)
        raw_by_name = defaultdict(lambda: defaultdict(Counter))

        for _ in range(rounds):
            response = self.llm_client.fetch_answer_str(prompt, model)
            round_votes, round_raw = self._response_to_vote_maps(response, batch)
            for norm_name, votes in round_votes.items():
                votes_by_name[norm_name].update(votes)
            for norm_name, by_synonym in round_raw.items():
                for norm_synonym, raw_counts in by_synonym.items():
                    raw_by_name[norm_name][norm_synonym].update(raw_counts)

        return [
            self._build_audit_entry(
                entry,
                votes_by_name[entry["norm_name"]],
                raw_by_name[entry["norm_name"]],
                min_votes,
            )
            for entry in batch
        ]

    def mine(
        self,
        preset_name=None,
        *,
        model=None,
        rounds=DEFAULT_ROUNDS,
        min_votes=DEFAULT_MIN_VOTES,
        batch_size=DEFAULT_BATCH_SIZE,
        max_workers=1,
        limit=None,
        audit_jsonl_fn=None,
        audit_md_fn=None,
    ):
        if rounds < 1:
            raise ValueError("rounds must be at least 1")
        if min_votes < 1:
            raise ValueError("min_votes must be at least 1")
        if batch_size < 1:
            raise ValueError("batch_size must be at least 1")

        model = model or self.models[0]
        audit_jsonl_fn = audit_jsonl_fn or self.default_audit_jsonl_fn
        audit_md_fn = audit_md_fn or self.default_audit_md_fn

        unmapped = self.collect_unmapped_names(preset_name, limit=limit)
        batches = list(batched(unmapped, batch_size))
        audit_entries = []

        for result in iter_parallel(
            self.llm_client,
            batches,
            lambda batch: self._mine_batch(batch, model, rounds, min_votes),
            self.logger,
            max_workers=max_workers,
            description="Mining raw reaction synonym candidates",
        ):
            audit_entries.extend(result)

        audit_entries.sort(key=lambda entry: (-entry["count"], entry["norm_name"]))
        self.write_audit(audit_entries, audit_jsonl_fn, audit_md_fn)
        return {
            "unmapped_names": len(unmapped),
            "audit_entries": len(audit_entries),
            "ready": sum(1 for entry in audit_entries if entry["status"] == "ready"),
            "audit_jsonl": audit_jsonl_fn,
            "audit_markdown": audit_md_fn,
        }

    def write_audit(self, audit_entries, audit_jsonl_fn=None, audit_md_fn=None):
        audit_jsonl_fn = audit_jsonl_fn or self.default_audit_jsonl_fn
        audit_md_fn = audit_md_fn or self.default_audit_md_fn
        for path in (audit_jsonl_fn, audit_md_fn):
            parent = os.path.dirname(path)
            if parent:
                os.makedirs(parent, exist_ok=True)
        self.store.write_jsonl(audit_entries, audit_jsonl_fn, backup=False, no_sort=True)

        with open(audit_md_fn, "w") as f:
            f.write("# Raw Reaction Synonym Candidates\n\n")
            for entry in audit_entries:
                f.write(
                    f"## {entry['name']} "
                    f"({entry['count']} occurrence{'s' if entry['count'] != 1 else ''})\n\n"
                )
                f.write(f"- Status: `{entry['status']}`\n")
                if entry.get("target_cid"):
                    f.write(f"- Target: CID {entry['target_cid']} `{entry['target_name']}`\n")
                    f.write(f"- Add synonym: `{entry['add_synonym']}`\n")
                variants = ", ".join(
                    f"{variant['name']} ({variant['count']})"
                    for variant in entry.get("raw_variants", [])[:5]
                )
                f.write(f"- Variants: {variants}\n")
                if entry.get("candidate_synonyms"):
                    f.write("- Candidate synonyms:\n")
                    for candidate in entry["candidate_synonyms"]:
                        target = ""
                        if candidate.get("matched_cid"):
                            target = (
                                f" -> CID {candidate['matched_cid']} "
                                f"`{candidate['matched_name']}`"
                            )
                        f.write(
                            f"  - `{candidate['synonym']}` "
                            f"({candidate['votes']} votes, {candidate['decision']})"
                            f"{target}\n"
                        )
                f.write("\n")

    def _addition_guard(self, target_cid, synonym):
        target = self.compounds.cid_chem_map.get(target_cid)
        if not target or target_cid <= 0:
            return "invalid_target"
        if len(target.get("cmpdsynonym") or []) >= self.compounds.max_synonyms_thr:
            return "max_synonyms"

        result = chem_names.classify_name(
            synonym,
            manual_filtered=self.compounds.filtered_synonyms_by_cid.get(target_cid, set()),
            unknown_name=self.compounds.unknown_name_ph,
        )
        if not result.accepted:
            return f"bad_name:{result.reason}"

        positive_cids = tuple(
            cid for cid in self.compounds.name_cids_map.get(result.norm_name, ())
            if cid > 0
        )
        if any(cid != target_cid for cid in positive_cids):
            return "would_create_conflict"
        if target_cid in positive_cids:
            return "already_known"
        return None

    def _apply_entry(self, entry, mapped_by_cid, pending_norm_targets):
        if entry.get("approved") is False:
            return {"status": "skipped", "reason": "not_approved", "entry": entry["name"]}
        if entry.get("status") != "ready":
            return {"status": "skipped", "reason": entry.get("status"), "entry": entry["name"]}

        target_cid = entry.get("target_cid")
        synonym = entry.get("add_synonym") or entry.get("name")
        guard = self._addition_guard(target_cid, synonym)
        if guard is not None:
            return {"status": "skipped", "reason": guard, "entry": entry["name"]}

        result = chem_names.classify_name(
            synonym,
            manual_filtered=self.compounds.filtered_synonyms_by_cid.get(target_cid, set()),
            unknown_name=self.compounds.unknown_name_ph,
        )
        pending_target = pending_norm_targets.get(result.norm_name)
        if pending_target is not None and pending_target != target_cid:
            return {"status": "skipped", "reason": "would_create_conflict", "entry": entry["name"]}
        if pending_target == target_cid:
            return {"status": "skipped", "reason": "already_known", "entry": entry["name"]}

        mapped_by_cid[target_cid]["cmpdsynonym"].append(result.clean_name)
        pending_norm_targets[result.norm_name] = target_cid
        return {
            "status": "applied",
            "cid": target_cid,
            "synonym": result.clean_name,
            "entry": entry["name"],
        }

    def apply(self, audit_jsonl_fn=None):
        audit_jsonl_fn = audit_jsonl_fn or self.default_audit_jsonl_fn
        audit_entries = self.store.load_jsonl(audit_jsonl_fn)
        mapped_by_cid = {
            chem["cid"]: {**chem, "cmpdsynonym": list(chem.get("cmpdsynonym") or [])}
            for chem in self.compounds.chems_mapped
        }
        pending_norm_targets = {}
        results = [
            self._apply_entry(entry, mapped_by_cid, pending_norm_targets)
            for entry in audit_entries
        ]

        if any(result["status"] == "applied" for result in results):
            self.compounds.update_mapped_chems(list(mapped_by_cid.values()))

        summary = Counter(result["status"] for result in results)
        skipped = Counter(
            result["reason"]
            for result in results
            if result["status"] == "skipped"
        )
        return {
            "processed": len(results),
            "applied": summary.get("applied", 0),
            "skipped": summary.get("skipped", 0),
            "skipped_by_reason": dict(sorted(skipped.items())),
        }
