import hashlib
import json
import math
import os

import numpy as np

from .thermo_solver import (
    ELEMENT_REFERENCE_SIGMA,
    solve_formation_thermo,
)


class ThermoExperimentOps:
    POLICY_ELEMENTS_ONLY = "elements-only"
    POLICY_BURCAT_ALL = "burcat-all"
    POLICY_BURCAT_SPLIT = "burcat-split"
    POLICY_CORE = "core"

    def __init__(self, data_dir, reactions, store, logger, properties):
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.properties = properties

        self.thermo_dir = os.path.join(data_dir, "thermo")
        self.evidence_dir = os.path.join(self.thermo_dir, "evidence")
        self.solves_dir = os.path.join(self.thermo_dir, "solves")
        self.reports_dir = os.path.join(self.thermo_dir, "reports")
        self.experiments_dir = os.path.join(self.thermo_dir, "experiments")
        os.makedirs(self.evidence_dir, exist_ok=True)
        os.makedirs(self.solves_dir, exist_ok=True)
        os.makedirs(self.reports_dir, exist_ok=True)
        os.makedirs(self.experiments_dir, exist_ok=True)

        self.reaction_estimates_fn = os.path.join(self.evidence_dir, "reaction_estimates.jsonl")
        self.formation_references_fn = os.path.join(self.evidence_dir, "formation_references.jsonl")

        self.llm_solve_dir = os.path.join(self.solves_dir, "llm")
        self.hybrid_solve_dir = os.path.join(self.solves_dir, "hybrid")
        os.makedirs(self.llm_solve_dir, exist_ok=True)
        os.makedirs(self.hybrid_solve_dir, exist_ok=True)

    def _element_formation_references(self):
        references = []
        for _symbol, el_entry in self.properties.compounds.symb_to_el.items():
            references.append(
                {
                    "cid": el_entry["cid"],
                    "source": "element",
                    "reference_type": "element",
                    "dHf": 0.0,
                    "dGf": 0.0,
                    "T": 298.15,
                    "sigma_dHf": ELEMENT_REFERENCE_SIGMA,
                    "sigma_dGf": ELEMENT_REFERENCE_SIGMA,
                    "match": "element",
                }
            )

        return references

    def _load_reference_groups(self):
        references = self.store.load_jsonl(self.formation_references_fn)
        elements = [entry for entry in references if entry.get("reference_type") == "element"]
        burcat = [entry for entry in references if entry.get("reference_type") != "element"]
        if not elements:
            elements = self._element_formation_references()

        elements = sorted(elements, key=lambda entry: entry["cid"])
        burcat = sorted(burcat, key=lambda entry: entry["cid"])
        return elements, burcat

    def _split_key(self, cid, seed):
        payload = f"{seed}:{cid}".encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def select_formation_references(self, policy, train_fraction=0.8, seed=42):
        elements, burcat = self._load_reference_groups()

        if policy == self.POLICY_ELEMENTS_ONLY:
            return elements, burcat

        if policy == self.POLICY_BURCAT_ALL:
            return elements + burcat, []

        if policy == self.POLICY_BURCAT_SPLIT:
            if not 0 <= train_fraction <= 1:
                raise ValueError("train_fraction must be between 0 and 1")

            ordered = sorted(burcat, key=lambda entry: self._split_key(entry["cid"], seed))
            train_count = int(len(ordered) * train_fraction)
            train = ordered[:train_count]
            test = ordered[train_count:]
            return elements + sorted(train, key=lambda entry: entry["cid"]), sorted(test, key=lambda entry: entry["cid"])

        raise ValueError(f"Unknown thermo reference policy: {policy}")

    def select_solve_reactions(self, reaction_estimates):
        estimate_rids = {entry["rid"] for entry in reaction_estimates}
        chem_reactions_occurence = self.properties.get_chems_reactions_occurence(self.properties.parsed_reactions)

        def reactions_filter(reaction):
            if reaction["rid"] not in estimate_rids:
                return False
            if reaction.get("source") == "ord":
                return False
            all_cids = self.reactions.get_all_reaction_cids(reaction)
            return not any(chem_reactions_occurence.get(cid, 0) < 3 for cid in all_cids)

        return list(filter(reactions_filter, self.properties.parsed_reactions_balanced))

    def _load_reaction_estimates(self, reaction_estimates=None):
        if reaction_estimates is not None:
            return reaction_estimates
        return self.store.load_jsonl(self.reaction_estimates_fn)

    def _solve(self, reaction_estimates, formation_references):
        reactions = self.select_solve_reactions(reaction_estimates)
        reaction_rids = {reaction["rid"] for reaction in reactions}
        reaction_estimates = [
            estimate for estimate in reaction_estimates
            if estimate["rid"] in reaction_rids
        ]
        result = solve_formation_thermo(
            reactions,
            self.reactions.reactions_balance,
            reaction_estimates,
            formation_references,
        )
        return result, reactions, reaction_estimates

    def _legacy_paths(self, mode):
        if mode == "llm":
            solve_dir = self.llm_solve_dir
        elif mode == "hybrid":
            solve_dir = self.hybrid_solve_dir
        else:
            raise ValueError(f"Unknown thermo solve mode: {mode}")

        return (
            os.path.join(solve_dir, "chems_formation.jsonl"),
            os.path.join(solve_dir, "reactions_corrected.jsonl"),
            os.path.join(solve_dir, "metadata.json"),
        )

    def run_legacy_solve(self, mode, reaction_estimates=None):
        policy = self.POLICY_ELEMENTS_ONLY if mode == "llm" else self.POLICY_BURCAT_ALL
        reaction_estimates = self._load_reaction_estimates(reaction_estimates)
        formation_references, _test_references = self.select_formation_references(policy)
        result, _reactions, _estimates = self._solve(reaction_estimates, formation_references)

        chems_fn, reactions_fn, metadata_fn = self._legacy_paths(mode)
        self.store.write_jsonl(result.formation, chems_fn)
        self.store.write_jsonl(result.corrected_reactions, reactions_fn)
        with open(metadata_fn, "w") as f:
            json.dump(result.metadata, f, indent=2)

        self.logger.log(
            f"Wrote {mode} thermo solve: {len(result.formation)} chems, "
            f"{len(result.corrected_reactions)} reactions"
        )
        return result

    def validate_legacy_solve(self, mode):
        chems_fn, _reactions_fn, _metadata_fn = self._legacy_paths(mode)
        solved = {entry["cid"]: entry for entry in self.store.load_jsonl(chems_fn)}
        references = [
            entry for entry in self.store.load_jsonl(self.formation_references_fn)
            if entry.get("reference_type") != "element"
        ]

        rows = self._comparison_rows(solved.values(), [], references)
        for row in rows:
            row.pop("subset", None)

        out_fn = os.path.join(self.reports_dir, f"{mode}_burcat_validation.jsonl")
        self.store.write_jsonl(rows, out_fn, no_sort=True)
        return rows

    def _experiment_name(self, policy, train_fraction, seed):
        if policy == self.POLICY_ELEMENTS_ONLY:
            return "elements_only"
        if policy == self.POLICY_BURCAT_ALL:
            return "burcat_all"
        if policy == self.POLICY_BURCAT_SPLIT:
            train_pct = int(round(train_fraction * 100))
            test_pct = int(round((1 - train_fraction) * 100))
            return f"burcat_split_{train_pct}_{test_pct}_seed{seed}"

        raise ValueError(f"Unknown thermo reference policy: {policy}")

    def _write_json(self, payload, filename):
        with open(filename, "w") as f:
            json.dump(payload, f, indent=2)

    def _comparison_rows(self, solved, train_references, test_references):
        solved_by_cid = {entry["cid"]: entry for entry in solved}
        rows = []
        for subset, references in (("train", train_references), ("test", test_references)):
            for ref in references:
                if ref.get("reference_type") == "element":
                    continue

                solved_entry = solved_by_cid.get(ref["cid"], {})
                row = {
                    "cid": ref["cid"],
                    "subset": subset,
                    "source": ref.get("source"),
                    "match": ref.get("match"),
                    "burcat_species": ref.get("burcat_species"),
                    "dHf_ref": ref.get("dHf"),
                    "dHf_solved": solved_entry.get("dHf"),
                    "dHf_error": None,
                    "dGf_ref": ref.get("dGf"),
                    "dGf_solved": solved_entry.get("dGf"),
                    "dGf_error": None,
                }
                if is_finite_number(row["dHf_ref"]) and is_finite_number(row["dHf_solved"]):
                    row["dHf_error"] = row["dHf_solved"] - row["dHf_ref"]
                if is_finite_number(row["dGf_ref"]) and is_finite_number(row["dGf_solved"]):
                    row["dGf_error"] = row["dGf_solved"] - row["dGf_ref"]
                rows.append(row)

        return rows

    def _prop_metrics(self, rows, prop):
        ref_key = f"{prop}_ref"
        error_key = f"{prop}_error"
        reference_count = sum(1 for row in rows if is_finite_number(row.get(ref_key)))
        errors = [
            row[error_key] for row in rows
            if is_finite_number(row.get(ref_key)) and is_finite_number(row.get(error_key))
        ]

        if not errors:
            return {
                "reference_count": reference_count,
                "count": 0,
                "missing_count": reference_count,
                "bias": None,
                "mae": None,
                "rmse": None,
                "median_abs_error": None,
                "p90_abs_error": None,
            }

        abs_errors = [abs(error) for error in errors]
        return {
            "reference_count": reference_count,
            "count": len(errors),
            "missing_count": reference_count - len(errors),
            "bias": float(np.mean(errors)),
            "mae": float(np.mean(abs_errors)),
            "rmse": float(math.sqrt(np.mean([error * error for error in errors]))),
            "median_abs_error": float(np.median(abs_errors)),
            "p90_abs_error": float(np.percentile(abs_errors, 90)),
        }

    def validation_metrics(self, rows):
        metrics = {}
        for subset in ("train", "test", "overall"):
            subset_rows = rows if subset == "overall" else [row for row in rows if row["subset"] == subset]
            metrics[subset] = {
                "dHf": self._prop_metrics(subset_rows, "dHf"),
                "dGf": self._prop_metrics(subset_rows, "dGf"),
            }

        return metrics

    def run_policy(self, policy, train_fraction=0.8, seed=42, name=None, reaction_estimates=None):
        if policy == self.POLICY_CORE:
            return self.run_core_policies(train_fraction=train_fraction, seed=seed, reaction_estimates=reaction_estimates)

        reaction_estimates = self._load_reaction_estimates(reaction_estimates)
        train_references, test_references = self.select_formation_references(policy, train_fraction, seed)
        result, reactions, filtered_estimates = self._solve(reaction_estimates, train_references)

        experiment_name = name or self._experiment_name(policy, train_fraction, seed)
        experiment_dir = os.path.join(self.experiments_dir, experiment_name)
        os.makedirs(experiment_dir, exist_ok=True)

        validation_rows = self._comparison_rows(result.formation, train_references, test_references)
        metrics = self.validation_metrics(validation_rows)
        config = {
            "name": experiment_name,
            "policy": policy,
            "train_fraction": train_fraction if policy == self.POLICY_BURCAT_SPLIT else None,
            "seed": seed if policy == self.POLICY_BURCAT_SPLIT else None,
            "reaction_estimates": len(filtered_estimates),
            "reactions": len(reactions),
            "formation_references_train": len(train_references),
            "formation_references_test": len(test_references),
        }

        self._write_json(config, os.path.join(experiment_dir, "config.json"))
        self.store.write_jsonl(train_references, os.path.join(experiment_dir, "formation_references_train.jsonl"), no_sort=True)
        self.store.write_jsonl(test_references, os.path.join(experiment_dir, "formation_references_test.jsonl"), no_sort=True)
        self.store.write_jsonl(result.formation, os.path.join(experiment_dir, "chems_formation.jsonl"), no_sort=True)
        self.store.write_jsonl(result.corrected_reactions, os.path.join(experiment_dir, "reactions_corrected.jsonl"), no_sort=True)
        self.store.write_jsonl(validation_rows, os.path.join(experiment_dir, "validation.jsonl"), no_sort=True)
        self._write_json({**result.metadata, "experiment": config}, os.path.join(experiment_dir, "metadata.json"))
        self._write_json(metrics, os.path.join(experiment_dir, "metrics.json"))

        self.logger.log(
            f"Wrote thermo experiment '{experiment_name}': "
            f"{len(result.formation)} chems, {len(result.corrected_reactions)} reactions"
        )
        return {
            "name": experiment_name,
            "dir": experiment_dir,
            "config": config,
            "result": result,
            "validation": validation_rows,
            "metrics": metrics,
        }

    def run_core_policies(self, train_fraction=0.8, seed=42, reaction_estimates=None):
        reaction_estimates = self._load_reaction_estimates(reaction_estimates)
        return [
            self.run_policy(self.POLICY_ELEMENTS_ONLY, train_fraction, seed, reaction_estimates=reaction_estimates),
            self.run_policy(self.POLICY_BURCAT_ALL, train_fraction, seed, reaction_estimates=reaction_estimates),
            self.run_policy(self.POLICY_BURCAT_SPLIT, train_fraction, seed, reaction_estimates=reaction_estimates),
        ]


def is_finite_number(value):
    return isinstance(value, (int, float)) and math.isfinite(value)
