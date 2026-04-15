import os
import re

import numpy as np

from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED


class ThermoLLMOps:
    def __init__(
        self,
        data_dir,
        reactions,
        store,
        logger,
        properties,
        thermo,
        llm_client,
        reaction_parser,
    ):
        self.reactions = reactions
        self.store = store
        self.logger = logger
        self.properties = properties
        self.thermo = thermo
        self.llm_client = llm_client
        self.reaction_parser = reaction_parser

        self.llm_thermo_dir = os.path.join(data_dir, 'thermo', 'llm')
        os.makedirs(self.llm_thermo_dir, exist_ok=True)
        self.reactions_thermo_llm_fn = os.path.join(self.llm_thermo_dir, 'reactions_thermo_llm.jsonl')
        self.chems_formation_thermo_llm_fn = os.path.join(
            self.llm_thermo_dir,
            'chems_formation_thermo_llm.jsonl',
        )
        self.corrected_reactions_thermo_llm_fn = os.path.join(
            self.llm_thermo_dir,
            'corrected_reactions_thermo_llm.jsonl',
        )

        store.register_sorting(self.reactions_thermo_llm_fn, 'rid')
        store.register_sorting(self.chems_formation_thermo_llm_fn, 'cid')
        store.register_sorting(self.corrected_reactions_thermo_llm_fn, 'rid')

    def _get_reactions_thermo_batch(self, reactions):
        if not reactions:
            return None

        instruct = (
            "You will be given a list of chemical reaction schemes. "
            "Your task is to estimate the enthalpy and free energy of each reaction based on your chemical "
            "knowledge. Assume standard conditions. Provide both values in kcal/mole as plain integers, "
            "separated by a comma, one reaction per line. Format: <enthalpy>, <free energy>\n"
            "Do not include anything other than the estimates."
        )

        model = self.reaction_parser.gpt_oss
        reactions_formatted_str = '\n'.join(
            [f"{i+1}. {self.reactions.get_reaction_as_str(react)}" for i, react in enumerate(reactions)]
        )

        for _ in range(2):
            response = self.llm_client.fetch_answer_str(f"{instruct}\n\n{reactions_formatted_str}", model)
            lines = response.strip().split('\n')
            if len(lines) != len(reactions):
                continue

            results = []
            for i, line in enumerate(lines):
                try:
                    dH, dG = re.sub(r'\s+', '', line).split(',')
                    dH_value = float(dH)
                    dG_value = float(dG)
                except (TypeError, ValueError):
                    continue

                results.append(
                    {
                        'rid': reactions[i]['rid'],
                        'estimates': {'dH': dH_value, 'dG': dG_value, 'source': model},
                    }
                )
            return results

        self.logger.log_err(f"Failed to parse LLM response for batch of {len(reactions)} reactions")
        return None

    def get_reactions_thermo(self, max_workers=1):
        current_thermo = self.store.load_jsonl(self.reactions_thermo_llm_fn)
        current_thermo_map = {entry['rid']: entry['estimates'] for entry in current_thermo}

        def save_results():
            self.logger.log("Writing results...")
            result_thermo = []
            for rid, estimates in current_thermo_map.items():
                dHs = [x['dH'] for x in estimates]
                dGs = [x['dG'] for x in estimates]
                mean = {'dH': np.mean(dHs), 'dG': np.mean(dGs)}
                std = {'dH': np.std(dHs), 'dG': np.std(dGs)}
                result_thermo.append({'rid': rid, 'estimates': estimates, 'mean': mean, 'std': std})

            self.store.write_jsonl(result_thermo, self.reactions_thermo_llm_fn)

        def missing_estimates():
            return [
                reaction
                for reaction in self.properties.parsed_reactions_balanced
                if len(current_thermo_map.get(reaction['rid'], [])) < 10
            ]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = set()
            reactions_staged = missing_estimates()
            iter_i = 0
            try:
                while reactions_staged:
                    iter_i += 1

                    for i in range(0, len(reactions_staged), 10):
                        batch = reactions_staged[i:i + 10]
                        futures.add(executor.submit(self._get_reactions_thermo_batch, batch))

                    with self.logger.progress(transient=True) as progress:
                        task = progress.add_task(f"Iteration {iter_i}", total=len(futures))
                        while futures:
                            done, futures = wait(futures, return_when=FIRST_COMPLETED)
                            for future in done:
                                results = self.llm_client.get_future_result(future, executor)
                                if not results:
                                    continue

                                for entry in results:
                                    current_thermo_map.setdefault(entry['rid'], []).append(entry['estimates'])

                            progress.update(task, advance=len(done))
                            progress.refresh()

                    reactions_staged = missing_estimates()
                    save_results()

            finally:
                save_results()
                self.logger.log("Done!")

    def get_llm_thermo_estimations_distribution(self, num_bins=20):
        import matplotlib.pyplot as plt

        thermo_entries = self.store.load_jsonl(self.reactions_thermo_llm_fn)
        thermo_entries = [entry for entry in thermo_entries if abs(entry['mean']['dH']) < 500]

        means_dH = np.array([entry['mean']['dH'] for entry in thermo_entries])
        stddevs_dH = np.array([entry['std']['dH'] for entry in thermo_entries])
        means_dG = np.array([entry['mean']['dG'] for entry in thermo_entries])
        stddevs_dG = np.array([entry['std']['dG'] for entry in thermo_entries])

        def bin_and_average(means, stddevs, bins_count):
            bins = np.linspace(min(means), max(means), bins_count + 1)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            avg_stddevs = []
            for i in range(bins_count):
                in_bin = (means >= bins[i]) & (means < bins[i + 1])
                avg_stddevs.append(np.mean(stddevs[in_bin]) if np.any(in_bin) else np.nan)
            return bin_centers, avg_stddevs

        bin_centers_dH, avg_stddevs_dH = bin_and_average(means_dH, stddevs_dH, num_bins)
        bin_centers_dG, avg_stddevs_dG = bin_and_average(means_dG, stddevs_dG, num_bins)

        plt.figure()
        plt.plot(bin_centers_dH, avg_stddevs_dH, marker='o', label='dH')
        plt.plot(bin_centers_dG, avg_stddevs_dG, marker='o', label='dG')
        plt.xlabel('Mean value')
        plt.ylabel('Average standard deviation')
        plt.title('Dependence of standard deviation on mean')
        plt.legend()
        plt.grid(True)
        plt.show()

    def filter_anomaly_llm_thermo(self):
        entries = self.store.load_jsonl(self.reactions_thermo_llm_fn)
        entries = [entry for entry in entries if abs(entry['mean']['dH']) < 500 and abs(entry['mean']['dG']) < 100]
        self.store.write_jsonl(entries, self.reactions_thermo_llm_fn)

    def process_llm_thermo_estimates(self):
        from scipy.sparse import lil_matrix
        from scipy.sparse.linalg import lsqr

        thermo_entries = self.store.load_jsonl(self.reactions_thermo_llm_fn)
        thermo_map = {entry['rid']: entry for entry in thermo_entries}
        chem_reactions_occurence = self.properties.get_chems_reactions_occurence(self.properties.parsed_reactions)

        def reactions_filter(reaction):
            if reaction['rid'] not in thermo_map:
                return False
            if reaction['source'] == 'ord':
                return False
            all_cids = self.reactions.get_all_reaction_cids(reaction)
            return not any(chem_reactions_occurence[cid] < 3 for cid in all_cids)

        reactions = list(filter(reactions_filter, self.properties.parsed_reactions_balanced))

        cid_to_index = {}
        for reaction in reactions:
            for cid in self.reactions.get_all_reaction_cids(reaction):
                cid_to_index.setdefault(cid, len(cid_to_index))

        def coeff_map(reaction):
            balance = self.reactions.reactions_balance[reaction['rid']]
            coeffs = {}
            for reagent in reaction['reagents']:
                coeffs[reagent['cid']] = -balance[reagent['cid']]
            for product in reaction['products']:
                coeffs[product['cid']] = balance[product['cid']]
            return coeffs

        def normalize_value(value):
            # Clamp near-zero magnitudes to 1 to avoid inflating weights for small estimates
            return value if abs(value) > 1 else 1

        reaction_coeffs = {reaction['rid']: coeff_map(reaction) for reaction in reactions}

        # Build design matrix M: one row per (reaction, estimate) pair.
        # Each row is scaled by 1/normalize(est) so lsqr minimises the weighted residual directly,
        # avoiding the numerical instability of forming normal equations explicitly.
        rows = [
            (reaction, est)
            for reaction in reactions
            for est in thermo_map[reaction['rid']]['estimates']
        ]

        n_rows, n_cols = len(rows), len(cid_to_index)
        M_dH = lil_matrix((n_rows, n_cols))
        M_dG = lil_matrix((n_rows, n_cols))
        b_dH = np.zeros(n_rows)
        b_dG = np.zeros(n_rows)

        for row_i, (reaction, est) in self.logger.track(enumerate(rows), "Building LSQR matrix", total=n_rows):
            coeffs = reaction_coeffs[reaction['rid']]
            w_dH = 1.0 / normalize_value(abs(est['dH']))
            w_dG = 1.0 / normalize_value(abs(est['dG']))

            for cid, coeff in coeffs.items():
                cid_i = cid_to_index[cid]
                M_dH[row_i, cid_i] = coeff * w_dH
                M_dG[row_i, cid_i] = coeff * w_dG

            b_dH[row_i] = est['dH'] * w_dH
            b_dG[row_i] = est['dG'] * w_dG

        dH = lsqr(M_dH.tocsr(), b_dH)[0]
        dG = lsqr(M_dG.tocsr(), b_dG)[0]
        cid_to_dH = {cid: dH[cid_i] for cid, cid_i in cid_to_index.items()}
        cid_to_dG = {cid: dG[cid_i] for cid, cid_i in cid_to_index.items()}

        formation_thermo = []
        for cid, cid_i in cid_to_index.items():
            dHf = self.thermo.compute_formation_value(cid, dH[cid_i], cid_to_dH)
            dGf = self.thermo.compute_formation_value(cid, dG[cid_i], cid_to_dG)
            formation_thermo.append({'cid': cid, 'dHf': dHf, 'dGf': dGf})

        self.store.write_jsonl(formation_thermo, self.chems_formation_thermo_llm_fn)

        def compute_reaction_thermo(reaction, cid_to_value):
            coeffs = reaction_coeffs[reaction['rid']]
            return sum(coeff * cid_to_value[cid] for cid, coeff in coeffs.items())

        corrected_thermo = [
            {
                'rid': reaction['rid'],
                'dH': compute_reaction_thermo(reaction, cid_to_dH),
                'dG': compute_reaction_thermo(reaction, cid_to_dG),
            }
            for reaction in reactions
        ]

        self.store.write_jsonl(corrected_thermo, self.corrected_reactions_thermo_llm_fn)
