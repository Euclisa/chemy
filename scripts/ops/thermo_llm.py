import os
import re

import numpy as np

from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

from .thermo_solver import llm_sigma


class ThermoLLMOps:
    THERMO_PROMPT_VERSION = 'thermo_v2'

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

        self.thermo_dir = os.path.join(data_dir, 'thermo')
        self.evidence_dir = os.path.join(self.thermo_dir, 'evidence')
        os.makedirs(self.evidence_dir, exist_ok=True)

        self.llm_thermo_dir = os.path.join(data_dir, 'thermo', 'llm')
        os.makedirs(self.llm_thermo_dir, exist_ok=True)
        self.reactions_thermo_llm_fn = os.path.join(self.llm_thermo_dir, 'reactions_thermo_llm.jsonl')

        self.reaction_estimates_fn = os.path.join(self.evidence_dir, 'reaction_estimates.jsonl')

        store.register_sorting(self.reactions_thermo_llm_fn, 'rid')
        store.register_sorting(self.reaction_estimates_fn, 'rid')

    def _parse_thermo_line(self, line):
        number = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
        match = re.fullmatch(
            rf'\s*(?:\d+[\.\)]\s*)?({number})\s*,\s*({number})\s*',
            line,
        )
        if not match:
            return None

        return float(match.group(1)), float(match.group(2))

    def _get_reactions_thermo_batch(self, reactions):
        if not reactions:
            return None

        instruct = (
            "You will be given a list of balanced chemical reaction schemes. "
            "Estimate the standard reaction enthalpy dHr°298 and standard reaction Gibbs free energy dGr°298 "
            "for each reaction exactly as written. "
            "Use 298.15 K and 1 bar. Assume pure substances in their usual standard states at 298.15 K; "
            "if a phase is not shown, infer the usual standard state. "
            "Report values in kcal per balanced reaction equation, not per mole of a selected reagent or product. "
            "Prefer known standard formation thermochemistry when available; otherwise make your best chemically "
            "reasonable estimate. "
            "Return exactly one line per reaction, in the same order. "
            "Each line must contain only two numbers separated by a comma: <dH>, <dG>. "
            "Floats are allowed. Do not include units, labels, numbering, explanations, or extra text."
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
                parsed = self._parse_thermo_line(line)
                if parsed is None:
                    results = None
                    break

                dH_value, dG_value = parsed

                results.append(
                    {
                        'rid': reactions[i]['rid'],
                        'estimates': {
                            'dH': dH_value,
                            'dG': dG_value,
                            'source': model,
                            'prompt': self.THERMO_PROMPT_VERSION,
                        },
                    }
                )

            if results is not None and len(results) == len(reactions):
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
        entries = [entry for entry in entries if abs(entry['mean']['dH']) < 500 and abs(entry['mean']['dG']) < 500]
        self.store.write_jsonl(entries, self.reactions_thermo_llm_fn)

    def normalize_reaction_estimates(self):
        thermo_entries = self.store.load_jsonl(self.reactions_thermo_llm_fn)

        estimates = []
        for entry in thermo_entries:
            samples = entry.get('estimates') or []
            if not samples:
                continue

            dHs = np.array([sample['dH'] for sample in samples], dtype=float)
            dGs = np.array([sample['dG'] for sample in samples], dtype=float)
            mean_dH = float(entry.get('mean', {}).get('dH', np.mean(dHs)))
            mean_dG = float(entry.get('mean', {}).get('dG', np.mean(dGs)))
            std_dH = float(entry.get('std', {}).get('dH', np.std(dHs)))
            std_dG = float(entry.get('std', {}).get('dG', np.std(dGs)))
            sources = sorted(set(sample.get('source') for sample in samples if sample.get('source')))

            estimates.append(
                {
                    'rid': entry['rid'],
                    'source': 'llm',
                    'model_sources': sources,
                    'dH': mean_dH,
                    'dG': mean_dG,
                    'std_dH': std_dH,
                    'std_dG': std_dG,
                    'sigma_dH': llm_sigma(mean_dH, std_dH),
                    'sigma_dG': llm_sigma(mean_dG, std_dG),
                    'sample_count': len(samples),
                }
            )

        self.store.write_jsonl(estimates, self.reaction_estimates_fn)
        return estimates
