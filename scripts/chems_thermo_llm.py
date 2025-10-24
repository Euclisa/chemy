import re
import os
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

from chems_thermo import ChemsThermo
from chems_llm_fetch import ChemsLLMFetch


class ChemsThermoLLM(ChemsThermo, ChemsLLMFetch):
    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.reactions_thermo_llm_fn = os.path.join(self.thermo_dir, 'llm', 'reactions_thermo_llm.jsonl')

        self._file_sorting_prefs[self.reactions_thermo_llm_fn] = 'rid'
    

    def __get_reactions_thermo(self, reactions):
        if not reactions:
            return None

        instruct = (
            "You will be given a list of chemical reaction schemes. "
            f"Your task is to estimate the enthalpy and free energy of each reaction based on your chemical knowledge. "
            "Assume standrd conditions. Provide both values in kcal/mole as plain integers, separated by a comma, one reaction per line. "
            "Format: <enthalpy>, <free energy>\n"
            "Do not include anything other than the estimates."
        )

        model = self.gpt_oss

        reactions_formatted_str = '\n'.join([f"{i+1}. {self._get_reaction_as_str(react)}" for i, react in enumerate(reactions)])

        tries_num = 2
        results = []
        for try_i in range(tries_num):            
            response = self._fetch_llm_answer_str( f"{instruct}\n\n{reactions_formatted_str}", model)
            response = response.strip().split('\n')
            if len(response) == len(reactions):
                def is_float(value):
                    try:
                        float(value)
                        return True
                    except (TypeError, ValueError):
                        return False

                for i, entry in enumerate(response):
                    dH, dG = re.sub(r'\s+', '', entry).split(',')
                    if is_float(dH) and is_float(dG):
                        results.append({'rid': reactions[i]['rid'], 'estimates': {'dH': float(dH), 'dG': float(dG), 'source': model}})
            
                return results

    

    def get_reactions_thermo(self, max_workers=1):
        current_thermo = self._load_jsonl(self.reactions_thermo_llm_fn)

        current_thermo_map = {x['rid']: x['estimates'] for x in current_thermo}

        reactions = [
            r for r in self._load_jsonl(self.reactions_parsed_fn)
            if r['balanced']
        ]

        TARGET_ESTIMATES_NUM = 10
        REACTIONS_BATCH_SIZE = 10

        def save_results():
            self.log(f"Writing results...")
            result_thermo = []
            for rid, est in current_thermo_map.items():
                dHs, dGs = [x['dH'] for x in est], [x['dG'] for x in est]
                mean = {'dH': np.mean(dHs), 'dG': np.mean(dGs)}
                std = {'dH': np.std(dHs), 'dG': np.std(dGs)}

                result_thermo.append({'rid': rid, 'estimates': est, 'mean': mean, 'std': std})

            self._write_jsonl(result_thermo, self.reactions_thermo_llm_fn)

        def missing_estimates():
            return [r for r in reactions if len(current_thermo_map.get(r['rid'], [])) < TARGET_ESTIMATES_NUM]

        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = set()
            reactions_staged = missing_estimates()
            iter_i = 0
            try:
                while reactions_staged:
                    iter_i += 1

                    futures.clear()

                    for i in range(0, len(reactions_staged), REACTIONS_BATCH_SIZE):
                        batch = reactions_staged[i:i+REACTIONS_BATCH_SIZE]
                        futures.add(executor.submit(self.__get_reactions_thermo, batch))

                    with self._rich_progress(transient=True) as progress:
                        task = progress.add_task(f"Iteration {iter_i}", total=len(futures))
                        while futures:
                            done, futures = wait(futures, return_when=FIRST_COMPLETED)
                            processed_cnt = 0
                            for future in done:
                                results = self._get_future_result(future, executor)
                                if not results:
                                    continue

                                for entry in results:
                                    current_thermo_map.setdefault(entry['rid'], []).append(entry['estimates'])
                                    processed_cnt += 1

                            progress.update(task, advance=len(done))
                            progress.refresh()
                    
                    reactions_staged = missing_estimates()
                    save_results()
    
            finally:
                save_results()
                executor.shutdown(wait=False, cancel_futures=True)
                self.log(f"Done!")


    
    def get_llm_thermo_estimations_distribuition(self, num_bins=20):
        thermo_entries = self._load_jsonl(self.reactions_thermo_llm_fn)
        thermo_entries = list(filter(lambda x: abs(x['mean']['dH']) < 500, thermo_entries))

        means_dH, stddevs_dH = [], []
        means_dG, stddevs_dG = [], []

        for entry in thermo_entries:
            means_dH.append(entry['mean']['dH'])
            stddevs_dH.append(entry['std']['dH'])
            means_dG.append(entry['mean']['dG'])
            stddevs_dG.append(entry['std']['dG'])

        means_dH = np.array(means_dH)
        stddevs_dH = np.array(stddevs_dH)
        means_dG = np.array(means_dG)
        stddevs_dG = np.array(stddevs_dG)

        def bin_and_average(means, stddevs, num_bins):
            bins = np.linspace(min(means), max(means), num_bins + 1)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            avg_stddevs = []

            for i in range(num_bins):
                in_bin = (means >= bins[i]) & (means < bins[i+1])
                if np.any(in_bin):
                    avg_stddevs.append(np.mean(stddevs[in_bin]))
                else:
                    avg_stddevs.append(np.nan)  # no points in this bin

            return bin_centers, avg_stddevs

        means_dH = np.array(means_dH)
        stddevs_dH = np.array(stddevs_dH)
        means_dG = np.array(means_dG)
        stddevs_dG = np.array(stddevs_dG)

        bin_centers_dH, avg_stddevs_dH = bin_and_average(means_dH, stddevs_dH, num_bins)
        bin_centers_dG, avg_stddevs_dG = bin_and_average(means_dG, stddevs_dG, num_bins)

        # Plot
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
        entries = self._load_jsonl(self.reactions_thermo_llm_fn)
        entries = list(filter(lambda x: abs(x['mean']['dH']) < 500 and abs(x['mean']['dG']) < 500, entries))
        self._write_jsonl(entries, self.reactions_thermo_llm_fn)
    


if __name__ == "__main__":
    thermo_llm = ChemsThermoLLM('data/')
    thermo_llm.filt()
