from openai import OpenAI
import threading
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import re
import random

from chems_llm import ChemsLLM



class ChemsLLMFetch(ChemsLLM):
    def __init__(self, data_dir, api_key=None):
        super().__init__(data_dir)

        if api_key is None:
            api_key = os.getenv("OPENROUTER_API_KEY")

        self.client = OpenAI(base_url="https://openrouter.ai/api/v1", api_key=api_key)

        self.completion_tokens_total = 0
        self.input_tokens_total = 0
        self.tokens_total_lock = threading.Lock()

    
    def __fetch_llm_answer(self, messages, model, reasoning_effort="medium"):
        completion = self.client.chat.completions.create(
            model=model,
            messages=messages,
            reasoning_effort=reasoning_effort
            )

        with self.tokens_total_lock:
            self.input_tokens_total += completion.usage.prompt_tokens
            self.completion_tokens_total += completion.usage.completion_tokens

        return completion.choices[0].message.content

    def __get_processed_entries(self, out_fn, key):
        return set(x[key] for x in self._load_jsonl(out_fn))


    def __get_reactions_from_response(self, response: str):
        reactions = []
        for line in response.split('\n'):
            if '->' in line:
                reactions.append(line)
        
        return reactions
    
    
    def __fetch_raw_reactions(self, chem, mode="documented_rp"):
        try:
            chem_name = chem['cmpdname']
            
            if mode == "documented_rp":
                compound_role="a reactant or a product"
                special_requirement = ""

            elif mode == "documented_less_common_rp":
                compound_role="a reactant or a product"
                special_requirement = "Include not only the most common reactions, but also less common or more unusual or exotic ones, as long as you are absolutely sure they are real and correct."

            elif mode == "documented_less_common_p":
                compound_role="a product"
                special_requirement = "Include not only the most common reactions, but also less common or more unusual or exotic ones, as long as you are absolutely sure they are real and correct."

            elif mode == "rare_rp":
                compound_role="a reactant or a product"
                special_requirement = "Include only the uncommon, rare and exotic reactions, but you must be absolutely sure they are real and correct."

            else:
                raise Exception(f"Invalid raw reactions generation mode: {mode}")

            instruct = \
            f"Please provide a comprehensive and diverse list of documented chemical reactions involving {chem_name}, where it appears as {compound_role}. " \
            f"{special_requirement} "
            "Write the reactions as schemes using only '->' and '+' symbols. Use the full chemical names of substances instead of formulas or generic terms. " \
            "Do not include balancing coefficients, comments, or any markup - only the reaction schemes themselves one per line. " \
            "If no such substance exists or no documented reactions are available, return 'None'."

            instruct_revalidate = \
            "Please, review the provided reactions. Identify any erroneous reactions and correct them where possible. Return revised reactions list that comply with the initial requirements."


            messages = [
                {"role": "system", "content": ""},
                {"role": "user", "content": instruct}
            ]

            models_schedule = [self.gpt_oss, self.qwen]
            
            for curr_model in models_schedule:
                model = curr_model
                response = self.__fetch_llm_answer(messages, model)
                reactions = self.__get_reactions_from_response(response)
                if reactions:
                    break
            else:
                self.log(f"Failed to fetch reactions for '{chem_name}'")
                return None
            
            #self.log(f"Got {len(reactions)} initial reactions for {chem_name}")
            
            response = '\n'.join(reactions)

            messages.append({"role": "assistant", "content": response})
            messages.append({"role": "user", "content": instruct_revalidate})

            response = self.__fetch_llm_answer(messages, model)
            reactions_revalid = self.__get_reactions_from_response(response)
            if not reactions_revalid:
                self.log(f"Failed to fetch revalidated reactions for '{chem_name}'. Assuming all are valid...")
                reactions_revalid = reactions
            
            self.log(f"Got {len(reactions_revalid)} reactions for {chem_name} with '{model}'; CTT: {self.completion_tokens_total}")
        
        except Exception as e:
            self.log(f"Exception in '__fetch_raw_reactions': {e}")
            return None
        
        return {'cid': chem['cid'], 'reactions': reactions_revalid}



    def get_raw_reactions(self, max_workers=1):        
        processed = self.__get_processed_entries(self.raw_reactions_fn, 'cid')
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.raw_reactions_fn, 'a') as f_out:
            futures = []
            for chem in self.chems:
                if chem['cid'] not in processed:
                    futures.append(executor.submit(self.__fetch_raw_reactions, chem))
            
            for future in as_completed(futures):
                res = future.result()
                if res is None:
                    continue
                f_out.write(json.dumps(res) + '\n')
                f_out.flush()
    

    def get_uncommon_raw_reactions_for_wiki_chems(self, max_workers=1):        
        with open(self.wiki_chems_fn) as f:
            wiki_chems_cids = set([json.loads(x)['cid'] for x in f.read().strip().split('\n')])
        
        processed = self.__get_processed_entries(self.wiki_raw_reactions_fn, 'cid')

        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.wiki_raw_reactions_fn, 'a') as f_out:
            futures = []
            for chem in self.chems:
                cid = chem['cid']
                if cid not in processed and cid in wiki_chems_cids:
                    futures.append(executor.submit(self.__fetch_raw_reactions, chem, "documented_less_common_rp"))
            
            for future in as_completed(futures):
                res = future.result()
                if res is None:
                    continue
                f_out.write(json.dumps(res) + '\n')
                f_out.flush()
    

    def get_rare_raw_reactions_for_top_chems(self, max_workers=1):        
        hazards_chems = dict()
        with open(self.hazards_chems_fn) as f:
            for line in f:
                hazard = json.loads(line)
                hazards_chems[hazard['cid']] = hazard
        
        processed = self.__get_processed_entries(self.top_rare_raw_reactions_fn, 'cid')

        def is_top_chem(hazards):
            for pic in hazards['pictograms']:
                if pic in {'GHS01', 'GHS03', 'GHS06'}:
                    return True
            return False

        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.top_rare_raw_reactions_fn, 'a') as f_out:
            futures = []
            for chem in self.chems:
                cid = chem['cid']
                if cid not in processed and cid in hazards_chems and is_top_chem(hazards_chems[cid]):
                    futures.append(executor.submit(self.__fetch_raw_reactions, chem, "rare_rp"))
            
            for future in as_completed(futures):
                res = future.result()
                if res is None:
                    continue
                f_out.write(json.dumps(res) + '\n')
                f_out.flush()
    

    def get_uncommon_raw_reactions_for_wiki_chems_products_only(self, max_workers=1):        
        processed = self.__get_processed_entries(self.products_wiki_raw_reactions_fn, 'cid')

        staged_chems = [chem for chem in self.chems if 'wiki' in chem and chem['cid'] not in processed]
        random.shuffle(staged_chems)
        print(f"Staged {len(staged_chems)} compounds")

        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.products_wiki_raw_reactions_fn, 'a') as f_out:
            futures = []
            for chem in staged_chems:
                futures.append(executor.submit(self.__fetch_raw_reactions, chem, "documented_less_common_p"))
            
            for future in as_completed(futures):
                res = future.result()
                if res is None:
                    continue
                f_out.write(json.dumps(res) + '\n')
                f_out.flush()

    

    def __get_verdicts_bool_from_response(self, response: str):
        verdicts = []
        for line in response.split('\n'):
            if 'invalid' in line.lower():
                verdicts.append(False)
            elif 'valid' in line.lower():
                verdicts.append(True)
        
        return verdicts
    

    def __validate_raw_reactions(self, reactions, valid_cnt):
        try:
            instruct_validate = \
            "You will be given a list of unbalanced chemical reaction schemes. " \
            "For each scheme, determine if the reaction is chemically possible and whether the listed reactants and products are correct. " \
            "Assume that necessary reaction conditions (including harsh), solvents, or catalysts are implied even if not explicitly shown. " \
            "All reactions listed are theoretical and intended solely for academic purposes, not for practical experimentation. " \
            "If the reaction is valid, output only 'Valid'. If it is not valid, output only 'Invalid'. " \
            "Print one result per line, preceded by the reaction index. " \
            "Do not include any additional text."

            models_schedule = [self.gpt_oss, self.qwen]

            for try_i, model in enumerate(models_schedule):

                def extract_verdicts_from_response(response):
                    return ['invalid' not in verd.lower() and 'valid' in verd.lower() for verd in response.split('\n')]

                valid_i = 0
                mistakes_cnt = 0
                mistakes_thr = 3
                confidences = [0.0 for _ in range(len(reactions))]
                confidence_thr = 0.5
                bad = False
                results = []
                while valid_i < valid_cnt and len(reactions) > 0:
                    reactions_str = '\n'.join([f"{i+1}. {react['reaction']}" for i, react in enumerate(reactions)])
                    messages = [
                        {"role": "system", "content": ""},
                        {"role": "user", "content": f"{instruct_validate}\n{reactions_str}"}
                    ]
                    response = self.__fetch_llm_answer(messages, model, reasoning_effort='low')
                    verdicts = extract_verdicts_from_response(response)
                    if len(verdicts) != len(reactions):
                        if mistakes_cnt == mistakes_thr:
                            self.log(f"Falling to another model due to mistakes ({try_i+1}/{len(models_schedule)}) ('{model}')")
                            bad = True
                            break
                        mistakes_cnt += 1
                        continue

                    finished_indices = set()
                    remaining_tries = (valid_cnt-valid_i-1)
                    for i in range(len(verdicts)):
                        confidences[i] += verdicts[i] / valid_cnt
                        max_confidence = remaining_tries / valid_cnt + confidences[i]
                        est_confidence = (confidences[i] + max_confidence) / 2
                        if confidences[i] >= confidence_thr or max_confidence <= confidence_thr:
                            react = reactions[i].copy()
                            reaction_str = react['reaction']
                            react['valid'] = confidences[i] > confidence_thr
                            react['confidence'] = confidences[i]
                            react['source'] = model
                            results.append(react)
                            finished_indices.add(i)
                            self.log(f"Processed reaction '{reaction_str}'; confidence: {est_confidence:.2f}; CTT: {self.completion_tokens_total}")

                    reactions = [reactions[i] for i in range(len(reactions)) if i not in finished_indices]
                    confidences = [confidences[i] for i in range(len(confidences)) if i not in finished_indices]
                    
                    valid_i += 1

                if bad:
                    continue
                
                return results
        
        except Exception as e:
            self.log(f"Exception in '__validate_raw_reaction': {e}")
            return None
    

    def validate_raw_reactions(self, raw_reactions_fn=None, max_workers=1):
        if raw_reactions_fn is None:
            raw_reactions_fn = self.raw_reactions_fn

        entries = self._load_jsonl(raw_reactions_fn)
        
        processed = self.__get_processed_entries(self.raw_reactions_verdict_fn, 'reaction')
        
        reactions = []

        for entry in entries:
            cid = entry['cid']
            reactions_curr = entry['reactions']
            for react in reactions_curr:
                if react not in processed:
                    if self._parse_reaction_str(react['reaction'])[0] is not None:
                        reactions.append({'cid': cid, 'reaction': react})
        
        reactions_batch_size = 10
        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.raw_reactions_verdict_fn, 'a') as f_out:
            futures = []
            i = 0
            while i < len(reactions):
                futures.append(executor.submit(self.__validate_raw_reactions, reactions[i:i+reactions_batch_size], 9))
                i += reactions_batch_size
            
            for future in as_completed(futures):
                res = future.result()
                if res:
                    for react in res:
                        f_out.write(json.dumps(react) + '\n')
                    f_out.flush()

    

    def __get_chem_description(self, chem, valid_cnt):
        chem_name = chem['cmpdname']

        try:
            instruct = \
            "Write an engaging and informative plain text description of the compound these synonyms refer to. " \
            "You may decide freely what aspects to include — composition, properties, uses, history, relevance, or anything else meaningful — " \
            "as long as the result is interesting to read and based on reliable chemical knowledge. " \
            "If nothing meaningful can be said about it, respond exactly with the word None.\n" \
            "Result should be in plain text format. Markdown or other formatting types are strictly forbidden." \
            "Guidelines:\n" \
            "- Everything you write must be true or based on well-accepted chemical knowledge.\n" \
            "- You may include lesser-known but reliable information, but do not guess or invent facts.\n" \
            "- The text should feel lively, colorful, and pleasant to read — avoid dull or purely technical writing.\n" \
            "- You may write at any length, as long as the text remains engaging and free of boring or redundant details.\n" \
            "- The output should be a smooth, uninterrupted narrative without explicit topic breaks or labeled sections.\n"

            instruct_validate = \
            "Decide if the following text is factually correct and does not contain misleading or speculative information. " \
            "If the description is fully correct or at least safely accurate within general chemical knowledge, return 'Valid'. " \
            "If it contains anything blatantly wrong, incorrect, or unverified, return 'Invalid'. " \
            "Do not write anything except 'Valid' or 'Invalid'."

            instruct_fix = \
            "Review the provided chemical compound description for errors or inconsistencies. " \
            "Identify any mistakes and return a corrected, accurate version of the description, preserving its original narrative tone. " \
            "Return only the corrected description, without any explanations or extra text."

            av_confidence = 0
            max_syns = 3
            synonyms = ', '.join(list(map(lambda x: f'"{x}"', chem['cmpdsynonym'][:max_syns])))
            models_schedule = [self.gpt_oss, self.gpt_oss, self.qwen, self.qwen]
            for try_i, model in enumerate(models_schedule):
                messages = [
                    {"role": "system", "content": ""},
                    {"role": "user", "content": f"{synonyms}\n\n{instruct}"}
                ]

                description = self.__fetch_llm_answer(messages, model)
                
                if len(description) < 20:
                    self.log(f"Failed to generate description for {chem_name} on try {try_i+1}/{len(models_schedule)} ('{model}')")
                    continue
                

                def validate_description(descr, confidence_thr):
                    nonlocal instruct_validate, model

                    messages = [
                        {"role": "system", "content": ""},
                        {"role": "user", "content": f"{instruct_validate}\n\n{descr}"}
                    ]
                    confidence = 0
                    for valid_i in range(valid_cnt):
                        verdict = self.__fetch_llm_answer(messages, model)
                        if 'invalid' not in verdict.lower():
                            if 'valid' in verdict.lower():
                                confidence += 1 / valid_cnt

                        remaining_tries = (valid_cnt-valid_i-1)
                        max_confidence = remaining_tries / valid_cnt + confidence
                        if max_confidence < confidence_thr or confidence >= confidence_thr:
                            confidence += remaining_tries / 2 / valid_cnt
                            break

                    return confidence
                
                confidence_thr = 0.49
                confidence = validate_description(description, confidence_thr)
                
                if confidence >= confidence_thr:
                    self.log(f"('{model}') Generated description for {chem_name} of length {len(description)}. confidence: {confidence}; CTT: {self.completion_tokens_total}")
                    return {'cid': chem['cid'], 'description': description, 'confidence': confidence, 'source': model}

                messages = [
                    {"role": "system", "content": ""},
                    {"role": "user", "content": f"{instruct_fix}\n\n{description}"}
                ]
                description = self.__fetch_llm_answer(messages, model)
                confidence = validate_description(description, confidence_thr)
                if confidence >= confidence_thr:
                    self.log(f"('{model}') Generated description for {chem_name} of length {len(description)} after fixing. confidence: {confidence}; CTT: {self.completion_tokens_total}")
                    return {'cid': chem['cid'], 'description': description, 'confidence': confidence, 'source': model}

                self.log(f"Low validation confidence for {chem_name} on try {try_i+1}/{len(models_schedule)} ('{model}'): {confidence}")
                av_confidence += confidence / len(models_schedule)
            
            self.log(f"Failed to generate description for {chem_name} due to low validation confidence: {av_confidence}")
                
        except Exception as e:
            self.log(f"Exception during description generation for '{chem_name}': {e}")

        return None
        
    

    def get_chems_descriptions(self, max_workers=1):
        chems_power = dict()
        for chem in self.chems:
            chems_power[chem['cid']] = 0

        with open(self.reactions_parsed_fn) as f:
            for line in f:
                reaction = json.loads(line)
                all_cids = [x['cid'] for x in reaction['reagents']] + [x['cid'] for x in reaction['products']]
                for cid in all_cids:
                    chems_power[cid] += 1
        
        self.chems.sort(key=lambda x: chems_power[x['cid']], reverse=True)

        processed = self.__get_processed_entries(self.chems_descriptions_fn, 'cid')

        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.chems_descriptions_fn, 'a') as f_out:
            futures = []
            for chem in self.chems:
                if chem['cid'] not in processed and 'wiki' in chem:
                    futures.append(executor.submit(self.__get_chem_description, chem, 6))
            
            for future in as_completed(futures):
                res = future.result()
                if res:
                    f_out.write(json.dumps(res) + '\n')
                    f_out.flush()
    

    def __get_reactions_description(self, reactions, valid_cnt):
        try:
            instruct = \
            "Please generate detailed short plain text descriptions for each of the chemical reactions provided below. " \
            "For each reaction, include the general type of reaction, its purpose if applicable, key transformations, typical reaction conditions, " \
            "and common solvents or catalysts if needed. Prefix each description with the corresponding reaction number. " \
            "Provide only the descriptions, do not add any extra text or commentary and do not restate or repeat the reactions themselves. Provide only information that you are confident is correct. " \
            "Do not use formatting, but plain text only. If you cannot provide reliable information for a specific reaction, simply write 'None' as the description."
            
            instruct_validate = \
            "You will be given pairs of chemical reaction schemes and their descriptions. Each reaction scheme shows reagents and products (unbalanced). " \
            "Your task is to validate whether each description accurately describes the corresponding reaction. " \
            "Return 'Valid' if the description is chemically accurate and contains no significant errors. " \
            "Return 'Invalid' if the description contains factual errors, incorrect mechanisms, impossible conditions, or misidentified reaction types. " \
            "Minor imprecision is acceptable if the core chemistry is correct. " \
            "Respond with only 'Valid' or 'Invalid' for each pair, one per line, in the same order as presented."

            instruct_fix = \
            "You will be given pairs of chemical reaction schemes and their descriptions. " \
            "Find and correct the chemical errors if any in each description while keeping all other text as similar as possible to the original. " \
            "Provide only the corrected descriptions, one per line, in the same order as presented. " \
            "If description is wrong and you cannot reliably correct a description, write 'None'."

            models_schedule = [self.gpt_oss, self.qwen]
            for try_i, model in enumerate(models_schedule):
                reactions_formatted_str = '\n'.join([f"{i}. {react['reaction']}" for i, react in enumerate(reactions)])
                messages = [
                    {"role": "system", "content": ""},
                    {"role": "user", "content": f"{instruct}\n\n{reactions_formatted_str}"}
                ]

                def extract_descriptions_from_response(response):
                    descriptions = re.split(r'\n(?=\d+[.)]\s)', response)
                    descriptions = [d.strip() for d in descriptions if d.strip()]
                    descriptions = [re.sub(r'^\d+[.)]\s+', '', d) for d in descriptions]
                    return descriptions
                
                response = self.__fetch_llm_answer(messages, model)
                descriptions = extract_descriptions_from_response(response)
                if len(descriptions) != len(reactions):
                    self.log(f"Failed to get descriptions ({len(descriptions)} != {len(reactions)}) ('{model}') ({try_i+1}/{len(models_schedule)})")
                    continue

                def filter_descriptions_reactions(descriptions, reactions):
                    indices = [i for i, desc in enumerate(descriptions) if len(desc) > 20]
                    return [descriptions[i] for i in indices], [reactions[i] for i in indices]

                descriptions, reactions = filter_descriptions_reactions(descriptions, reactions)
                
                def extract_verdicts_from_response(response):
                    return ['invalid' not in verd.lower() and 'valid' in verd.lower() for verd in response.split('\n')]

                valid_i = 0
                mistakes_cnt = 0
                mistakes_thr = 2
                confidences = [0 for _ in range(len(descriptions))]
                confidence_thr = 0.39
                results = []
                while valid_i < valid_cnt and len(descriptions) > 0:
                    formatted_descriptions_str = []
                    for i in range(len(descriptions)):
                        react = reactions[i]['reaction']
                        formatted_descriptions_str.append(f'{i+1}. Scheme: "{react}"\nDescription: "{descriptions[i]}"')
                    formatted_descriptions_str = '\n\n'.join(formatted_descriptions_str)
                    messages = [
                        {"role": "system", "content": ""},
                        {"role": "user", "content": f"{instruct_validate}\n\n{formatted_descriptions_str}"}
                    ]
                    response = self.__fetch_llm_answer(messages, model)
                    verdicts = extract_verdicts_from_response(response)
                    if len(verdicts) != len(descriptions):
                        if mistakes_cnt == mistakes_thr:
                            self.log("Returning early due to mistakes at validation phase")
                            return results
                        mistakes_cnt += 1
                        continue
                    
                    finished_indices = set()
                    remaining_tries = (valid_cnt-valid_i-1)
                    for i in range(len(verdicts)):
                        confidences[i] += verdicts[i] / valid_cnt
                        max_confidence = remaining_tries / valid_cnt + confidences[i]
                        est_confidence = (confidences[i] + max_confidence) / 2
                        if confidences[i] >= confidence_thr:
                            reaction_str = reactions[i]['reaction']
                            self.log(f"Generated description for '{reaction_str}'; confidence: {est_confidence}")
                            results.append({'rid': reactions[i]['rid'], 'description': descriptions[i], 'confidence': est_confidence, 'source': model})
                            finished_indices.add(i)
                        elif max_confidence < confidence_thr:
                            finished_indices.add(i)
                    reactions = [reactions[i] for i in range(len(reactions)) if i not in finished_indices]
                    descriptions = [descriptions[i] for i in range(len(descriptions)) if i not in finished_indices]
                    confidences = [confidences[i] for i in range(len(confidences)) if i not in finished_indices]

                    valid_i += 1
                
                return results

        
        except Exception as e:
            self.log(f"Exception during reactions description generation: {e}")
        
        return None


    def get_reactions_descriptions(self, max_workers=1):
        chems_power = dict()
        for chem in self.chems:
            chems_power[chem['cid']] = 0

        reactions = []
        with open(self.reactions_parsed_fn) as f:
            for line in f:
                reaction = json.loads(line)
                all_cids = [x['cid'] for x in reaction['reagents']] + [x['cid'] for x in reaction['products']]
                for cid in all_cids:
                    chems_power[cid] += 1
                reactions.append(reaction)

        reactions_power = dict()
        for react in reactions:
            all_cids = [x['cid'] for x in reaction['reagents']] + [x['cid'] for x in reaction['products']]
            reactions_power[react['rid']] = sum(chems_power[x] for x in all_cids) / len(all_cids)
        
        processed = self.__get_processed_entries(self.reactions_descriptions_fn, 'rid')
        reactions = list(filter(lambda x: x['rid'] not in processed and x['source'] != 'ord', reactions))
        reactions.sort(key=lambda x: reactions_power[x['rid']], reverse=True)

        reactions_batch_size = 6
        with ThreadPoolExecutor(max_workers=max_workers) as executor, open(self.reactions_descriptions_fn, 'a') as f_out:
            futures = []
            for i in range(0, len(reactions), reactions_batch_size):
                reactions_arg = [{'reaction': self._get_reaction_as_str(x), 'rid': x['rid']} for x in reactions[i:i + reactions_batch_size]]
                futures.append(executor.submit(self.__get_reactions_description, reactions_arg, 6))
            
            for future in as_completed(futures):
                res = future.result()
                if res:
                    self.log(f"\nLost {reactions_batch_size-len(res)}/{reactions_batch_size}; CTT: {self.completion_tokens_total}\n")
                    for entry in res:
                        f_out.write(json.dumps(entry) + '\n')
                    f_out.flush()




if __name__ == "__main__":
    chemsllm = ChemsLLMFetch("data/", "chemistry")
    #chemsllm.get_raw_reactions(max_workers=20)
    #chemsllm.process_raw_reactions('raw_reactions.jsonl')
    #chemsllm.deduplicate_raw_reactions()
    #chemsllm.validate_raw_reactions(max_workers=20)
    #chemsllm.fix_broken_raw_reactions(max_workers=1)
    #print(chemsllm.find_all_unicode_chars_in_raw_reactions())
    #chemsllm.parse_raw_llm_reactions()
    #chemsllm.fetch_names_from_pubchem()
    #chemsllm.organize_chems_file()
    #chemsllm.balance_parsed_reactions()
    #chemsllm.find_unbalancing_chems()
    #chemsllm.get_uncommon_raw_reactions_for_wiki_chems(max_workers=20)
    #chemsllm.validate_raw_reactions(raw_reactions_fn="data/top_rare_raw_reactions.jsonl", max_workers=20)
    #chemsllm.generate_chems_structures_svg()
    #chemsllm.get_rare_raw_reactions_for_top_chems(max_workers=20)
    #chemsllm.generate_organic_marks_for_chems()
    #chemsllm.generate_edges()
    #chemsllm.filter_chems_with_invalid_smiles()
    #chemsllm.generate_chems_fingerprints()
    #chemsllm.merge_wiki_chems()
    #chemsllm.make_chems_smiles_canonic()
    #chemsllm.generate_chems_bertz_complexity()
    #chemsllm.show_bertz_complexity_diff()
    #chemsllm.parse_raw_ord_reactions('cleaned_ord.jsonl')
    #chemsllm.fetch_smiles_from_pubchem()
    #chemsllm.fetch_chems_cids_from_pubchem('cids.txt')
    #chemsllm.merge_parsed_files("data/reactions_parsed/merged_reactions_parsed.jsonl", "data/reactions_parsed/reactions_parsed_ord.jsonl", "data/reactions_parsed/reactions_parsed.jsonl")
    #chemsllm.balance_parsed_reactions()
    #chemsllm.generate_edges()
    #chemsllm.populate_db()
    #chemsllm.deduplicate_chems_rebind_reactions()
    #chemsllm.fix_details()
    #chemsllm.fix_reactions()
    #chemsllm.test()
    #chemsllm.get_uncommon_raw_reactions_for_wiki_chems_products_only(max_workers=20)
    #chemsllm.get_chems_descriptions(max_workers=20)
    #chemsllm.get_reactions_descriptions(max_workers=20)
    #chemsllm.validate_raw_reactions(raw_reactions_fn="data/wiki_products_raw_reactions.jsonl", max_workers=20)
    #chemsllm.clean_data_populate_tables(rehash_required=True)
    #chemsllm.extract_chems_cas_numbers()
    #chemsllm.get_background_substances(20)
    #chemsllm.get_commonnes_chems_sorting()
    #chemsllm.filter_ord_reactions()
    #chemsllm.extract_radicals_list('radicals.jsonl')
    #chemsllm.clean_ord_reactions_from_radicals()
    #chemsllm.filter_ord_reactions()
    #chemsllm.fetch_elements()
    #chemsllm.compute_reactions_thermo_xtb('reaction_enthalpies.jsonl')
    #chemsllm.test()
    #chemsllm.merge_details()
    #chemsllm.merge_reactions()
    #chemsllm.get_conflicting_synonyms('conflict.txt')
    #chemsllm.map_crc_chems_to_cids()
    #chemsllm.fetch_names_from_pubchem('data/crc_unmapped_names.txt')
    #chemsllm.get_commonnes_chems_sorting()
    #chemsllm.extract_pubchem_dump_to_chems('pubchem.json', override=True)
    #chemsllm.get_chems_cids('cids.txt')
    #chemsllm.resolve_conflicting_synonyms()
    #chemsllm.test()

    

