import json
import pdfplumber as pdfp
import re
import os

from tempfile import NamedTemporaryFile

from chems_pubchem_parse import ChemsParsePubchem


class ChemsCRC(ChemsParsePubchem):

    def __init__(self, data_dir):
        super().__init__(data_dir)

        self.crc_assests_dir = os.path.join(self.data_dir, 'assets', 'crc_handbook')

        self.crc_flammability_pdf_fn = os.path.join(self.crc_assests_dir, 'flammability.pdf')
        self.crc_inorganic_constants_pdf_fn = os.path.join(self.crc_assests_dir, 'inorganic_constants.pdf')
        self.crc_organic_constants_pdf_fn = os.path.join(self.crc_assests_dir, 'organic_constants.pdf')

        self.crc_inorganic_abbreviations_map_fn = os.path.join(self.crc_assests_dir, 'inorganic_abbreviations_map.txt')
        self.crc_organic_abbreviations_map_fn = os.path.join(self.crc_assests_dir, 'organic_abbreviations_map.txt')

        self.crc_inorganic_constants_fn = os.path.join(self.data_dir, 'crc_handbook', 'inorganic_constants.jsonl')
        self.crc_organic_constants_fn = os.path.join(self.data_dir, 'crc_handbook', 'organic_constants.jsonl')
        self.crc_flammability_fn = os.path.join(self.data_dir, 'crc_handbook', 'flammability.jsonl')

        self.crc_unmapped_names_fn = os.path.join(self.data_dir, 'crc_handbook', 'crc_unmapped_names.txt')

        self._file_sorting_prefs[self.crc_flammability_fn] = 'name'
        self._file_sorting_prefs[self.crc_inorganic_constants_fn] = 'name'
        self._file_sorting_prefs[self.crc_organic_constants_fn] = 'name'

        self.solubility_cats = {'insoluble', 'miscible', 'soluble', 'slightly soluble', 'very soluble'}


    def __parse_page_raw(self, words, chars, column_headers_fields_map: dict, headers_cutoff_thr_bottom: float, headers_min_top: float, mandatory_field: str, row_min_height):
        column_headers = list(column_headers_fields_map.keys())
        header_to_i = {h: i for i, h in enumerate(column_headers)}
        column_placeholder = 2**30
        column_left_x = [column_placeholder for _ in column_headers]
        for w in words:
            if w['top'] < headers_min_top:
                continue
            if w['bottom'] > headers_cutoff_thr_bottom:
                continue
            text = w['text']

            if text in header_to_i:
                curr_value = column_left_x[header_to_i[text]]
                if curr_value > w['x0']:
                    column_left_x[header_to_i[text]] = w['x0']
        
        if any(x == column_placeholder for x in column_left_x):
            print(column_left_x)
            raise Exception("Failed to find all defined column headers")

        tops = sorted(w['top'] for w in words if w['bottom'] > headers_cutoff_thr_bottom)
        rows_top_y = [tops[0]]
        for top in tops[1:]:
            if top - rows_top_y[-1] > row_min_height:
                rows_top_y.append(top)
        
        table = [[[] for c in column_left_x] for r in rows_top_y]

        def binary_search_index(arr, value):
            left, right = 0, len(arr) - 1
            result = -1
            
            while left <= right:
                mid = left + (right - left) // 2
                
                if arr[mid] < value:
                    result = mid
                    left = mid + 1
                else:
                    right = mid - 1
            
            return result

        char_height_thr = 5.0
        for c in chars:
            if c['height'] < char_height_thr:
                continue
            if c['bottom'] < headers_cutoff_thr_bottom:
                continue
            text = c['text']
            row_i = binary_search_index(rows_top_y, c['bottom'])
            col_i = binary_search_index(column_left_x, c['x1'])
            table[row_i][col_i].append(c)
        
        for row_i in range(len(table)):
            for col_i in range(len(table[row_i])):
                table[row_i][col_i].sort(key=lambda c: (round(c['top'], 1), c['x0']))
                table[row_i][col_i] = ''.join([c['text'] for c in table[row_i][col_i]])
        

        def process_entry(entry):
            entry = entry.strip()
            return entry if entry else None
        
        results = []
        for row_i in range(len(table)):
            curr_res = dict()
            for h_i, header in enumerate(column_headers):
                res_field = column_headers_fields_map[header]
                if not res_field:
                    continue
                curr_res[res_field] = process_entry(table[row_i][h_i])
            results.append(curr_res)
        

        results_merged = []
        for entry in results:
            if entry[mandatory_field] is not None:
                results_merged.append(entry)
                continue

            if all(entry[field] is None for field in entry.keys()):
                continue

            for field in entry.keys():
                if entry[field] is not None and results_merged[-1][field] is not None:
                    # Filtering invisible garbage text
                    if '.indb' in entry[field]:
                        continue
                    results_merged[-1][field] += f" {entry[field]}"
        

        return results_merged




    def _parse_organic_constants_raw(self, out_fn, start_page=1, last_page=-1):
        with pdfp.open(self.crc_organic_constants_pdf_fn) as pdf:
            page_i = start_page-1
            last_page_i = last_page if last_page != -1 else len(pdf.pages)
            column_headers_fields_map = {"No.": 'ind',
                                        "Name": 'name',
                                        "Synonym": 'synonym',
                                        "Mol.": None,
                                        "CAS": 'cas',
                                        "Wt.": None,
                                        "Physical": 'physical_form',
                                        "mp/\u02daC": 'mp',
                                        "bp/\u02daC": 'bp',
                                        "den/": 'density',
                                        "n": 'refractive_index',
                                        "Solubility": 'solubility'
                                        }
            headers_cutoff_thr_bottom = 80.0
            headers_min_top = 45.0
            row_min_height = 5.0
            results = []
            for page_i in range(start_page-1, last_page_i, 2):
                page = pdf.pages[page_i]
                words = page.extract_words()
                chars = page.chars

                try:
                    page_result = self.__parse_page_raw(words, chars, column_headers_fields_map, headers_cutoff_thr_bottom, headers_min_top, 'ind', row_min_height)
                except Exception as e:
                    print(f"Exception: '{e}'")
                    return

                results += page_result
            
            self._write_jsonl(results, out_fn)



    def _parse_inorganic_constants_raw(self, out_fn, start_page=1, last_page=-1):
        with pdfp.open(self.crc_inorganic_constants_pdf_fn) as pdf:
            start_page_i = start_page-1
            page_i = start_page_i
            last_page_i = last_page if last_page != -1 else len(pdf.pages)
            row_min_height = 5.0
            column_headers_fields_map = {"No.": 'ind',
                                        "Name": 'name',
                                        "Formula": None,
                                        "CAS": 'cas',
                                        "Mol.": None,
                                        "Physical": 'physical_form',
                                        "mp/°C": 'mp',
                                        "bp/°C": 'bp',
                                        "Density": 'density',
                                        "g/100": 'solubility_aq',
                                        "Qualitative": 'solubility'
                                        }
            
            results = []
            for page_i in range(start_page-1, last_page_i):
                page = pdf.pages[page_i]
                words = page.extract_words()
                chars = page.chars
                is_first_page = page_i==start_page_i
                headers_cutoff_thr_bottom = 550.0 if is_first_page else 80.0
                headers_min_top = 525 if is_first_page else 60.0

                try:
                    page_result = self.__parse_page_raw(words, chars, column_headers_fields_map, headers_cutoff_thr_bottom, headers_min_top, 'ind', row_min_height)
                except Exception as e:
                    print(f"Exception: '{e}'")
                    return

                results += page_result

            self._write_jsonl(results, out_fn)


    def _parse_flammability_raw(self, out_fn, start_page=1, last_page=-1):
        with pdfp.open(self.crc_flammability_pdf_fn) as pdf:
            start_page_i = start_page-1
            page_i = start_page_i
            last_page_i = last_page if last_page != -1 else len(pdf.pages)
            row_min_height = 5.0
            column_headers_fields_map = {"mol.": 'formula',
                                        "name": 'name',
                                        "t": None,
                                        "fp/\u00b0C": 'flash_point',
                                        "fl.": "flash_limits",
                                        "it/\u00b0C": 'ignition_temp'
                                        }
            
            results = []
            for page_i in range(start_page-1, last_page_i):
                page = pdf.pages[page_i]
                words = page.extract_words()
                chars = page.chars
                is_first_page = page_i==start_page_i
                headers_cutoff_thr_bottom = 430.0 if is_first_page else 96.0
                headers_min_top = 409 if is_first_page else 85.0

                page_result = self.__parse_page_raw(words, chars, column_headers_fields_map, headers_cutoff_thr_bottom, headers_min_top, 'formula', row_min_height)
                try:
                    pass
                except Exception as e:
                    print(f"Exception: '{e}'")
                    return

                page_result = list(filter(lambda x: x['name'] and x['formula'], page_result))
                for res in page_result:
                    res.pop('formula')

                results += page_result

            self._write_jsonl(results, out_fn)
    

    def __clean_solubility(self, sol_str, ab_map):
        sol_str = re.sub(r',;', '', sol_str)
        sol_str = re.sub(r'\s+', ' ', sol_str)
        match = re.findall(r'(\bi-|\b[a-zA-Z0-9]+\b)', sol_str, re.ASCII)
        cleaned = dict()
        curr_sol_cat = None
        for word in match:
            if word in ab_map:
                word = word.replace(word, ab_map[word])
                if word in self.solubility_cats:
                    curr_sol_cat = word
                else:
                    if curr_sol_cat is None:
                        return None
                    cleaned[word] = curr_sol_cat
        
        return cleaned


    def __clean_physical_form(self, phys_str, ab_map):
        phys_str = re.sub(r'\s+', ' ', phys_str)
        match = re.findall(r'(\bi-|\b[a-zA-Z0-9]+\b|\W)', phys_str, re.ASCII)
        cleaned = ''
        for word in match:
            if word in ab_map:
                word = word.replace(word, ab_map[word])
            cleaned += word
        
        return cleaned


    def __check_approx(self, string):
        return True if re.search(r'[\u2248<>]', string) else False

    def __clean_mp(self, mp_str):
        dec = 'dec' in mp_str
        approx = self.__check_approx(mp_str)
        value = self.__clean_float_value(mp_str)
        cleaned = {'decomposes': dec, 'value': value, 'approx': approx}

        return cleaned if dec or (value is not None) else None


    def __clean_bp(self, bp_str):
        sub = 'sub' in bp_str or 'sp' in bp_str
        dec = 'dec' in bp_str
        approx =self. __check_approx(bp_str)
        value = self.__clean_float_value(bp_str)
        cleaned = {'sublimes': sub, 'decomposes': dec, 'value': value, 'approx': approx}

        return cleaned if dec or sub or (value is not None) else None


    def __clean_float_value(self, val_str):
        value_match = re.findall(r'-?\d+(?:\.\d+)?', val_str)
        if len(value_match) != 1:
            value = None
        else:
            value = float(value_match[0])
        
        return value


    def __clean_unicode(self, string):
        if not string:
            return None

        unicode_map = {
            '\ufb02': 'fl',
            '\ufb01': 'fi',
            '\u2019': "'",
            '\u2013': '-'
        }

        return ''.join(map(lambda c: unicode_map[c] if c in unicode_map else c, string))


    def __clean_flash_point(self, fp_str):
        fp_str = self.__clean_unicode(fp_str)
        approx = self.__check_approx(fp_str)
        value = self.__clean_float_value(fp_str)
        cleaned = {'approx': approx, 'value': value}

        return cleaned if value is not None else None


    def __clean_flash_limits(self, fl_str):
        fl_str = self.__clean_unicode(fl_str)
        fl_str = re.sub(r'\s+', '', fl_str)

        return fl_str


    def __clean_ignition_temp(self, it_str):
        approx = self.__check_approx(it_str)
        value = self.__clean_float_value(it_str)
        cleaned = {'approx': approx, 'value': value}

        return cleaned if value is not None else None


    def _clean_organic_constants(self, input_fn, output_fn=None):
        if output_fn is None:
            output_fn = input_fn

        entries = self._load_jsonl(input_fn)

        with open(self.crc_organic_abbreviations_map_fn) as f:
            ab_map = [x.split('|') for x in f.read().strip().split('\n')]
            ab_map = dict(ab_map)

        for entry in entries:
            if not entry['name']:
                continue

            if entry['physical_form']:
                entry['physical_form'] = self.__clean_physical_form(entry['physical_form'], ab_map)
            
            if entry['solubility']:
                entry['solubility'] = self.__clean_solubility(entry['solubility'], ab_map)
            
            if entry['mp']:
                entry['mp'] = self.__clean_mp(entry['mp'])

            if entry['bp']:
                entry['bp'] = self.__clean_bp(entry['bp'])
            
            if entry['density']:
                entry['density'] = self.__clean_float_value(entry['density'])
            
            if entry['refractive_index']:
                entry['refractive_index'] = self.__clean_float_value(entry['refractive_index'])

            entry['name'] = self.__clean_unicode(entry['name'])
            entry['synonym'] = self.__clean_unicode(entry['synonym'])

            entry.pop('ind')


        self._write_jsonl(entries, output_fn)


    def _clean_inorganic_constants(self, input_fn, output_fn=None):
        if output_fn is None:
            output_fn = input_fn

        entries = self._load_jsonl(input_fn)

        with open(self.crc_inorganic_abbreviations_map_fn) as f:
            ab_map = [x.split('|') for x in f.read().strip().split('\n')]
            ab_map = dict(ab_map)

        for entry in entries:
            if not entry['name']:
                continue

            if entry['physical_form']:
                entry['physical_form'] = self.__clean_physical_form(entry['physical_form'], ab_map)
            
            if entry['solubility']:
                entry['solubility'] = self.__clean_solubility(entry['solubility'], ab_map)
            
            if entry['mp']:
                entry['mp'] = self.__clean_mp(entry['mp'])

            if entry['bp']:
                entry['bp'] = self.__clean_bp(entry['bp'])
            
            if entry['density']:
                entry['density'] = self.__clean_float_value(entry['density'])
            
            if entry['solubility_aq']:
                entry['solubility_aq'] = self.__clean_float_value(entry['solubility_aq'])

            entry['name'] = self.__clean_unicode(entry['name'])

            entry.pop('ind')


        self._write_jsonl(entries, output_fn)



    def _map_inorganic_organic_constants(self, input_fn, output_fn=None):
        if output_fn is None:
            output_fn = input_fn

        entries = self._load_jsonl(input_fn)

        parsed_entries = []
        processed_cids = set()
        for entry in entries:
            name = self._normalize_chem_name(entry['name'], is_clean=True)
            cas = entry['cas']
            cid = self.cas_cid_map[cas] if cas in self.cas_cid_map else self.name_cid_map.get(name)
            if cid is None or cid in processed_cids:
                continue

            solubility_parsed = []
            solubility_raw = entry.get('solubility', None)
            if not solubility_raw:
                continue

            for solvent, sol_cat in solubility_raw.items():
                norm_solvent = self._normalize_chem_name(solvent, is_clean=True)
                if norm_solvent not in self.name_cid_map:
                    continue
                solvent_cid = self.name_cid_map[norm_solvent]
                if sol_cat in self.solubility_cats:
                    solubility_parsed.append({"solvent": solvent, "cid": solvent_cid, 'solubility': sol_cat})

            entry['solubility'] = solubility_parsed
            entry['cid'] = cid

            processed_cids.add(cid)
            parsed_entries.append(entry)
        
        self._write_jsonl(parsed_entries, output_fn)



    def _clean_flammability(self, input_fn, output_fn=None):
        if output_fn is None:
            output_fn = input_fn

        entries = self._load_jsonl(input_fn)
        
        for entry in entries:
            if entry['flash_point']:
                entry['flash_point'] = self.__clean_flash_point(entry['flash_point'])
            
            if entry['flash_limits']:
                entry['flash_limits'] = self.__clean_flash_limits(entry['flash_limits'])
            
            if entry['ignition_temp']:
                entry['ignition_temp'] = self.__clean_ignition_temp(entry['ignition_temp'])

            entry['name'] = self.__clean_unicode(entry['name'])
        

        self._write_jsonl(entries, output_fn)
    

    def _map_flammability(self, input_fn, output_fn=None):
        if output_fn is None:
            output_fn = input_fn

        entries = self._load_jsonl(input_fn)

        parsed_entries = []
        processed_cids = set()
        for entry in entries:
            name = self._normalize_chem_name(entry['name'], is_clean=True)
            cid = self.name_cid_map.get(name)
            if cid is None or cid in processed_cids:
                continue

            entry['cid'] = cid

            processed_cids.add(cid)
            parsed_entries.append(entry)
        
        self._write_jsonl(parsed_entries, output_fn)
    

    def parse_crc(self):
        with NamedTemporaryFile(suffix='.jsonl', delete=True) as tmp, self.no_warnings():
            tmp_fn = tmp.name
            self._parse_inorganic_constants_raw(tmp_fn)
            self._clean_inorganic_constants(tmp_fn)
            self._map_inorganic_organic_constants(tmp_fn, self.crc_inorganic_constants_fn)

            self._parse_organic_constants_raw(tmp_fn)
            self._clean_organic_constants(tmp_fn)
            self._map_inorganic_organic_constants(tmp_fn, self.crc_organic_constants_fn)

            self._parse_flammability_raw(tmp_fn)
            self._clean_flammability(tmp_fn)
            self._map_flammability(tmp_fn, self.crc_flammability_fn)
    

    def map_crc_chems_to_cids(self):
        names = [x['name'] for x in self._load_jsonl(self.crc_inorganic_constants_fn)]
        names += [x['name'] for x in self._load_jsonl(self.crc_organic_constants_fn)]

        with open(self.crc_unmapped_names_fn, 'w') as f:
            for name in names:
                norm_name = self._normalize_chem_name(name, is_clean=True)
                if norm_name not in self.name_cid_map:
                    f.write(f"{norm_name}||{name}||0" + '\n')


if __name__ == "__main__":
    crc = ChemsCRC('data/')
    crc.parse_crc()