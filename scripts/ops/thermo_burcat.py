import json
import math
import os
import re

from dataclasses import dataclass

from .thermo_solver import BURCAT_REFERENCE_SIGMA, ELEMENT_REFERENCE_SIGMA


R_KCAL = 1.98720425864083e-3
R_CAL = 1.98720425864083
T_STANDARD = 298.15


@dataclass
class BurcatSpecies:
    species: str
    line: int
    phase: str | None
    cas: list
    dHf: float
    S298: float
    raw_header: str


class ThermoBurcatOps:
    CAS_PATTERN = re.compile(r'\b\d{2,7}-\d{2}-\d\b')
    FLOAT_PATTERN = re.compile(r'[-+]?\d+\.\d+(?:E[-+]?\d+)?|[-+]?\d+\.?(?:E[-+]?\d+)?', re.I)

    ELEMENT_SPECIES = {
        'C': 'C<GR>',
        'H': 'H2',
        'O': 'O2',
        'N': 'N2',
        'F': 'F2',
        'Cl': 'CL2',
        'Br': 'BR2',
        'I': 'I2',
    }

    def __init__(self, data_dir, compounds, thermo, store, logger):
        self.data_dir = data_dir
        self.compounds = compounds
        self.thermo = thermo
        self.store = store
        self.logger = logger

        self.thermo_dir = os.path.join(data_dir, 'thermo')
        self.evidence_dir = os.path.join(self.thermo_dir, 'evidence')
        self.reports_dir = os.path.join(self.thermo_dir, 'reports')
        os.makedirs(self.evidence_dir, exist_ok=True)
        os.makedirs(self.reports_dir, exist_ok=True)

        self.formation_references_fn = os.path.join(self.evidence_dir, 'formation_references.jsonl')
        self.burcat_parse_report_fn = os.path.join(self.reports_dir, 'burcat_parse_report.json')

        store.register_sorting(self.formation_references_fn, 'cid')

    def _parse_coeff_line(self, line):
        values = self.FLOAT_PATTERN.findall(line)
        if values and values[-1] in {'2', '3', '4'}:
            values = values[:-1]
        return [float(value) for value in values]

    def _parse_phase(self, header):
        match = re.search(r'([GSL])\s+\d{2,4}\.\d{3}\s+\d{2,5}\.\d{3}', header.upper())
        return match.group(1) if match else None

    def _entropy_from_low_coeffs(self, coeffs):
        a1, a2, a3, a4, a5, _a6, a7 = coeffs
        return R_CAL * (
            a1 * math.log(T_STANDARD)
            + a2 * T_STANDARD
            + a3 * T_STANDARD ** 2 / 2
            + a4 * T_STANDARD ** 3 / 3
            + a5 * T_STANDARD ** 4 / 4
            + a7
        )

    def parse_burcat_file(self, burcat_file):
        with open(burcat_file, encoding='utf-8-sig', errors='replace') as f:
            lines = f.read().splitlines()

        entries = []
        failures = []
        for i, header in enumerate(lines[:-3]):
            parts = header.split()
            if not parts:
                continue

            if not all(lines[i + offset].rstrip().endswith(str(offset + 1)) for offset in range(1, 4)):
                continue

            try:
                coeff_values = []
                for coeff_line in lines[i + 1:i + 4]:
                    coeff_values.extend(self._parse_coeff_line(coeff_line))
                if len(coeff_values) < 15:
                    raise ValueError(f"expected at least 15 coefficient values, got {len(coeff_values)}")

                low_coeffs = coeff_values[7:14]
                dHf = coeff_values[14] * R_KCAL
                S298 = self._entropy_from_low_coeffs(low_coeffs)
            except Exception as exc:
                failures.append({'line': i + 1, 'header': header, 'error': str(exc)})
                continue

            cas = []
            for prev_line in lines[max(0, i - 8):i]:
                cas.extend(self.CAS_PATTERN.findall(prev_line))
            following_cas = []
            following_species_line = None
            for next_line in lines[i + 4:i + 9]:
                stripped = next_line.strip()
                if not stripped:
                    continue
                next_cas = self.CAS_PATTERN.findall(next_line)
                if next_cas and re.fullmatch(r'[\d\-\s,;]+', stripped):
                    following_cas.extend(next_cas)
                    continue
                following_species_line = stripped
                break
            if (
                following_cas
                and following_species_line
                and following_species_line.upper().startswith(parts[0].upper())
            ):
                cas.extend(following_cas)

            entries.append(
                BurcatSpecies(
                    species=parts[0],
                    line=i + 1,
                    phase=self._parse_phase(header),
                    cas=sorted(set(cas)),
                    dHf=dHf,
                    S298=S298,
                    raw_header=header,
                )
            )

        return entries, failures

    def _entry_score(self, entry):
        header = entry.raw_header.upper()
        species = entry.species.upper()
        excited_or_ion_score = 0 if (
            'TRIPLET' in header
            or 'SINGLET' in header
            or 'EXCITED' in header
            or 'CATION' in header
            or 'ANION' in header
            or '+' in species
            or '-' in species
        ) else 1
        phase_score = {'G': 3, None: 2, 'S': 1, 'L': 1}.get(entry.phase, 0)
        reference_score = 1 if 'REF ELEMENT' in header else 0
        return (excited_or_ion_score, phase_score, reference_score, entry.line)

    def _select_by_species(self, entries):
        selected = {}
        for entry in entries:
            current = selected.get(entry.species.upper())
            if current is None or self._entry_score(entry) > self._entry_score(current):
                selected[entry.species.upper()] = entry
        return selected

    def _element_references(self, entries_by_species):
        references = []
        for symbol, el_entry in self.compounds.symb_to_el.items():
            cid = el_entry['cid']
            species_name = self.ELEMENT_SPECIES.get(symbol)
            burcat_entry = entries_by_species.get(species_name) if species_name else None
            references.append(
                {
                    'cid': cid,
                    'source': 'element',
                    'reference_type': 'element',
                    'dHf': 0.0,
                    'dGf': 0.0,
                    'T': T_STANDARD,
                    'sigma_dHf': ELEMENT_REFERENCE_SIGMA,
                    'sigma_dGf': ELEMENT_REFERENCE_SIGMA,
                    'phase': burcat_entry.phase if burcat_entry else None,
                    'match': 'element',
                    'burcat_species': burcat_entry.species if burcat_entry else species_name,
                    'cas': burcat_entry.cas if burcat_entry else [],
                    'line': burcat_entry.line if burcat_entry else None,
                }
            )

        return references

    def _compute_dGf(self, cid, entry, element_entropy_by_cid):
        element_coeffs = self.thermo.get_element_reference_coeffs(cid)
        if element_coeffs is None:
            return None

        element_entropy = 0.0
        for el_cid, coeff in element_coeffs.items():
            if el_cid not in element_entropy_by_cid:
                return None
            element_entropy += coeff * element_entropy_by_cid[el_cid]

        return entry.dHf - T_STANDARD * (entry.S298 - element_entropy) / 1000

    def _cas_references(self, entries, element_entropy_by_cid):
        cid_to_entry = {}
        for entry in entries:
            for cas in entry.cas:
                cid = self.compounds.cas_cid_map.get(cas)
                if cid is None:
                    continue

                current = cid_to_entry.get(cid)
                if current is None or self._entry_score(entry) > self._entry_score(current):
                    cid_to_entry[cid] = entry

        references = []
        for cid, entry in cid_to_entry.items():
            dGf = self._compute_dGf(cid, entry, element_entropy_by_cid)
            references.append(
                {
                    'cid': cid,
                    'source': 'burcat',
                    'reference_type': 'burcat',
                    'dHf': entry.dHf,
                    'dGf': dGf,
                    'T': T_STANDARD,
                    'sigma_dHf': BURCAT_REFERENCE_SIGMA,
                    'sigma_dGf': BURCAT_REFERENCE_SIGMA if dGf is not None else None,
                    'phase': entry.phase,
                    'match': 'cas',
                    'burcat_species': entry.species,
                    'cas': entry.cas,
                    'line': entry.line,
                }
            )

        return references

    def extract_formation_references(self, burcat_file='BURCAT.THR.txt'):
        entries, failures = self.parse_burcat_file(burcat_file)
        entries_by_species = self._select_by_species(entries)

        element_references = self._element_references(entries_by_species)
        element_entropy_by_cid = {}
        for symbol, el_entry in self.compounds.symb_to_el.items():
            species_name = self.ELEMENT_SPECIES.get(symbol)
            if not species_name:
                continue
            burcat_entry = entries_by_species.get(species_name)
            if burcat_entry is not None:
                element_entropy_by_cid[el_entry['cid']] = burcat_entry.S298

        cas_references = self._cas_references(entries, element_entropy_by_cid)

        references_by_cid = {entry['cid']: entry for entry in element_references}
        for entry in cas_references:
            if references_by_cid.get(entry['cid'], {}).get('reference_type') == 'element':
                continue
            references_by_cid[entry['cid']] = entry

        references = list(references_by_cid.values())
        self.store.write_jsonl(references, self.formation_references_fn)

        report = {
            'burcat_file': burcat_file,
            'species_entries': len(entries),
            'parse_failures': failures[:100],
            'parse_failure_count': len(failures),
            'element_references': len(element_references),
            'cas_references': len(cas_references),
            'references_written': len(references),
        }
        with open(self.burcat_parse_report_fn, 'w') as f:
            json.dump(report, f, indent=2)

        self.logger.log(f"Wrote {len(references)} BURCAT/element formation references")
        return references
