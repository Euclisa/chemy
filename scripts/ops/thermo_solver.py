import math

import numpy as np

from dataclasses import dataclass


ELEMENT_REFERENCE_SIGMA = 0.01
BURCAT_REFERENCE_SIGMA = 1.0
LLM_RELATIVE_SIGMA_FRACTION = 0.25
LLM_ABSOLUTE_SIGMA_FLOOR = 15.0
LLM_FALLBACK_SIGMA = 25.0


@dataclass
class ThermoSolveResult:
    formation: list
    corrected_reactions: list
    metadata: dict


def is_finite_number(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def llm_sigma(value, std=None):
    if not is_finite_number(value):
        return LLM_FALLBACK_SIGMA

    candidates = [
        abs(value) * LLM_RELATIVE_SIGMA_FRACTION,
        LLM_ABSOLUTE_SIGMA_FLOOR,
    ]
    if is_finite_number(std):
        candidates.append(abs(std))

    return max(candidates)


def reference_sigma(entry, prop):
    sigma_key = f'sigma_{prop}'
    if is_finite_number(entry.get(sigma_key)) and entry[sigma_key] > 0:
        return entry[sigma_key]

    if entry.get('reference_type') == 'element':
        return ELEMENT_REFERENCE_SIGMA

    return BURCAT_REFERENCE_SIGMA


def reaction_coeff_map(reaction, balance):
    coeffs = {}
    for reagent in reaction['reagents']:
        cid = reagent['cid']
        coeffs[cid] = coeffs.get(cid, 0) - balance[cid]
    for product in reaction['products']:
        cid = product['cid']
        coeffs[cid] = coeffs.get(cid, 0) + balance[cid]
    return coeffs


def compute_from_coeffs(coeffs, cid_to_value):
    result = 0.0
    for cid, coeff in coeffs.items():
        value = cid_to_value.get(cid)
        if value is None:
            return None
        result += coeff * value
    return result


def _add_unknown(cid, cid_to_index):
    if cid not in cid_to_index:
        cid_to_index[cid] = len(cid_to_index)


def _valid_weighted_row(value, sigma):
    return is_finite_number(value) and is_finite_number(sigma) and sigma > 0


def solve_formation_thermo(reactions, reaction_balances, reaction_estimates, formation_references):
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import lsqr

    reaction_by_rid = {reaction['rid']: reaction for reaction in reactions}
    reaction_coeffs = {}
    cid_to_index = {}
    skipped = []

    for reaction in reactions:
        rid = reaction['rid']
        balance = reaction_balances.get(rid)
        if balance is None:
            skipped.append({'kind': 'reaction', 'rid': rid, 'reason': 'missing_balance'})
            continue

        coeffs = reaction_coeff_map(reaction, balance)
        reaction_coeffs[rid] = coeffs
        for cid in coeffs:
            _add_unknown(cid, cid_to_index)

    for ref in formation_references:
        _add_unknown(ref['cid'], cid_to_index)

    def build_rows(prop):
        rows = []

        for estimate in reaction_estimates:
            rid = estimate['rid']
            if rid not in reaction_by_rid:
                skipped.append({'kind': 'reaction_estimate', 'rid': rid, 'reason': 'unknown_reaction'})
                continue
            if rid not in reaction_coeffs:
                skipped.append({'kind': 'reaction_estimate', 'rid': rid, 'reason': 'missing_coeffs'})
                continue

            value = estimate.get(prop)
            sigma = estimate.get(f'sigma_{prop}')
            if sigma is None:
                sigma = llm_sigma(value, estimate.get(f'std_{prop}'))
            if not _valid_weighted_row(value, sigma):
                skipped.append({'kind': 'reaction_estimate', 'rid': rid, 'prop': prop, 'reason': 'invalid_value_or_sigma'})
                continue

            rows.append(
                {
                    'coeffs': reaction_coeffs[rid],
                    'value': value,
                    'sigma': sigma,
                    'source': estimate.get('source'),
                    'kind': 'reaction_estimate',
                }
            )

        for ref in formation_references:
            value = ref.get(f'{prop}f')
            sigma = reference_sigma(ref, f'{prop}f')
            if not _valid_weighted_row(value, sigma):
                skipped.append({'kind': 'formation_reference', 'cid': ref['cid'], 'prop': prop, 'reason': 'invalid_value_or_sigma'})
                continue

            rows.append(
                {
                    'coeffs': {ref['cid']: 1.0},
                    'value': value,
                    'sigma': sigma,
                    'source': ref.get('source'),
                    'kind': 'formation_reference',
                }
            )

        return rows

    def solve_prop(prop):
        rows = build_rows(prop)
        if not rows or not cid_to_index:
            return {}, {
                'rows': len(rows),
                'istop': None,
                'iterations': 0,
                'residual_norm': None,
                'status': 'empty',
            }

        matrix = lil_matrix((len(rows), len(cid_to_index)))
        values = np.zeros(len(rows))

        for row_i, row in enumerate(rows):
            weight = 1.0 / row['sigma']
            for cid, coeff in row['coeffs'].items():
                matrix[row_i, cid_to_index[cid]] = coeff * weight
            values[row_i] = row['value'] * weight

        solved = lsqr(matrix.tocsr(), values)
        vector = solved[0]

        return (
            {cid: vector[cid_i] for cid, cid_i in cid_to_index.items()},
            {
                'rows': len(rows),
                'istop': solved[1],
                'iterations': solved[2],
                'residual_norm': solved[3],
                'status': 'ok',
            },
        )

    cid_to_dHf, dH_meta = solve_prop('dH')
    cid_to_dGf, dG_meta = solve_prop('dG')

    formation = []
    for cid in sorted(cid_to_index):
        formation.append(
            {
                'cid': cid,
                'dHf': cid_to_dHf.get(cid),
                'dGf': cid_to_dGf.get(cid),
                'T': 298.15,
            }
        )

    corrected_reactions = []
    for reaction in reactions:
        rid = reaction['rid']
        coeffs = reaction_coeffs.get(rid)
        if coeffs is None:
            continue
        corrected_reactions.append(
            {
                'rid': rid,
                'dH': compute_from_coeffs(coeffs, cid_to_dHf),
                'dG': compute_from_coeffs(coeffs, cid_to_dGf),
            }
        )

    metadata = {
        'chemicals': len(cid_to_index),
        'reactions': len(reaction_coeffs),
        'reaction_estimates': len(reaction_estimates),
        'formation_references': len(formation_references),
        'dH': dH_meta,
        'dG': dG_meta,
        'skipped': skipped,
        'weights': {
            'element_reference_sigma': ELEMENT_REFERENCE_SIGMA,
            'burcat_reference_sigma': BURCAT_REFERENCE_SIGMA,
            'llm_relative_sigma_fraction': LLM_RELATIVE_SIGMA_FRACTION,
            'llm_absolute_sigma_floor': LLM_ABSOLUTE_SIGMA_FLOOR,
            'llm_fallback_sigma': LLM_FALLBACK_SIGMA,
        },
    }

    return ThermoSolveResult(formation, corrected_reactions, metadata)
