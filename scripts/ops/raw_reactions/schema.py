from .prompts import VALID_PHASES


REACTION_KEYS = frozenset({'reagents', 'products', 'solvent', 'catalyst'})
COMPOUND_KEYS = frozenset({'name', 'phase', 'note'})


def is_valid_compound(compound):
    if not isinstance(compound, dict):
        return False
    if not set(compound).issubset(COMPOUND_KEYS):
        return False
    if not isinstance(compound.get('name'), str) or not compound['name'].strip():
        return False
    if compound.get('phase') not in VALID_PHASES:
        return False
    if 'note' in compound and (
        not isinstance(compound['note'], str) or not compound['note'].strip()
    ):
        return False
    return True


def is_valid_reaction_obj(obj):
    if not isinstance(obj, dict):
        return False
    if not set(obj).issubset(REACTION_KEYS):
        return False
    for key in ('reagents', 'products'):
        compounds = obj.get(key)
        if not isinstance(compounds, list) or not compounds:
            return False
        if not all(is_valid_compound(c) for c in compounds):
            return False
    for key in ('solvent', 'catalyst'):
        if key in obj and obj[key] is not None and (
            not isinstance(obj[key], str) or not obj[key].strip()
        ):
            return False
    return True


def normalize_reaction_obj(obj):
    for key in ('reagents', 'products'):
        for compound in obj[key]:
            compound['name'] = compound['name'].strip()
            if 'note' in compound:
                compound['note'] = compound['note'].strip()
    if not isinstance(obj.get('solvent'), str):
        obj['solvent'] = None
    else:
        obj['solvent'] = obj['solvent'].strip()
    if not isinstance(obj.get('catalyst'), str):
        obj['catalyst'] = None
    else:
        obj['catalyst'] = obj['catalyst'].strip()
    return obj
