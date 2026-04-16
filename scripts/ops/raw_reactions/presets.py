import os
from dataclasses import dataclass
from typing import Callable


POSITION_ANY = 'any'
POSITION_REAGENT = 'reagent'
POSITION_PRODUCT = 'product'

POSITION_SUFFIXES = {
    'rp': POSITION_ANY,
    'r': POSITION_REAGENT,
    'p': POSITION_PRODUCT,
}


@dataclass(frozen=True)
class RawReactionsBasePreset:
    name: str
    scope: str           # "documented" | "documented_less_common" | "rare"
    build_criteria: Callable  # (compounds) -> Callable[[dict], bool]


@dataclass(frozen=True)
class RawReactionsPreset:
    name: str
    base_name: str
    position: str        # "any" | "reagent" | "product"
    scope: str           # "documented" | "documented_less_common" | "rare"
    build_criteria: Callable  # (compounds) -> Callable[[dict], bool]


def _build_dangerous_criteria(compounds):
    return _build_wiki_hazard_criteria({'GHS01', 'GHS03', 'GHS06'})(compounds)


def _mapped(chem):
    return chem['cid'] > 0


def _has_wiki(chem):
    return bool(chem.get('wiki'))


def _has_annotation(chem):
    return bool(chem.get('annotation'))


def _load_cids(store, filename):
    return frozenset(entry['cid'] for entry in store.load_jsonl(filename))


def _load_wiki_hazard_cids(compounds, pictograms):
    pictograms = set(pictograms)
    hazard_cids = set()
    for hazard in compounds.store.load_jsonl(compounds.chems_hazards_wiki_fn):
        if pictograms & set(hazard.get('pictograms') or []):
            hazard_cids.add(hazard['cid'])
    return frozenset(hazard_cids)


def _build_wiki_hazard_criteria(pictograms):
    def build(compounds):
        hazard_cids = _load_wiki_hazard_cids(compounds, pictograms)
        return lambda chem: _mapped(chem) and _has_wiki(chem) and chem['cid'] in hazard_cids

    return build


def _build_crc_cids(compounds):
    crc_dir = os.path.join(compounds.data_dir, 'chems_properties', 'crc_handbook')
    crc_inorganic = _load_cids(compounds.store, os.path.join(crc_dir, 'inorganic_constants.jsonl'))
    crc_organic = _load_cids(compounds.store, os.path.join(crc_dir, 'organic_constants.jsonl'))
    return crc_inorganic | crc_organic


def _build_wiki_crc_criteria(compounds):
    crc_cids = _build_crc_cids(compounds)
    return lambda chem: _mapped(chem) and _has_wiki(chem) and chem['cid'] in crc_cids


def _build_simple_inorganic_wiki_criteria(compounds):
    return lambda chem: (
        _mapped(chem)
        and _has_wiki(chem)
        and not chem.get('organic')
        and chem.get('heavy_count', 0) <= 12
    )


def _build_simple_organic_wiki_criteria(compounds):
    return lambda chem: (
        _mapped(chem)
        and _has_wiki(chem)
        and bool(chem.get('organic'))
        and chem.get('heavy_count', 0) <= 12
    )


BASE_PRESETS = {
    'default': RawReactionsBasePreset(
        name='default',
        scope='documented',
        build_criteria=lambda compounds: _mapped,
    ),
    'wiki_crc': RawReactionsBasePreset(
        name='wiki_crc',
        scope='documented_less_common',
        build_criteria=_build_wiki_crc_criteria,
    ),
    'simple_inorganic_wiki': RawReactionsBasePreset(
        name='simple_inorganic_wiki',
        scope='documented_less_common',
        build_criteria=_build_simple_inorganic_wiki_criteria,
    ),
    'simple_organic_wiki': RawReactionsBasePreset(
        name='simple_organic_wiki',
        scope='documented_less_common',
        build_criteria=_build_simple_organic_wiki_criteria,
    ),
    'oxidizers_wiki': RawReactionsBasePreset(
        name='oxidizers_wiki',
        scope='documented_less_common',
        build_criteria=_build_wiki_hazard_criteria({'GHS03'}),
    ),
    'corrosive_wiki': RawReactionsBasePreset(
        name='corrosive_wiki',
        scope='documented_less_common',
        build_criteria=_build_wiki_hazard_criteria({'GHS05'}),
    ),
    'wiki_uncommon': RawReactionsBasePreset(
        name='wiki_uncommon',
        scope='documented_less_common',
        build_criteria=lambda compounds: lambda chem: _mapped(chem) and _has_wiki(chem),
    ),
    'top_rare': RawReactionsBasePreset(
        name='top_rare',
        scope='rare',
        build_criteria=_build_dangerous_criteria,
    ),
    'annotated_uncommon': RawReactionsBasePreset(
        name='annotated_uncommon',
        scope='documented_less_common',
        build_criteria=lambda compounds: lambda chem: (
            _mapped(chem) and _has_annotation(chem) and not _has_wiki(chem)
        ),
    ),
}

DEFAULT_PRESET_NAMES = (
    'default_rp',
    'wiki_crc_rp',
    'simple_inorganic_wiki_rp',
    'simple_organic_wiki_rp',
    'oxidizers_wiki_rp',
    'corrosive_wiki_rp',
    'wiki_uncommon_rp',
    'top_rare_rp',
    'wiki_uncommon_p',
    'annotated_uncommon_p',
)


def list_base_preset_names():
    return tuple(BASE_PRESETS)


def list_position_suffixes():
    return tuple(POSITION_SUFFIXES)


def _format_preset_error(name):
    bases = ', '.join(list_base_preset_names())
    suffixes = ', '.join(list_position_suffixes())
    return (
        f"Unknown preset '{name}'. Expected '<base>_<suffix>' where "
        f"base is one of [{bases}] and suffix is one of [{suffixes}]."
    )


def get_preset(name: str) -> RawReactionsPreset:
    try:
        base_name, suffix = name.rsplit('_', 1)
    except ValueError:
        raise ValueError(_format_preset_error(name))

    base = BASE_PRESETS.get(base_name)
    position = POSITION_SUFFIXES.get(suffix)
    if base is None or position is None:
        raise ValueError(_format_preset_error(name))

    return RawReactionsPreset(
        name=name,
        base_name=base.name,
        position=position,
        scope=base.scope,
        build_criteria=base.build_criteria,
    )
