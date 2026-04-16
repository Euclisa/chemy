import os
from dataclasses import dataclass
from typing import Callable


@dataclass
class RawReactionsPreset:
    name: str
    compound_role: str   # "reactant or product" | "product"
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


PRESETS = {
    "default": RawReactionsPreset(
        name="default",
        compound_role="reactant or product",
        scope="documented",
        build_criteria=lambda compounds: _mapped,
    ),
    "wiki_crc_rp": RawReactionsPreset(
        name="wiki_crc_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=_build_wiki_crc_criteria,
    ),
    "simple_inorganic_wiki_rp": RawReactionsPreset(
        name="simple_inorganic_wiki_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=_build_simple_inorganic_wiki_criteria,
    ),
    "simple_organic_wiki_rp": RawReactionsPreset(
        name="simple_organic_wiki_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=_build_simple_organic_wiki_criteria,
    ),
    "oxidizers_wiki_rp": RawReactionsPreset(
        name="oxidizers_wiki_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=_build_wiki_hazard_criteria({'GHS03'}),
    ),
    "corrosive_wiki_rp": RawReactionsPreset(
        name="corrosive_wiki_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=_build_wiki_hazard_criteria({'GHS05'}),
    ),
    "wiki_uncommon_rp": RawReactionsPreset(
        name="wiki_uncommon_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=lambda compounds: lambda chem: _mapped(chem) and _has_wiki(chem),
    ),
    "top_rare_rp": RawReactionsPreset(
        name="top_rare_rp",
        compound_role="reactant or product",
        scope="rare",
        build_criteria=_build_dangerous_criteria,
    ),
    "wiki_uncommon_p": RawReactionsPreset(
        name="wiki_uncommon_p",
        compound_role="product",
        scope="documented_less_common",
        build_criteria=lambda compounds: lambda chem: _mapped(chem) and _has_wiki(chem),
    ),
    "annotated_uncommon_p": RawReactionsPreset(
        name="annotated_uncommon_p",
        compound_role="product",
        scope="documented_less_common",
        build_criteria=lambda compounds: lambda chem: (
            _mapped(chem) and _has_annotation(chem) and not _has_wiki(chem)
        ),
    ),
}


def get_preset(name: str) -> RawReactionsPreset:
    if name not in PRESETS:
        raise ValueError(f"Unknown preset '{name}'. Available: {list(PRESETS)}")
    return PRESETS[name]
