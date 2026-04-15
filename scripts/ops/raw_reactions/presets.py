import json
from dataclasses import dataclass
from typing import Callable


@dataclass
class RawReactionsPreset:
    name: str
    compound_role: str   # "reactant or product" | "product"
    scope: str           # "documented" | "documented_less_common" | "rare"
    build_criteria: Callable  # (compounds) -> Callable[[dict], bool]


def _build_dangerous_criteria(compounds):
    hazards_chems = {}
    with open(compounds.chems_hazards_wiki_fn) as f:
        for line in f:
            hazard = json.loads(line)
            hazards_chems[hazard['cid']] = hazard

    def criteria(chem):
        cid = chem['cid']
        if cid not in hazards_chems:
            return False
        return any(pic in {'GHS01', 'GHS03', 'GHS06'} for pic in hazards_chems[cid]['pictograms'])

    return criteria


PRESETS = {
    "default": RawReactionsPreset(
        name="default",
        compound_role="reactant or product",
        scope="documented",
        build_criteria=lambda compounds: lambda _: True,
    ),
    "wiki_uncommon_rp": RawReactionsPreset(
        name="wiki_uncommon_rp",
        compound_role="reactant or product",
        scope="documented_less_common",
        build_criteria=lambda compounds: lambda chem: chem['wiki'],
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
        build_criteria=lambda compounds: lambda chem: chem['wiki'],
    ),
    "annotated_uncommon_p": RawReactionsPreset(
        name="annotated_uncommon_p",
        compound_role="product",
        scope="documented_less_common",
        build_criteria=lambda compounds: lambda chem: 'annotation' in chem and not chem['wiki'],
    ),
}


def get_preset(name: str) -> RawReactionsPreset:
    if name not in PRESETS:
        raise ValueError(f"Unknown preset '{name}'. Available: {list(PRESETS)}")
    return PRESETS[name]
