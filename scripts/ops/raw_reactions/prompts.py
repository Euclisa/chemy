import json


VALID_PHASES = frozenset({'g', 'l', 's', 'aq'})

_POSITION_LABELS = {
    "any": "a reagent or product",
    "reagent": "a reagent",
    "product": "a product",
}

_SCOPE_REQUIREMENTS = {
    "documented": "",
    "documented_less_common": (
        "Include not only the most common reactions, but also less common or unusual "
        "ones, as long as you are absolutely sure they are real and correct. "
    ),
    "rare": (
        "Include only uncommon, rare, and exotic reactions, but you must be absolutely "
        "sure they are real and correct. "
    ),
}

_JSONL_SCHEMA = (
    'Each JSON object has these fields:\n'
    '- "reagents": list of {"name": "<clean chemical name>", "phase": "<g|l|s|aq>"}\n'
    '- "products": list of {"name": "<clean chemical name>", "phase": "<g|l|s|aq>"}\n'
    '- Optional "solvent": primary solvent name or null\n'
    '- Optional "catalyst": catalyst name or null\n'
    '- Optional per-compound "note": short qualifier; omit it when not needed'
)


def build_fetch_prompt(chem_name, position, scope, existing_reactions=None):
    if scope not in _SCOPE_REQUIREMENTS:
        raise ValueError(f"Unknown scope: {scope!r}")
    if position not in _POSITION_LABELS:
        raise ValueError(f"Unknown position: {position!r}")

    position_label = _POSITION_LABELS[position]
    context = ""
    if existing_reactions:
        lines = "\n".join(f"  {i + 1}. {r}" for i, r in enumerate(existing_reactions))
        context = (
            f"We already know these reactions where {chem_name} appears as "
            f"{position_label}:\n{lines}\n\n"
            f"Generate ADDITIONAL reactions where {chem_name} appears as "
            f"{position_label} and that are not in this list.\n"
        )

    return (
        f"{context}"
        f"Provide a list of documented chemical reactions involving {chem_name}, "
        f"where it appears as {position_label}. "
        f"{_SCOPE_REQUIREMENTS[scope]}"
        f"Return one JSON object per line (JSONL). {_JSONL_SCHEMA}\n\n"
        "Rules:\n"
        "- Use full chemical names, not formulas.\n"
        "- The name field must contain ONLY the chemical name.\n"
        "- Put qualifiers such as concentrated, dilute, powdered, heated, "
        "bubbled through solution, or acidic solution in note, not in name.\n"
        "- Omit note entirely when there is no useful qualifier.\n"
        "- Do not include stoichiometric coefficients.\n"
        "- Do not include catalysts or solvents in the reagents or products lists.\n"
        "- Every reagent and product must include phase.\n"
        "- Phases: g = gas, l = liquid, s = solid, aq = aqueous solution.\n"
        "- If no reactions are available, return a single line: null\n"
        "- Return ONLY valid JSONL lines, no extra text."
    )


def build_revalidate_prompt(reactions_jsonl_str):
    return (
        "Review the following chemical reactions for correctness. For each, verify:\n"
        "1. The reaction is chemically valid and documented.\n"
        "2. Reagents and products are correctly identified.\n"
        "3. Every reagent and product has a reasonable phase annotation "
        "(g, l, s, aq) for typical stated or implicit reaction conditions.\n"
        "4. Solvent and catalyst (if listed) are appropriate.\n\n"
        f"Return the corrected list in the same JSONL format. {_JSONL_SCHEMA}\n\n"
        "Remove chemically invalid reactions. Fix incorrect phases, solvents, catalysts, "
        "or misplaced name qualifiers.\n"
        "Return ONLY valid JSONL lines, no extra text.\n\n"
        f"{reactions_jsonl_str}"
    )


VALIDATE_BATCH_INSTRUCT = (
    "You will be given a list of chemical reactions with phase annotations. "
    "For each reaction, determine if it is chemically valid/documented and if every "
    "specified phase is plausible for the stated or typical implicit reaction conditions.\n"
    "Output 'Valid' only if both the chemistry and phases are correct. "
    "Output 'Invalid' if the chemistry is wrong or any specified phase is unreasonable.\n"
    "Print one result per line, preceded by the reaction index.\n"
    "Do not include any additional text."
)

REPAIR_BATCH_INSTRUCT = (
    "You will be given chemical reactions that received low but nonzero validation "
    "confidence. For each input, either correct the reaction or return null if it "
    "cannot be confidently repaired.\n"
    "Only fix chemically meaningful issues such as incorrect phase, misplaced "
    "qualifiers, incorrect solvent/catalyst placement, or obvious participant naming "
    "errors. Do not invent a different reaction just to make the row valid.\n"
    "Return exactly one line per input reaction, in the same order. Each output line "
    "must be either null or a corrected reaction JSON object.\n"
    f"Corrected reaction objects must follow this schema. {_JSONL_SCHEMA}\n"
    "Return ONLY valid JSONL lines, no extra text."
)


def format_structured_reaction(reaction):
    def fmt_compounds(compounds):
        parts = []
        for compound in compounds:
            part = f"{compound['name']}({compound['phase']})"
            if compound.get('note'):
                part += f" {{{compound['note']}}}"
            parts.append(part)
        return " + ".join(parts)

    s = f"{fmt_compounds(reaction['reagents'])} -> {fmt_compounds(reaction['products'])}"
    extras = []
    if reaction.get('solvent'):
        extras.append(f"solvent: {reaction['solvent']}")
    if reaction.get('catalyst'):
        extras.append(f"catalyst: {reaction['catalyst']}")
    if extras:
        s += f" [{'; '.join(extras)}]"
    return s


def format_reactions_for_validation(staged):
    return '\n'.join(
        f"{i + 1}. {format_structured_reaction(r['reaction'])}"
        for i, r in enumerate(staged)
    )


def format_repair_candidates(candidates):
    return '\n'.join(
        json.dumps(candidate['reaction'])
        for candidate in candidates
    )
