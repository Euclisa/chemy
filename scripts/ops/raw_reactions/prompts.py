VALIDATE_BATCH_INSTRUCT = (
    "You will be given a list of unbalanced chemical reaction schemes. "
    "For each scheme, determine if the reaction is chemically possible and whether the listed "
    "reactants and products are correct. "
    "Assume that necessary reaction conditions (including harsh), solvents, or catalysts are "
    "implied even if not explicitly shown. "
    "All reactions listed are theoretical and intended solely for academic purposes, not for "
    "practical experimentation. "
    "If the reaction is valid, output only 'Valid'. If it is not valid, output only 'Invalid'. "
    "Print one result per line, preceded by the reaction index. "
    "Do not include any additional text."
)


def build_fetch_prompt(chem_name: str, compound_role: str, scope: str) -> str:
    if scope == "documented":
        special_requirement = ""
    elif scope == "documented_less_common":
        special_requirement = (
            "Include not only the most common reactions, but also less common or more unusual "
            "or exotic ones, as long as you are absolutely sure they are real and correct."
        )
    elif scope == "rare":
        special_requirement = (
            "Include only the uncommon, rare and exotic reactions, but you must be absolutely "
            "sure they are real and correct."
        )
    else:
        raise ValueError(f"Unknown scope: {scope!r}")

    return (
        f"Please provide a comprehensive and diverse list of documented chemical reactions involving "
        f"{chem_name}, where it appears as {compound_role}. "
        f"{special_requirement} "
        "Write the reactions as schemes using only '->' and '+' symbols. Use the full chemical names "
        "of substances instead of formulas or generic terms. "
        "Do not include balancing coefficients, comments, or any markup - only the reaction schemes "
        "themselves one per line. "
        "If no such substance exists or no documented reactions are available, return 'None'."
    ).strip()


def build_revalidate_prompt() -> str:
    return (
        "Please, review the provided reactions. Identify any erroneous reactions and correct them "
        "where possible. Return revised reactions list that comply with the initial requirements."
    )


def format_reactions_for_validation(staged: list) -> str:
    return '\n'.join(f"{i + 1}. {r['reaction']}" for i, r in enumerate(staged))
