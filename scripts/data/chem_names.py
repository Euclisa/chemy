from dataclasses import dataclass
import re


UNKNOWN_NAME = "<Unknown>"
CAS_PATTERN = r'\d{2,7}-\d{2}-\d'


@dataclass(frozen=True)
class NameFilterResult:
    accepted: bool
    reason: str | None
    clean_name: str | None
    norm_name: str | None


UNICODE_MAP = {
    '\u2011': '-',
    '\u2013': '-',
    '\u2019': "'",
    '\u03b1': 'alpha',
    '\u0391': 'alpha',
    '\u03b2': 'beta',
    '\u0392': 'beta',
    '\u03b3': 'gamma',
    '\u0393': 'gamma',
    '\u03b4': 'delta',
    '\u0394': 'delta',
}


SYNONYM_TAGS = [
    'USP', 'EP', 'HSDB', 'EP MONOGRAPH', 'WHO-DD', 'USAN', 'MI',
    'USP MONOGRAPH', 'VANDF', 'CZECH', 'ACGIH', 'NDIPA', 'German',
    'VAN', 'INN', '9CI', '8CI', 'JAN', 'French', 'Latin', 'ISO', 'Standard',
    'natural', 'HPUS', 'IARC', 'NF', 'INCI', 'TN', 'Dutch', 'Italian',
    'FCC', 'DCIT', 'BAN', 'ORANGE BOOK', 'Spanish', 'IUPAC', 'FHFI', 'Polish',
]


DISCARD_WORDS_PART = [
    'oil', 'solid', 'solids', 'liquid', 'dry', 'powder', 'nanopowder',
    'beads', 'impurity', 'grade', 'intermediate', 'title', 'desired',
    'material', 'solution', 'syrup', 'crystals', 'residue', 'compound',
    'product', 'titled', 'mixture', 'foam', 'needles', 'crude', 'resin',
    'gum',
]


DISCARD_WORDS_WHOLE = [
    'acetate', 'acid', 'salt', 'phosphonate', 'ester', 'amide', 'nitrile',
    'aldehyde', 'amine', 'alcohol', 'ketone', 'diacid', 'gas', 'cyano',
    'azide', 'metal', 'ether', 'oxime', 'imine', 'alkyne', 'alkene',
    'hydrochloride', 'i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii',
    'ix', 'x',
]


DISCARD_WORD_PART_PATTERN = (
    '(' + '|'.join([r'\b' + re.escape(word) + r'\b' for word in DISCARD_WORDS_PART]) + ')'
)
DISCARD_WORD_WHOLE_PATTERN = (
    '(' + '|'.join([r'^' + re.escape(word) + r'$' for word in DISCARD_WORDS_WHOLE]) + ')'
)


BAD_PATTERN_RULES = [
    r'[:=%<>@/\\_.#&*";?!]',
    r'[-,.]$',
    r'^[-,.]',
    r'\(\d:\d\)',
    r'\s{2,}',
    r'\b[nm]m\b',
    r'-,',
    r'^[a-z]?[0-9-()\s]+[a-z]?$',
    r'\(\s+.+\s+\)',
]


def _has_unsupported_non_ascii(value):
    return any((not char.isascii()) and char not in UNICODE_MAP for char in str(value))


def chem_name_to_ascii(chem_name_raw):
    chem_name_ascii = ''
    for char in str(chem_name_raw):
        if not char.isascii():
            char = UNICODE_MAP.get(char, '')
        chem_name_ascii += char
    return chem_name_ascii


def clean_chem_name(chem_name_raw, is_clean=False, unknown_name=UNKNOWN_NAME):
    if chem_name_raw == unknown_name:
        return unknown_name

    chem_name = chem_name_to_ascii(chem_name_raw)
    chem_name = chem_name.strip()
    chem_name = re.sub(r'\s+', ' ', chem_name)

    if not is_clean:
        chem_name = chem_name.strip('`\'".,;:')
        chem_name = re.sub(r'^\d+ ', '', chem_name)

    return chem_name


def normalize_chem_name(chem_name_raw, is_clean=False, unknown_name=UNKNOWN_NAME):
    if chem_name_raw == unknown_name:
        return unknown_name

    chem_name = clean_chem_name(chem_name_raw, is_clean=is_clean, unknown_name=unknown_name)
    chem_name = chem_name.lower()
    chem_name = chem_name.strip()
    chem_name = chem_name.replace('aluminum', 'aluminium')

    if not is_clean:
        chem_name = re.sub(r' \([^\d]+\)$', '', chem_name)
        chem_name = chem_name.replace(' vapor', '')
        chem_name = chem_name.replace(' dust', '')
        chem_name = chem_name.replace('solution', '')
        chem_name = chem_name.replace('concentrated', '')
        chem_name = chem_name.replace('dilute ', '')
        chem_name = chem_name.replace('fuming ', '')
        chem_name = chem_name.replace('solid', '')
        chem_name = chem_name.replace('glacial ', '')
        chem_name = chem_name.replace('elemental', '')
        chem_name = chem_name.replace(' metal', '')
        chem_name = chem_name.replace('aqueous', '')
        chem_name = chem_name.replace(' gas', '')
        chem_name = chem_name.replace('hot ', '')
        chem_name = chem_name.replace('uv light', 'light')
        chem_name = chem_name.replace('blue light', 'light')
        chem_name = chem_name.replace('ultraviolet light', 'light')

        if 'catalyst' in chem_name or 'raney nickel' in chem_name:
            chem_name = 'catalyst'

    return re.sub(r'\s+', '', chem_name)


def clean_synonym(name, unknown_name=UNKNOWN_NAME):
    if name == unknown_name:
        return unknown_name

    name = clean_chem_name(name, is_clean=True, unknown_name=unknown_name)
    pattern = r'[\(\[]\s*(?:' + '|'.join(re.escape(tag) for tag in SYNONYM_TAGS) + r')\s*[\)\]]$'
    name = re.sub(pattern, '', name, flags=re.IGNORECASE).strip()
    name = re.sub(r'\bfume\b', '', name, flags=re.IGNORECASE).strip()
    return name


def classify_name(name, manual_filtered=None, unknown_name=UNKNOWN_NAME):
    manual_filtered = set() if manual_filtered is None else set(manual_filtered)

    if name == unknown_name:
        return NameFilterResult(True, None, unknown_name, unknown_name)

    if name is None:
        return NameFilterResult(False, 'empty', None, None)

    clean_name = clean_synonym(name, unknown_name=unknown_name)
    norm_name = normalize_chem_name(clean_name, is_clean=True, unknown_name=unknown_name)

    if norm_name in manual_filtered:
        return NameFilterResult(False, 'manual_decision', clean_name, norm_name)

    if not clean_name or not norm_name:
        return NameFilterResult(False, 'empty', clean_name, norm_name)

    if _has_unsupported_non_ascii(name):
        return NameFilterResult(False, 'non_ascii', clean_name, norm_name)

    if len(clean_name) <= 1:
        return NameFilterResult(False, 'generic_word', clean_name, norm_name)

    name_lower = clean_name.strip().lower()

    if re.search(CAS_PATTERN, name_lower):
        return NameFilterResult(False, 'cas', clean_name, norm_name)

    if (
        'unii-' in name_lower
        or re.fullmatch(r'[a-z0-9]{10}', name_lower)
        and re.search(r'[abdefgijklmqrtuvwxyz]', name_lower)
        and re.search(r'\d', name_lower)
        or re.search(r'[a-z]{14}-[a-z]{10}-[a-z]', name_lower)
        or re.search(r'\d{3,}', name_lower)
        or re.fullmatch(r'\d+', name_lower)
    ):
        return NameFilterResult(False, 'identifier', clean_name, norm_name)

    if re.search(DISCARD_WORD_PART_PATTERN, name_lower) or re.search(
        DISCARD_WORD_WHOLE_PATTERN,
        name_lower,
    ):
        return NameFilterResult(False, 'generic_word', clean_name, norm_name)

    if any(re.search(pattern, name_lower) for pattern in BAD_PATTERN_RULES):
        return NameFilterResult(False, 'bad_pattern', clean_name, norm_name)

    return NameFilterResult(True, None, clean_name, norm_name)
