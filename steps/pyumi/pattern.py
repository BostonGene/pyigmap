import math
import re

from logger import set_logger

logger = set_logger(name=__file__)

NORMAL_NUCLEOTIDES = 'ATGC'
IUPAC_WILDCARDS = 'RYSWKMBDHVN'
ALLOWED_LETTERS_IN_UMI = NORMAL_NUCLEOTIDES + 'N'

ADAPTER_PATTERN_REGEX = rf"(?<!\[)\b[{ALLOWED_LETTERS_IN_UMI}]+\b(?!\])"


class ValidationError(Exception):
    pass


def parse_umi_length(pattern: str, barcode_type='UMI') -> int:
    """Parses the length of the specified barcode type from the pattern.
    Example: for ^(UMI:N{12}) returns 12
    """
    match = re.search(rf'{barcode_type}:?(?:N?{{(\d+)}}|([^)]+))', pattern)
    if match:
        barcode_body = match.group(1) or match.group(2)
        return int(barcode_body) if barcode_body.isdigit() else len(barcode_body)
    return 0


def add_nucleotide_cost(pattern: str, max_error=2) -> str:
    """Adds mismatch cost for fuzzy matched patterns (adapters)"""
    new_pattern = pattern
    for adapter in set(re.findall(ADAPTER_PATTERN_REGEX, pattern)):
        adapter_max_error = math.ceil(len(adapter) * max_error / 10)
        new_pattern = re.sub(rf'(?<![\w\[(])({adapter})(?![\w\])])',
                             rf'(\1){{s<={adapter_max_error}}}', new_pattern)
    return new_pattern


def replace_barcode_type_to_regex_group(pattern: str, barcode_type: str) -> str:
    """Replaces barcode type to a named regex group in the given pattern.
    Example: ^(UMI:N{12}) -> ^(?P<UMI>N{12})
    """
    return re.sub(rf'{barcode_type}(:)?', rf'?P<{barcode_type}>', pattern)


def add_brackets_around_barcode(pattern: str, barcode_type: str) -> str:
    """Adds brackets around the specified barcode type in the pattern if not already present.
    Example: ^UMI:N{12} -> ^(UMI:N{12})
    """
    if not re.search(rf'\({barcode_type}:[A-Z0-9{{}}]+\)', pattern):
        pattern = re.sub(rf'({barcode_type}:[A-Z0-9{{}}]+)', r'(\1)', pattern)
    return pattern


def validate_pattern(pattern: str, umi_len: int):
    """Validates the pattern to ensure it contains a UMI placeholder and has a valid length."""
    if 'UMI' not in pattern:
        raise ValidationError(f'UMI placeholder does not found in {pattern}, exiting...')
    if not umi_len:
        raise ValidationError(f'UMI length in the pattern {pattern} should be > 0, exiting...')
    # TODO! invalid patterns: '^UMI{12}' -> ^(UMI:N{12}); '^UMI:N{12}' -> '^(UMI:N{12})'


def get_prepared_pattern_and_umi_len(pattern: str, max_error=2) -> tuple[str, int]:
    """Prepares a pattern for barcode matching by applying various transformations."""

    if not pattern:
        return '', 0

    logger.info(f"Converting '{pattern}' into regex format...")

    umi_len = parse_umi_length(pattern)
    validate_pattern(pattern, umi_len)
    pattern = pattern.upper()
    pattern = add_brackets_around_barcode(pattern, 'UMI')
    pattern = replace_barcode_type_to_regex_group(pattern, 'UMI')
    pattern = pattern.replace('N', f'[{ALLOWED_LETTERS_IN_UMI}]')
    pattern = add_nucleotide_cost(pattern, max_error)
    pattern = pattern.replace('{*}', '*')

    logger.info(f"Pattern has been converted into '{pattern}'...")

    return pattern, umi_len
