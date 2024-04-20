import sys
import re

from logger import set_logger

logger = set_logger(name=__file__)

NORMAL_NUCLEOTIDES = 'ATGC'
IUPAC_WILDCARDS = 'RYSWKMBDHVN'
ALLOWED_LETTERS_IN_UMI = NORMAL_NUCLEOTIDES + 'N'


class ValidationError(Exception):
    pass

def parse_umi_length(pattern: str, barcode_type='UMI') -> int:
    """Parses the length of the specified barcode type from the pattern.
    Example: ^(UMI:N{12}) -> 12
    """
    match = re.search(rf'{barcode_type}:(?:N{{(\d+)}}|([^)]+))', pattern)
    if match:
        barcode_body = match.group(1) or match.group(2)
        return int(barcode_body) if barcode_body.isdigit() else barcode_body.count('N')
    return 0

def add_nucleotide_cost(pattern: str, wildcard_cost=1, normal_cost=2, max_error=2) -> str:
    """Adds costs to nucleotide patterns in the given pattern to allow for a certain number of mismatches.

    Example:
    ^[ATGCN]{0:2}TGGTATCAACGCAGAGT(?P<UMI>[ATGCN]{14}) -> ^N{0:2}(TGGTATCAACGCAGAGT){2s<=2}(?P<UMI>[ATGCN]{14}),
    where 2s<=2: s - single nucleotide substitution, 2 - cost of one substitution, 2 - max mismatches count
    """
    pattern = re.sub(rf'(?<![\w([])([{IUPAC_WILDCARDS}]+)(?![\w\])])',
                     rf'(\1){{{wildcard_cost}s<={max_error}}}', pattern)
    pattern = re.sub(rf'(?<![\w([])([{NORMAL_NUCLEOTIDES}]+)(?![\w\])])',
                     rf'(\1){{{normal_cost}s<={max_error}}}', pattern)
    return pattern

def replace_nucleotide_patterns(pattern: str) -> str:
    """Replaces 'N' nucleotide patterns in the given pattern with a character class of allowed letters in the UMI.
    Example: ^(UMI:N{12}) -> ^(UMI:[ATGCN]{12})
    """
    return pattern.replace('N', f'[{ALLOWED_LETTERS_IN_UMI}]')

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


def get_prepared_pattern_and_umi_len(pattern: str, max_error=2, wildcard_cost=1, normal_cost=2) -> tuple[str, int]:
    """Prepares a pattern for barcode matching by applying various transformations."""
    if not pattern:
        return '', 0

    umi_len = parse_umi_length(pattern)
    pattern = pattern.upper()
    validate_pattern(pattern, umi_len)
    pattern = add_brackets_around_barcode(pattern, 'UMI')
    pattern = replace_barcode_type_to_regex_group(pattern, 'UMI')
    pattern = replace_nucleotide_patterns(pattern)
    pattern = add_nucleotide_cost(pattern, wildcard_cost, normal_cost, max_error)
    pattern = pattern.replace('{*}', '*')

    return pattern, umi_len