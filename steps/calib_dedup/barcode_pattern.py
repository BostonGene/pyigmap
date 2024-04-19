import sys
import re

from logger import set_logger

logger = set_logger(name=__file__)

ALLOWED_LETTERS_IN_UMI = 'ATGCN'

NORMAL_NUCLEOTIDES = ['A', 'T', 'G', 'C']
IUPAC_WILDCARDS = [
    'R',  # A or G
    'Y',  # C or T
    'S',  # G or C
    'W',  # A or T
    'K',  # G or T
    'M',  # A or C
    'B',  # C or G or T
    'D',  # A or G or T
    'H',  # A or C or T
    'V',  # A or C or G
    'N'   # any base
]


class BarcodePattern:
    def __init__(self, pattern: str, max_error=10, wildcard_cost=4, normal_cost=6):
        self.max_error = max_error
        self.pattern = pattern
        self.wildcard_cost = wildcard_cost
        self.normal_cost = normal_cost
        self.umi_len = self._parse_barcode_length()

    def get_prepared_pattern(self) -> str:

        if not self.pattern:
            return ''

        logger.info(f'Preparing pattern {self.pattern}...')

        pattern = self.pattern.upper()
        self._validate_pattern(pattern, self.umi_len)
        pattern = self._add_brackets_around_barcode(pattern, barcode_type='UMI')
        pattern = self._replace_barcode_type_to_regex_group(pattern, barcode_type='UMI')
        pattern = self._replace_nucleotide_patterns(pattern)
        pattern = self._add_nucleotide_cost(pattern, IUPAC_WILDCARDS, mismatch_cost=self.wildcard_cost)
        pattern = self._add_nucleotide_cost(pattern, NORMAL_NUCLEOTIDES, mismatch_cost=self.normal_cost)
        pattern = pattern.replace('{*}', '*')

        logger.info(f'Pattern prepared and converted into {pattern}')
        return pattern

    def _add_nucleotide_cost(self, pattern: str, nucleotides: list[str], mismatch_cost: int):
        """
        Adds mismatch cost for fuzzy subpattern (lower-case letters)

        Example:
        ^N{0:2}tggtatcaacgcagagt(UMI:N{14}) -> ^N{0:2}(tggtatcaacgcagagt){2s<=10}(UMI:N{14}),
        where 2s<=10: s - single nucleotide substitution, 2 - cost of one substitution, 10 - max mismatches count
        """
        nucleotides = ''.join(nucleotides)
        return re.sub(rf'(?<![\w([])([{nucleotides}]+)(?![\w\])])', rf'(\1){{{mismatch_cost}s<={self.max_error}}}', pattern)

    def _parse_barcode_length(self, barcode_type='UMI') -> int:
        """
        Returns barcode length from pattern. Example: for the pattern ^(UMI:N{12}) returns 12.
        """
        if not self.pattern:
            return 0
        match = re.search(rf'{barcode_type}:N{{(\d+)}}|{barcode_type}:([^)]+)', self.pattern)
        if not match:
            match = re.search(r'{:?([^}]*)}', self.pattern)
        barcode_body = match.group(1)
        if barcode_body.isnumeric():
            return int(barcode_body)
        return len(barcode_body)

    def get_umi_length(self) -> int:
        return self.umi_len

    def _replace_barcode_type_to_regex_group(self, pattern: str, barcode_type: str) -> str:
        """
        Replaces ^(UMI:N{12}) -> ^(?P<UMI>N{12}).
        """
        return re.sub(rf'{barcode_type}(:)?', rf'?P<{barcode_type}>', pattern)

    def _replace_nucleotide_patterns(self, pattern: str) -> str:
        """
        Replaces N{} (or n{}) -> [ATGCN]{}.
        """
        return re.sub(r'[Nn](?={)', f'[{ALLOWED_LETTERS_IN_UMI}]', pattern)

    def _add_brackets_around_barcode(self, pattern: str, barcode_type: str):
        """
        Replaces ^UMI:N{12} -> ^(UMI:N{12})
        """
        if re.search(rf'\({barcode_type}:[A-Z0-9{{}}]+\)', pattern):
            return pattern
        return re.sub(rf'({barcode_type}:[A-Z0-9{{}}]+)', r'(\1)', pattern)

    def _validate_pattern(self, pattern: str, umi_len: int):
        if pattern and 'UMI' not in pattern:
            logger.critical(f'UMI placeholder does not found in {pattern}, exiting...')
            sys.exit(1)
        if not umi_len and pattern:
            logger.critical(f'UMI length in the pattern {pattern} should be > 0, exiting...')
            sys.exit(1)
