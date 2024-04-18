import os
import sys
import tempfile
import pyfastx
import regex
import re
import multiprocessing
from collections import namedtuple

from logger import set_logger

logger = set_logger(name=__file__)

BarcodeMarkup = namedtuple('BarcodeMarkup', 'sequence quality')
PatternMarkup = namedtuple('PatternMarkup', 'barcode start end')

TRANSLATION_TABLE = bytes.maketrans(b"ATGCRYSWKMBDHVN", b"TACGAAAAAAAAAAA")

ALLOWED_LETTERS_IN_UMI = 'ATGCN'

NORMAL_NUCLEOTIDES = ['a', 't', 'g', 'c']
IUPAC_WILDCARDS = [
    'r',  # A or G
    'y',  # C or T
    's',  # G or C
    'w',  # A or T
    'k',  # G or T
    'm',  # A or C
    'b',  # C or G or T
    'd',  # A or G or T
    'h',  # A or C or T
    'v',  # A or C or G
    'n'  # any base
]


class BarcodePattern:
    def __init__(self, pattern: str, max_error=10, wildcard_cost=1, normal_cost=2):
        self.max_error = max_error
        self.pattern = pattern
        self.wildcard_cost = wildcard_cost
        self.normal_cost = normal_cost
        self.umi_len = self._parse_barcode_length()

    def get_prepared_pattern(self) -> str:

        if not self.pattern:
            return ''

        logger.info(f'Preparing pattern {self.pattern}...')

        self._validate_pattern(self.pattern, self.umi_len)
        pattern = self._add_brackets_around_barcode(self.pattern, barcode_type='UMI')
        pattern = self._replace_barcode_type_to_regex_group(pattern, barcode_type='UMI')
        pattern = self._replace_nucleotide_patterns(pattern)
        pattern = self._add_nucleotide_cost(pattern, ''.join(IUPAC_WILDCARDS), mismatch_cost=self.wildcard_cost)
        pattern = self._add_nucleotide_cost(pattern, ''.join(NORMAL_NUCLEOTIDES), mismatch_cost=self.normal_cost)
        pattern = pattern.replace('{*}', '*')
        pattern = self._replace_lowercase_to_uppercase(pattern)

        logger.info(f'Pattern prepared and converted into {pattern}')
        return pattern

    def _add_nucleotide_cost(self, pattern: str, nucleotides: str, mismatch_cost: int):
        """
        Adds mismatch cost for fuzzy subpattern (lower-case letters)

        Example:
        ^N{0:2}tggtatcaacgcagagt(UMI:N{14}) -> ^N{0:2}(tggtatcaacgcagagt){2s<=10}(UMI:N{14}),
        where 2s<=10: s - single nucleotide substitution, 2 - cost of one substitution, 10 - max mismatches count
        """
        return re.sub(rf'([{nucleotides}]+)', rf'(\1){{{mismatch_cost}s<={self.max_error}}}', pattern)

    def _replace_lowercase_to_uppercase(self, pattern):
        """
        Replaces lowercase nucleotide letters to uppercase

        ^N{0:2}(tggtatcaacgcagagt){2s<=10}(UMI:N{14}) -> ^N{0:2}(TGGTATCAACGCAGAGT){2s<=10}(UMI:N{14})
        """
        all_nucleotide_letters = ''.join(NORMAL_NUCLEOTIDES) + ''.join(IUPAC_WILDCARDS)
        return re.sub(rf'([{all_nucleotide_letters}]+)', lambda m: m.group(1).upper(), pattern)

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


class BarcodeExtractor:
    def __init__(self, in_fq1_path: str, in_fq2_path: str, read1_pattern: str, read2_pattern: str,
                 out_fq1_path=None, out_fq2_path=None, add_to_the_header=False, keep_reads_without_adapter=False,
                 find_umi_in_rc=False):

        self.read1_pattern = read1_pattern
        self.read2_pattern = read2_pattern
        self.find_umi_in_rc = find_umi_in_rc
        self.add_to_the_header = add_to_the_header
        self.keep_reads_without_adapter = keep_reads_without_adapter
        self.initial_reads_count, self.final_reads_count = 0, 0
        self.in_fq1_path, self.in_fq2_path = in_fq1_path, in_fq2_path
        self.out_fq1_path = out_fq1_path or tempfile.NamedTemporaryFile().name
        self.out_fq2_path = out_fq2_path or tempfile.NamedTemporaryFile().name

    def _process_umi_in_read(self, read1: list[str], read2: list[str], pattern: str, is_read2=False):
        pattern_markup, read1, read2 = self._match_barcode(read1, read2, pattern, is_read2)
        if pattern_markup.barcode.sequence or self.keep_reads_without_adapter:
            read_header = self._add_barcode_to_the_header(read1[0], pattern_markup.barcode.sequence) \
                if self.add_to_the_header else read1[0]
            read_seq, read_quality = self._replace_umi_to_the_seq_start(read1[1], read1[2], pattern_markup)
            return self._create_new_read(read_header, read_seq, read_quality)
        return ''

    def _create_new_read(self, header: str, seq: str, quality: str) -> str:
        return f"@{header}\n" \
               f"{seq}\n" \
               f"+\n" \
               f"{quality}\n"

    def process_barcode(self, read1: list[str], read2: list[str]) -> tuple[str, str]:
        new_read1, new_read2 = self._create_new_read(*read1), self._create_new_read(*read2)
        if self.read1_pattern:
            new_read1 = self._process_umi_in_read(read1, read2, self.read1_pattern)
        if self.read2_pattern:
            new_read2 = self._process_umi_in_read(read2, read1, self.read2_pattern, is_read2=True)
        return new_read1, new_read2

    def get_fastq_reads(self) -> tuple[list[str], list[str]]:
        logger.info(f'Reading fastq files...')

        reads1 = [seq for seq in pyfastx.Fastq(self.in_fq1_path, build_index=False, full_name=True)]
        reads2 = [seq for seq in pyfastx.Fastq(self.in_fq2_path, build_index=False, full_name=True)]

        self.initial_reads_count = len(reads1)
        logger.info(f'Initial read count in fastq: {self.initial_reads_count}')
        logger.info(f'Fastq files have been read.')

        return reads1, reads2

    def _get_reverse_complement(self, sequence: str) -> str:
        """
        Returns reverse complement sequence. Example: TTTAAAGGG -> CCCTTTAAA.
        """
        return sequence.translate(TRANSLATION_TABLE)[::-1]

    def _write(self, reads: list[str], output_path: str):
        with open(output_path, 'w') as f:
            for read in reads:
                f.write(read)

    def process_in_parallel(self, reads1: list[str], reads2: list[str]) -> tuple[str, str]:

        with multiprocessing.Pool(processes=os.cpu_count()) as pool:
            args_list = zip(reads1, reads2)
            processed_reads = pool.starmap(self.process_barcode, args_list)

        processed_reads1, processed_reads2 = [], []
        for read1, read2 in processed_reads:
            if read1 and read2:
                processed_reads1.append(read1)
                processed_reads2.append(read2)

        self._write(processed_reads1, self.out_fq1_path)
        self._write(processed_reads2, self.out_fq2_path)

        self.final_reads_count = len(processed_reads1)
        logger.info(f'Final read count in fastq: {self.final_reads_count}')

        return self.out_fq1_path, self.out_fq2_path

    def _add_barcode_to_the_header(self, header: str, umi_seq: str, barcode_type='UMI'):
        """
        Adds UMI to the header.

        Example:
        @M01691:10:000000000-A7F7L:1:1101:13747:1534 2:N:0:1 UMI:CCGTGAGCTGTTT
        """
        return ' '.join([header, f'{barcode_type}:{umi_seq}']) if umi_seq else header

    def _match_in_rc(self, read_seq, pattern):
        read_seq_rc = self._get_reverse_complement(read_seq)
        match = regex.search(pattern, read_seq_rc, regex.BESTMATCH)
        return match

    def _get_barcode_fields(self, read, match, barcode_type='UMI'):
        barcode_seq = match.group(barcode_type).replace('N', 'A')
        barcode_start, barcode_end = match.span(barcode_type)
        barcode_quality = read[2][barcode_start:barcode_end]
        pattern_match_start, pattern_match_end = match.span()
        return PatternMarkup(BarcodeMarkup(barcode_seq, barcode_quality), pattern_match_start, pattern_match_end)

    def _match_barcode(self, read1: list[str], read2: list[str],
                       pattern: str, is_read2: bool) -> tuple[PatternMarkup, list[str], list[str]]:
        match = regex.search(pattern, read1[1], regex.BESTMATCH)
        if not match:
            match = regex.search(pattern, read2[1], regex.BESTMATCH)
            if match:
                read1, read2 = read2, read1  # replace reads if pattern is matched in another read
        if not match and (self.find_umi_in_rc or not is_read2):
            read_seq_rc = self._get_reverse_complement(read1[1])
            match = regex.search(pattern, read_seq_rc, regex.BESTMATCH)
        if not match and self.find_umi_in_rc:
            read_seq_rc = self._get_reverse_complement(read2[1])
            match = regex.search(pattern, read_seq_rc, regex.BESTMATCH)
        if match:
            return self._get_barcode_fields(read1, match), read1, read2
        return PatternMarkup(BarcodeMarkup('', ''), 0, 0), read1, read2

    def _replace_umi_to_the_seq_start(self, read_seq: str, read_quality: str,
                                      pattern: PatternMarkup) -> tuple[str, str]:
        if not read_seq:
            return read_seq, read_quality
        new_read_seq = pattern.barcode.sequence + self._remove_subseq(read_seq, pattern.start, pattern.end)
        new_read_quality = pattern.barcode.quality + self._remove_subseq(read_quality, pattern.start, pattern.end)
        return new_read_seq, new_read_quality

    def _remove_subseq(self, sequence: str, subseq_start: int, subseq_end: int) -> str:
        if subseq_start == 0 and subseq_end == 0:
            return sequence
        return sequence[:subseq_start] + sequence[subseq_end:]

    def get_initial_reads_count(self):
        """
        Return count of all reads in fastq
        """
        return self.initial_reads_count

    def get_final_reads_count(self) -> int:
        """
        Returns count of reads in which patter is found.
        """
        return self.final_reads_count
