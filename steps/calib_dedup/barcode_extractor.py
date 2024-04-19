import os
import tempfile
import pyfastx
import regex
import multiprocessing
from collections import namedtuple

from logger import set_logger

logger = set_logger(name=__file__)

BarcodeMarkup = namedtuple('BarcodeMarkup', 'sequence quality')
PatternMarkup = namedtuple('PatternMarkup', 'barcode start end')

TRANSLATION_TABLE = bytes.maketrans(b"ATGCRYSWKMBDHVN", b"TACGAAAAAAAAAAA")



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
