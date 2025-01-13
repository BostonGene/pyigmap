import os
import tempfile
import pyfastx
import regex
import multiprocessing
from collections import namedtuple
from itertools import cycle

from logger import set_logger

logger = set_logger(name=__file__)

BarcodeMarkup = namedtuple('BarcodeMarkup', ['sequence', 'quality'])
PatternMarkup = namedtuple('PatternMarkup', ['barcode', 'start', 'end'])

TRANSLATION_TABLE = bytes.maketrans(b"ATGCRYSWKMBDHVN", b"TACGAAAAAAAAAAA")


def get_reverse_complement(sequence: str) -> str:
    """
    Returns reverse complement sequence. Example: TTTAAAGGG -> CCCTTTAAA.
    """
    return sequence.translate(TRANSLATION_TABLE)[::-1]


def write_new_reads_to_file(reads: list[str]) -> str:
    with tempfile.NamedTemporaryFile(mode='w', dir='tmp', delete=False) as temp_file:
        temp_file.writelines(reads)
        return temp_file.name


def add_barcode_to_the_header(header: str, umi_seq: str, barcode_type='UMI'):
    """Adds UMI to the header.

    Example:
    @M01691:10:000000000-A7F7L:1:1101:13747:1534 2:N:0:1 UMI:CCGTGAGCTGTTT
    """
    return ' '.join([header, f'{barcode_type}:{umi_seq}']) if umi_seq else header


def get_barcode_fields(read, match, barcode_type='UMI'):
    barcode_seq = ''.join(match.captures(barcode_type)).replace('N', 'A')
    barcode_start_end_indices = match.spans(barcode_type)
    barcode_quality = ''
    for (barcode_start, barcode_end) in barcode_start_end_indices:
        barcode_quality += read[2][barcode_start:barcode_end]
    pattern_match_start, pattern_match_end = match.span()
    # print(BarcodeMarkup(barcode_seq, barcode_quality), pattern_match_start, pattern_match_end)
    return PatternMarkup(BarcodeMarkup(barcode_seq, barcode_quality), pattern_match_start, pattern_match_end)


def replace_umi_to_the_seq_start(read_seq: str, read_quality: str,
                                 pattern: PatternMarkup) -> tuple[str, str]:
    if not read_seq:
        return read_seq, read_quality
    new_read_seq = pattern.barcode.sequence + remove_subseq(read_seq, pattern.start, pattern.end)
    new_read_quality = pattern.barcode.quality + remove_subseq(read_quality, pattern.start, pattern.end)
    return new_read_seq, new_read_quality


def remove_subseq(sequence: str, subseq_start: int, subseq_end: int) -> str:
    if subseq_start == 0 and subseq_end == 0:
        return sequence
    return sequence[:subseq_start] + sequence[subseq_end:]


def match_in_reverse_complement(read_seq, pattern):
    read_seq_rc = get_reverse_complement(read_seq)
    match = regex.search(pattern, read_seq_rc, regex.BESTMATCH)
    return match


def match_barcode(read1: list[str], read2: list[str], pattern: str,
                  is_read2: bool, find_umi_in_rc=True) -> tuple[PatternMarkup, list[str], list[str]]:

    def find_match(read):
        # print(read[1], pattern)
        match = regex.search(pattern, read[1], regex.BESTMATCH)
        if not match and (find_umi_in_rc or not is_read2):
            match = match_in_reverse_complement(read[1], pattern)
        return match

    match1 = find_match(read1)
    match2 = find_match(read2)

    if not match1 and match2:
        read1, read2 = read2, read1

    return (get_barcode_fields(read1, match1) if match1 else
            PatternMarkup(BarcodeMarkup('', ''), 0, 0)), read1, read2


def process_umi_in_read(read1: list[str], read2: list[str], pattern: str, find_umi_in_rc=True, is_read2=False,
                        keep_reads_without_adapter=False, add_to_the_header=False):
    pattern_markup, read1, read2 = match_barcode(read1, read2, pattern, is_read2, find_umi_in_rc)
    if pattern_markup.barcode.sequence or keep_reads_without_adapter:
        read_header = add_barcode_to_the_header(read1[0], pattern_markup.barcode.sequence) \
            if add_to_the_header else read1[0]
        read_seq, read_quality = replace_umi_to_the_seq_start(read1[1], read1[2], pattern_markup)
        return create_new_read(read_header, read_seq, read_quality)
    return ''


def create_new_read(header: str, seq: str, quality: str) -> str:
    return f"@{header}\n" \
           f"{seq}\n" \
           "+\n" \
           f"{quality}\n"


def process_barcode(read1: list[str], read2: list[str], read1_pattern: str,
                    read2_pattern: str, find_umi_in_rc: bool) -> tuple[str, str]:
    new_read1, new_read2 = create_new_read(*read1), create_new_read(*read2)
    if read1_pattern:
        new_read1 = process_umi_in_read(read1, read2, read1_pattern, find_umi_in_rc=find_umi_in_rc)
    if read2_pattern:
        new_read2 = process_umi_in_read(read2, read1, read2_pattern, find_umi_in_rc=find_umi_in_rc, is_read2=True)
    return new_read1, new_read2


def get_processed_reads(reads1: list[str], reads2: list[str], read1_pattern: str,
                        read2_pattern: str, find_umi_in_rc: bool) -> tuple[list, list]:
    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        args_list = list(zip(reads1, reads2, cycle([read1_pattern]), cycle([read2_pattern]), cycle([find_umi_in_rc])))
        processed_reads = [read for read in pool.starmap(process_barcode, args_list)]
    return map(list, zip(*processed_reads))


def get_fastq_reads(in_fq1_path: str, in_fq2_path: str) -> tuple[list[str], list[str]]:
    logger.info(f'Reading FASTQ files...')

    reads1 = [seq for seq in pyfastx.Fastq(in_fq1_path, build_index=False, full_name=True)]
    reads2 = [seq for seq in pyfastx.Fastq(in_fq2_path, build_index=False, full_name=True)]

    logger.info(f'Initial read count in FASTQ: {len(reads1)}')
    logger.info(f'FASTQ files have been read.')

    return reads1, reads2


def get_processed_fastqs(raw_fq1_path: str, raw_fq2_path: str, read1_pattern: str, read2_pattern: str,
                         find_umi_in_rc: bool) -> tuple[str, str, int, int]:
    reads1, reads2 = get_fastq_reads(raw_fq1_path, raw_fq2_path)
    processed_reads1, processed_reads2 = get_processed_reads(reads1, reads2, read1_pattern, read2_pattern, find_umi_in_rc)

    out_fq1_path = write_new_reads_to_file(processed_reads1)
    out_fq2_path = write_new_reads_to_file(processed_reads2)

    logger.info(f'Final read count in FASTQ: {len(processed_reads1)}')

    return out_fq1_path, out_fq2_path, len(reads1), len(processed_reads1)
