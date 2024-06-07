import gzip
import itertools
import tempfile

import pyfastx
import os
import multiprocessing

from logger import set_logger
logger = set_logger(name=__file__)

TRANSLATION_TABLE = bytes.maketrans(b"ATGCRYSWKMBDHVN", b"TACGAAAAAAAAAAA")


def get_fastq_reads(fq1_path: str, fq2_path: str) -> tuple[list[str], list[str]]:
    logger.info("Reading FASTQ files...")

    reads1 = [seq for seq in pyfastx.Fastq(fq1_path, build_index=False, full_name=True)]
    reads2 = [seq for seq in pyfastx.Fastq(fq2_path, build_index=False, full_name=True)]

    logger.info("FASTQ files have been read.")

    return reads1, reads2


def save_fastq_reads_to_file(reads: list[str], output_path: str):
    with gzip.open(output_path, "wb") as f:
        for read in reads:
            f.write(read.encode())


def get_reverse_complement(read_sequence: str) -> str:
    """Returns reverse complement of nucleotide sequence"""
    return read_sequence.translate(TRANSLATION_TABLE)[::-1]


def mock_merge_one_reads_pair(read1: list[str], read2: list[str], insert_size: int) -> str:
    """Perform a mock merging for non-overlapping reads"""
    new_header = f"{read1[0]} mock_merged_{len(read1[1])}_{len(read2[1])}"
    new_read_sequence = read1[1] + "N" * insert_size + get_reverse_complement(read2[1])
    new_read_quality = read1[2] + "#" * insert_size + read2[2]
    return f"@{new_header}\n" \
           f"{new_read_sequence}\n" \
           f"+\n" \
           f"{new_read_quality}\n"


def mock_merge_reads(fq1_path: str, fq2_path: str, insert_size) -> str:
    """Run mock merging in parallel"""
    reads1, reads2 = get_fastq_reads(fq1_path, fq2_path)

    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        args_list = zip(reads1, reads2, itertools.repeat(insert_size, len(reads1)))
        processed_reads = pool.starmap(mock_merge_one_reads_pair, args_list)

    merged_fq_path = tempfile.NamedTemporaryFile().name

    save_fastq_reads_to_file(processed_reads, merged_fq_path)

    return merged_fq_path
