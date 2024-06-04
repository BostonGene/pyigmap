import itertools
import tempfile

import pyfastx
import os
import multiprocessing

from logger import set_logger
logger = set_logger(name=__file__)


def get_fastq_reads(fq1_path: str, fq2_path: str) -> tuple[list[str], list[str]]:
    logger.info("Reading FASTQ files...")

    reads1 = [seq for seq in pyfastx.Fastq(fq1_path, build_index=False, full_name=True)]
    reads2 = [seq for seq in pyfastx.Fastq(fq2_path, build_index=False, full_name=True)]

    logger.info("FASTQ files have been read.")

    return reads1, reads2


def save_fastq_reads_to_file(reads: list[str], output_path: str):
    with open(output_path, 'w') as f:
        for read in reads:
            f.write(read)


def _create_new_read(header: str, seq: str, quality: str) -> str:
    return f"@{header}\n" \
           f"{seq}\n" \
           f"+\n" \
           f"{quality}\n"


def mock_merge_one_reads_pair(read1: list[str], read2: list[str], insert_size: int) -> str:
    new_header = f"{read1[0]} merged_{len(read1[1])}_{len(read2[1])}"
    new_read_sequence = read1[1] + "N" * insert_size + read2[1]
    new_read_quality = read1[2] + "#" * insert_size + read2[2]
    return f"@{new_header}\n" \
           f"{new_read_sequence}\n" \
           f"+\n" \
           f"{new_read_quality}\n"


def mock_merge_reads(fq1_path: str, fq2_path: str, insert_size) -> str:
    reads1, reads2 = get_fastq_reads(fq1_path, fq2_path)

    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        args_list = zip(reads1, reads2, itertools.repeat(insert_size, len(reads1)))
        processed_reads = pool.starmap(mock_merge_one_reads_pair, args_list)

    merged_fq_path = tempfile.NamedTemporaryFile().name

    save_fastq_reads_to_file(processed_reads, merged_fq_path)

    return merged_fq_path
