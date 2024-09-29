import typing
from collections import namedtuple
import gzip
from itertools import repeat
import os
import multiprocessing

from logger import set_logger

logger = set_logger(name=__file__)

CPU_COUNT = os.cpu_count()

TRANSLATION_TABLE = str.maketrans("ATGCRYSWKMBDHVN", "TACGAAAAAAAAAAA")

FastqRead = namedtuple('FastqRead', ['header', 'sequence', 'quality'])


def read_fastq_file_chunk(file_obj: typing.TextIO, reads_chunk_size: int) -> list[FastqRead]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    reads_list = []
    while len(reads_list) < reads_chunk_size and (read_header := file_obj.readline().strip()):
        read_sequence = file_obj.readline().strip()
        file_obj.readline()
        read_quality = file_obj.readline().strip()
        reads_list.append(FastqRead(read_header, read_sequence, read_quality))
    return reads_list


def get_reverse_complement(read_sequence: str) -> str:
    """Returns reverse complement of nucleotide sequence"""
    return read_sequence.translate(TRANSLATION_TABLE)[::-1]


def mock_merge_one_reads_pair(read1: FastqRead, read2: FastqRead, inner_distance_size: int) -> str:
    """Perform a mock merging for non-overlapping reads"""
    new_header = f"{read1.header} mock_merged_{len(read1.sequence)}_{len(read2.sequence)}"
    new_read_sequence = read1.sequence + "N" * inner_distance_size + get_reverse_complement(read2.sequence)
    new_read_quality = read1.quality + "#" * inner_distance_size + read2.quality
    return f"{new_header}\n" \
           f"{new_read_sequence}\n" \
           f"+\n" \
           f"{new_read_quality}\n"


def process_reads_chunk_in_parallel(fq1_reads_chunk_list: list[tuple[str, str, str]],
                                    fq2_reads_chunk_list: list[tuple[str, str, str]],
                                    inner_distance_size: int) -> tuple[bytes]:
    """Merges read chunk in parallel"""
    with multiprocessing.Pool(processes=CPU_COUNT) as pool:
        processed_reads_chunk = pool.starmap(mock_merge_one_reads_pair,
                                             zip(fq1_reads_chunk_list, fq2_reads_chunk_list,
                                                 repeat(inner_distance_size)))
    return tuple(read.encode() for read in processed_reads_chunk)


def append_merged_reads_to_fastq_file(out_fq12_file_obj: gzip.GzipFile, mock_merged_read_chunk: tuple[bytes]) -> None:
    """Append mock merged reads into fq12 file"""
    out_fq12_file_obj.writelines(mock_merged_read_chunk)


def mock_merge_by_chunks(fq1_path: str, fq2_path: str, inner_distance_size: int,
                         reads_chunk_size: int, out_fq12: str) -> None:
    """Mock merge reads by chunk and append merged read to fq12 file"""
    logger.info("Going to perform mock merge reads for FASTQ1 and FASTQ2 files...")
    with (gzip.open(fq1_path, "rt") as fq1_file_obj,
          gzip.open(fq2_path, "rt") as fq2_file_obj,
          gzip.open(out_fq12, "ab") as out_fq12_file_obj):

        fq1_reads_chunk_list = read_fastq_file_chunk(fq1_file_obj, reads_chunk_size)
        fq2_reads_chunk_list = read_fastq_file_chunk(fq2_file_obj, reads_chunk_size)

        while fq1_reads_chunk_list and fq2_reads_chunk_list:
            mock_merged_read_chunk = process_reads_chunk_in_parallel(fq1_reads_chunk_list, fq2_reads_chunk_list,
                                                                     inner_distance_size)

            append_merged_reads_to_fastq_file(out_fq12_file_obj, mock_merged_read_chunk)

            fq1_reads_chunk_list = read_fastq_file_chunk(fq1_file_obj, reads_chunk_size)
            fq2_reads_chunk_list = read_fastq_file_chunk(fq2_file_obj, reads_chunk_size)
            logger.info(f"Successfully processed chunk with '{len(mock_merged_read_chunk)}' reads.")

    logger.info("All reads successfully merged.")
