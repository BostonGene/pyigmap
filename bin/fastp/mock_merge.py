import mgzip

import os
import multiprocessing

from logger import set_logger

logger = set_logger(name=__file__)

TRANSLATION_TABLE = str.maketrans("ATGCRYSWKMBDHVN", "TACGAAAAAAAAAAA")


def read_fastq_file_chunk(file_obj, reads_chunk_size: int) -> list[tuple[str, str, str]]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    reads_list = []
    while len(reads_list) < reads_chunk_size and (read_header := file_obj.readline().strip()):
        read_sequence = file_obj.readline().strip()
        file_obj.readline()
        read_quality = file_obj.readline().strip()
        reads_list.append((read_header, read_sequence, read_quality))
    return reads_list


def get_reverse_complement(read_sequence: str) -> str:
    """Returns reverse complement of nucleotide sequence"""
    return read_sequence.translate(TRANSLATION_TABLE)[::-1]


def mock_merge_one_reads_pair(read1: list[str], read2: list[str], inner_distance_size: int) -> str:
    """Perform a mock merging for non-overlapping reads"""
    new_header = f"{read1[0]} mock_merged_{len(read1[1])}_{len(read2[1])}"
    new_read_sequence = read1[1] + "N" * inner_distance_size + get_reverse_complement(read2[1])
    new_read_quality = read1[2] + "#" * inner_distance_size + read2[2]
    return f"{new_header}\n" \
           f"{new_read_sequence}\n" \
           f"+\n" \
           f"{new_read_quality}\n"


def process_reads_chunk_in_parallel(fq1_reads_chunk_list, fq2_reads_chunk_list, inner_distance_size):
    """Merges read chunk in parallel"""
    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        processed_reads_chunk = pool.starmap(mock_merge_one_reads_pair,
                                             zip(fq1_reads_chunk_list, fq2_reads_chunk_list,
                                                 [inner_distance_size] * len(fq1_reads_chunk_list)))

    return (read.encode() for read in processed_reads_chunk)


def mock_merge_by_chunks(fq1_path: str, fq2_path: str, inner_distance_size: int, reads_chunk_size: int, out_fq12: str):
    """Mock merge reads by chunk and append merged read to fq12 file"""
    with (mgzip.open(fq1_path, "rt", thread=os.cpu_count()) as fq1_file_obj,
          mgzip.open(fq2_path, "rt", thread=os.cpu_count()) as fq2_file_obj,
          mgzip.open(out_fq12, "ab", thread=os.cpu_count()) as out_fq12_file_obj):

        fq1_reads_chunk_list = read_fastq_file_chunk(fq1_file_obj, reads_chunk_size)
        fq2_reads_chunk_list = read_fastq_file_chunk(fq2_file_obj, reads_chunk_size)

        while fq1_reads_chunk_list and fq2_reads_chunk_list:
            mock_merged_read_chunk = process_reads_chunk_in_parallel(fq1_reads_chunk_list, fq2_reads_chunk_list,
                                                                     inner_distance_size)

            # append mock merged reads into fq12 file
            out_fq12_file_obj.writelines(mock_merged_read_chunk)

            fq1_reads_chunk_list = read_fastq_file_chunk(fq1_file_obj, reads_chunk_size)
            fq2_reads_chunk_list = read_fastq_file_chunk(fq2_file_obj, reads_chunk_size)
