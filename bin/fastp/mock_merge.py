import gzip

import os
import multiprocessing

from logger import set_logger
logger = set_logger(name=__file__)

TRANSLATION_TABLE = str.maketrans("ATGCRYSWKMBDHVN", "TACGAAAAAAAAAAA")
READS_CHUNK_SIZE = 1000000


def read_fastq_file_chunk(file) -> list[tuple[str, str, str]]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    reads_list = []
    while header := file.readline().strip():  # read header string
        if not header or len(reads_list) >= READS_CHUNK_SIZE:
            break
        sequence = file.readline().strip()  # read sequence string
        _ = file.readline().strip()  # read separator string
        quality = file.readline().strip()  # read quality string
        reads_list.append((header, sequence, quality))
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


def append_non_empty_file(src_path: str, dest_path: str):
    """Appended lines of source file to destination file"""
    with gzip.open(src_path, "rb") as src_file:
        with gzip.open(dest_path, "ab") as dest_file:
            dest_file.write(src_file.read())


def run_mock_merge_by_read_chunks(fq1_path: str, fq2_path: str, inner_distance_size: int, out_fq12: str):
    """Run mock merging in parallel by chunks"""
    with (gzip.open(fq1_path, "rt") as fq1_file_obj,
          gzip.open(fq2_path, "rt") as fq2_file_obj,
          gzip.open(out_fq12, "ab") as out_fq12_file_obj):

        fq1_reads_list = read_fastq_file_chunk(fq1_file_obj)
        fq2_reads_list = read_fastq_file_chunk(fq2_file_obj)

        while fq1_reads_list and fq2_reads_list:
            with multiprocessing.Pool(processes=os.cpu_count()) as pool:
                processed_reads = pool.starmap(mock_merge_one_reads_pair,
                                               zip(fq1_reads_list, fq2_reads_list, [inner_distance_size] * len(fq1_reads_list)))

            out_fq12_file_obj.writelines((read.encode("utf-8") for read in processed_reads))

            fq1_reads_list = read_fastq_file_chunk(fq1_file_obj)
            fq2_reads_list = read_fastq_file_chunk(fq2_file_obj)


def run_mock_merge_reads(fq1_path: str, fq2_path: str, inner_distance_size: int, out_fq12: str):
    if os.path.getsize(fq1_path) == 0:
        append_non_empty_file(fq2_path, out_fq12)
    elif os.path.getsize(fq2_path) == 0:
        append_non_empty_file(fq1_path, out_fq12)
    else:
        run_mock_merge_by_read_chunks(fq1_path, fq2_path, inner_distance_size, out_fq12)
