import json
import os
import sys
from typing import Optional, Union

import shutil
import subprocess
import tempfile

from barcode_pattern import BarcodePattern
from barcode_extractor import BarcodeExtractor
from logger import set_logger

logger = set_logger(name=__file__)

TMP_DIR = '/tmp'
FASTQ_CHUNK_SIZE = 2_000_000  # reads count in one fastq chunk


def save_total_read_count(total_reads_count: int, out_json_path: str):
    json_content = json.dumps({
        "summary": {
            "before_filtering": {
                "total_reads": int(total_reads_count)
            }
        }
    })
    save_to_file(json_content, out_json_path)
    logger.info(f'Total read count has been saved into {out_json_path}.')


def exit_with_error(message: Optional[str]):
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


def run_command(command: list[str], stdin=None, stdout=False) -> Union[str, None]:
    logger.info(f'Running {command}...')

    try:
        command_process = subprocess.run(command, text=True, capture_output=True, stdin=stdin, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        exit_with_error(f"Failed to run {command}")
    except Exception as e:
        exit_with_error(f"Undefined error: {e}")

    if stdout:
        return command_process.stdout


def concat_files(files: list[str]) -> str:
    output_file = tempfile.NamedTemporaryFile().name
    with open(output_file, 'w') as out_f:
        for file in files:
            with open(file, 'r') as f:
                out_f.write(f.read())
    remove(*files)
    return output_file


def remove(*file: str):
    for f in file:
        os.remove(f)


def save_to_file(data: str, file_path=None) -> str:
    file_path = file_path or tempfile.NamedTemporaryFile().name
    with open(file_path, 'w') as f:
        f.write(data)
    return file_path


def compress(file: str):
    logger.info(f'Compressing {file} into {file}.gz...')

    pigz_cmd = ['pigz', file]
    run_command(pigz_cmd)

    file_gz = file + '.gz'
    check_if_exist(file_gz)

    logger.info(f'{file} has been compressed.')

    return file_gz


def replace_file(src_file: str, dst_file: str):
    """
    Replaces file from source to the destination

    :param src_file: path to the file, that we need to move
    :param dst_file: path to destination file
    """
    shutil.move(src_file, dst_file)
    logger.info(f'{src_file} replaced to {dst_file}')


def check_if_exist(file: str):
    """
    Checks if file exists

    :param file: path to the file
    """
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    if os.path.getsize(file) == 0:
        logger.critical(f'Output file {file} is an empty, exiting..')
        sys.exit(1)
    logger.info(f"Expected file {file} found.")


def save_results(read1_file: str, read2_file: str,
                 read1_out_file: str, read2_out_file: str):
    for file, out_file in [(read1_file, read1_out_file),
                           (read2_file, read2_out_file)]:
        check_if_exist(file)
        file_gz = compress(file)
        replace_file(file_gz, out_file)


def generate_consensus(in_fq1: str, in_fq2: str, cluster_file: str,
                       min_reads_per_cluster: int, max_reads_per_cluster: int) -> tuple[str, str]:
    logger.info('Generating consensus of clustered reads...')

    fq1_cons_prefix, fq2_cons_prefix = tempfile.NamedTemporaryFile().name, tempfile.NamedTemporaryFile().name
    calib_cons_cmd = ['calib_cons',
                      '--fastq', in_fq1, in_fq2,
                      '--cluster', cluster_file,
                      '--output-prefix', fq1_cons_prefix, fq2_cons_prefix,
                      '--threads', str(os.cpu_count()),
                      '--min-reads-per-cluster', str(min_reads_per_cluster),
                      '--max-reads-per-cluster', str(max_reads_per_cluster)]

    run_command(calib_cons_cmd)

    logger.info('Consensus generation has been done.')

    fq1_cons = fq1_cons_prefix + '.fastq'
    fq2_cons = fq2_cons_prefix + '.fastq'

    return fq1_cons, fq2_cons


def cluster_umi(in_fq1: str, in_fq2: str, fq1_umi_len: int, fq2_umi_len: int,
                kmer_size: int, minimizer_count: int, minimizer_threshold: int, error_tolerance: int) -> str:
    output_prefix = tempfile.NamedTemporaryFile().name
    logger.info(f'Clustering umi in {in_fq1} and {in_fq2} into {output_prefix}...')

    calib_cmd = ['calib',
                 '--input-forward', in_fq1,
                 '--input-reverse', in_fq2,
                 '--output-prefix', output_prefix,
                 '--threads', str(os.cpu_count()),
                 '--barcode-length-1', str(fq1_umi_len),
                 '--barcode-length-2', str(fq2_umi_len),
                 '--kmer-size', str(kmer_size),
                 '--minimizer-count', str(minimizer_count),
                 '--minimizer-threshold', str(minimizer_threshold),
                 '--error-tolerance', str(error_tolerance)]

    run_command(calib_cmd)

    cluster_file = output_prefix + 'cluster'
    check_if_exist(cluster_file)

    logger.info('UMI clustering has been done.')

    return cluster_file


def extract_umi(fq12_chunks: list[str], read1_pattern: str,
                read2_pattern: str, find_umi_in_rc: bool) -> tuple[str, str, int]:
    total_reads_count, initial_reads_count = 0, 0
    processed_fq1_chunks, processed_fq2_chunks = [], []

    logger.info(f'Extracting UMI...')
    for fq1_chunk, fq2_chunk in fq12_chunks:
        umi_extractor = BarcodeExtractor(fq1_chunk, fq2_chunk, read1_pattern,
                                         read2_pattern, find_umi_in_rc=find_umi_in_rc)

        chunk_reads1, chunk_reads2 = umi_extractor.get_fastq_reads()
        processed_fq1_chunk, processed_fq2_chunk = umi_extractor.process_in_parallel(chunk_reads1, chunk_reads2)

        processed_fq1_chunks.append(processed_fq1_chunk)
        processed_fq2_chunks.append(processed_fq2_chunk)

        total_reads_count += umi_extractor.get_initial_reads_count()
        initial_reads_count += umi_extractor.get_final_reads_count()

        remove(fq1_chunk, fq2_chunk)

    logger.info('UMI successfully extracted.')

    processed_fq1 = concat_files(processed_fq1_chunks)
    processed_fq2 = concat_files(processed_fq2_chunks)

    return processed_fq1, processed_fq2, total_reads_count


def split_by_chunks(fq1_path: str, fq2_path: str) -> list[str]:
    logger.info(f'Splitting {fq1_path} and {fq2_path} into chunks by {FASTQ_CHUNK_SIZE} reads...')

    fq1_outdir = os.path.join(TMP_DIR, 'chunks', 'fq1')
    fq2_outdir = os.path.join(TMP_DIR, 'chunks', 'fq2')

    for fq_path, output_dir in [(fq1_path, fq1_outdir), (fq2_path, fq2_outdir)]:
        cmd = ['seqkit', 'split2', fq_path, '--by-size', str(FASTQ_CHUNK_SIZE),
               '--by-size-prefix', '', '--out-dir', output_dir]
        run_command(cmd)

    fq1_chunks = sorted([entry.path for entry in os.scandir(fq1_outdir)])
    fq2_chunks = sorted([entry.path for entry in os.scandir(fq2_outdir)])

    logger.info(f'Splitting {fq1_path} and {fq2_path} has been done.')

    return zip(fq1_chunks, fq2_chunks)


def keep_only_paired_reads(fq1: str, fq2: str, clear=False):
    logger.info('Filter out unpaired reads...')

    outdir = os.path.join(TMP_DIR, 'paired')
    cmd = ['seqkit', 'pair', '-1', fq1, '-2', fq2, '--out-dir', outdir]

    run_command(cmd)

    new_fq1 = os.path.join(outdir, os.path.basename(fq1))
    new_fq2 = os.path.join(outdir, os.path.basename(fq2))

    if clear:
        remove(fq1, fq2)

    logger.info('Filtering unpaired reads has been done.')

    return new_fq1, new_fq2


def check_error_tolerance_size(umi_len: int, error_tolerance: int):
    if umi_len and error_tolerance > umi_len:
        logger.critical(
            f'Error tolerance size (={error_tolerance}) should be no larger than the umi length (={umi_len}), exiting...'
        )
        sys.exit(1)


def get_prepared_pattern_and_umi_len(pattern: str, max_error: int, error_tolerance: int) -> tuple[str, int]:
    barcode_pattern = BarcodePattern(pattern, max_error=max_error)
    umi_length = len(barcode_pattern)
    check_error_tolerance_size(umi_length, error_tolerance)
    prepared_pattern = str(barcode_pattern)
    return prepared_pattern, umi_length