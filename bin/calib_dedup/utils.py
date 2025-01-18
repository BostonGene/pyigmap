import os
import sys
from typing import Optional, Union
import shutil
import subprocess
import tempfile

from logger import set_logger

logger = set_logger(name=__file__)

TMP_DIR = tempfile.gettempdir()


def exit_with_error(message: Optional[str]):
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


def run_command(command: Union[list[str], str], stdin=None, shell=False):
    logger.info(f'Running {command}...')

    try:
        _ = subprocess.run(command, text=True, capture_output=True, shell=shell, stdin=stdin, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        exit_with_error(f"Failed to run {command}")
    except Exception as e:
        exit_with_error(f"Undefined error: {e}")


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


def check_error_tolerance_size(umi_len: int, error_tolerance: int):
    if umi_len and error_tolerance > umi_len:
        logger.critical(
            f'Error tolerance size (={error_tolerance}) should be no larger than the umi length (={umi_len}), exiting...'
        )
        sys.exit(1)


def decompress(gz_file: str) -> str:
    """Decompresses gzipped file"""
    output_file = tempfile.NamedTemporaryFile().name
    cmd = ' '.join(['pigz', '-dck', gz_file, '>', output_file])
    run_command(cmd, shell=True)
    return output_file
