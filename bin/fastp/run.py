import argparse
import shutil
import subprocess
import os
import sys
import tempfile

from mock_merge import run_mock_merge_reads

from logger import set_logger

logger = set_logger(name=__file__)


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input fastq.gz file, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2', help='Input fastq.gz file, PE pair 2')
    parser.add_argument('--trimq',
                        help='Reads trimming: minimal base quality for 5'' & 3'' read ends to keep them', type=int)
    parser.add_argument('--disable-filters', help='Mode(s) to disable',
                        choices=["length_filtering", "quality_filtering", "adapter_trimming"], nargs='+',
                        action='extend')
    parser.add_argument('--merge-reads', help='Enable merge reads of overlapped reads', action='store_true')
    parser.add_argument('--mock-merge-reads', help='Enable mock merging of non-overlapped reads', action='store_true')
    parser.add_argument('--inner-distance-size', help='Insert size for mock merging', type=int)
    parser.add_argument('--out-fq1', help='Output fastq file, SE or PE pair 1')
    parser.add_argument('--out-fq2', help='Output fastq file, PE pair 2')
    parser.add_argument('--out-fq12', help='Output merged fastq file, PE pairs 1 and 2')
    parser.add_argument('--html', help='Output html file', required=True)
    parser.add_argument('--json', help='Output json file', required=True)

    return vars(parser.parse_args())


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


def run_command(command: list[str], stdin=None):
    logger.info(f'Running {command}...')

    try:
        command_process = subprocess.run(command, text=True, capture_output=True, input=stdin, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        logger.critical(f"Failed to run '{' '.join(command)}'")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"Undefined error: {e}")
        sys.exit(1)

    logger.info(command_process.stdout)
    logger.info(command_process.stderr)


def run_fastp(in_fq1: str, in_fq2: str, trimq: int, disable_filters: list[str], merge: bool, json: str, html: str,
              out_fq1: str, out_fq2: str, out_fq12: str):
    cmd = ['fastp', '-i', in_fq1]
    cmd += ['-o', out_fq1]

    if in_fq2:
        cmd += ['-I', in_fq2, '-O', out_fq2, '--detect_adapter_for_pe']

    if merge:
        cmd += ['--merge', '--merged_out', out_fq12]

    cmd += ['--cut_tail', '--cut_front', '--cut_window_size', str(1), '--cut_mean_quality',
            str(trimq)] if trimq else []
    cmd += [f'--disable_{mode}' for mode in disable_filters] if disable_filters else []
    cmd += ['--thread', str(os.cpu_count()), '--html', html, '--json', json]

    run_command(cmd)


def save_final_fastq_by_mode(merge_reads: bool, mock_merge_reads: bool, in_fq2_path: str, processed_fq1_path: str,
                             processed_fq2_path: str, processed_fq12_path: str, out_fq1_path: str, out_fq2_path: str,
                             out_fq12_path: str):
    """Saves final FASTQ files by provided modes"""
    if not mock_merge_reads:
        replace_file(processed_fq1_path, out_fq1_path)

    if in_fq2_path and not mock_merge_reads:
        replace_file(processed_fq2_path, out_fq2_path)

    if merge_reads or mock_merge_reads:
        replace_file(processed_fq12_path, out_fq12_path)


def run(in_fq1: str, in_fq2: str, trimq: int, disable_filters: list[str], merge_reads: bool, mock_merge_reads: bool,
        inner_distance_size: int, out_fq1: str, out_fq2: str, out_fq12: str, html: str, json: str) -> None:
    merge_reads = True if mock_merge_reads else merge_reads

    processed_fq1_path = tempfile.NamedTemporaryFile(suffix=".gz").name
    processed_fq2_path = tempfile.NamedTemporaryFile(suffix=".gz").name
    processed_fq12_path = tempfile.NamedTemporaryFile(suffix=".gz").name

    run_fastp(in_fq1, in_fq2, trimq, disable_filters, merge_reads, json, html,
              processed_fq1_path, processed_fq2_path, processed_fq12_path)

    if mock_merge_reads:
        run_mock_merge_reads(processed_fq1_path, processed_fq2_path, inner_distance_size, processed_fq12_path)

    save_final_fastq_by_mode(merge_reads, mock_merge_reads, in_fq2, processed_fq1_path, processed_fq2_path,
                             processed_fq12_path,
                             out_fq1, out_fq2, out_fq12)


if __name__ == "__main__":
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {args}")

    run(**args)
    logger.info("Run is completed successfully.")
