import argparse
import shutil
import subprocess
import os
import sys
import tempfile
from typing import Optional

from mock_merge import mock_merge_by_chunks

from logger import set_logger

logger = set_logger(name=__file__)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if args.merge_reads and args.mock_merge_reads:
        msg_list += ["--merge-reads and --mock-merge-reads cannot be used at the same time"]
    if args.inner_distance_size and not args.mock_merge_reads:
        msg_list += ["--inner-distance-size cannot be provided without --mock-merge-reads"]
    if args.mock_merge_reads and (args.out_fq1 or args.out_fq2):
        msg_list += ["--out-fq1 or --out-fq2 cannot be used with --mock-merge-reads"]
    if args.reads_chunk_size and not args.mock_merge_reads:
        msg_list += ["--insert-size cannot be provided without --mock-merge-reads"]
    return msg_list


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input fastq.gz file, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2', help='Input fastq.gz file, PE pair 2')
    parser.add_argument('--disable-filters', help='Mode(s) to disable',
                        choices=["length_filtering", "quality_filtering", "adapter_trimming"], nargs='+',
                        action='extend')
    parser.add_argument('--merge-reads', help='Enable merge reads of overlapped reads', action='store_true')
    parser.add_argument('--mock-merge-reads', help='Enable mock merging of non-overlapped reads', action='store_true')
    parser.add_argument('--reads-chunk-size', help='Read chunk size used in mock merging', type=int)
    parser.add_argument('--inner-distance-size', help='Insert size for mock merging', type=int)
    parser.add_argument('--out-fq1', help='Output fastq file, SE or PE pair 1')
    parser.add_argument('--out-fq2', help='Output fastq file, PE pair 2')
    parser.add_argument('--out-fq12', help='Output merged fastq file, PE pairs 1 and 2')
    parser.add_argument('--html', help='Output html file', required=True)
    parser.add_argument('--json', help='Output json file', required=True)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error("\n".join(error_message_list))

    return vars(args)


def move_file(src_file: str, dst_file: str):
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


def run_fastp(in_fq1: str, in_fq2: str, disable_filters: list[str], merge: bool, json: str, html: str) -> tuple[str, str, str]:

    out_fq1, out_fq2, out_fq12 = (tempfile.NamedTemporaryFile(suffix=".gz").name for _ in range(3))

    cmd = ['fastp', '-i', in_fq1]
    cmd += ['-o', out_fq1]

    if in_fq2:
        cmd += ['-I', in_fq2, '-O', out_fq2, '--detect_adapter_for_pe']

    if merge:
        cmd += ['--merge', '--merged_out', out_fq12]

    cmd += [f'--disable_{mode}' for mode in disable_filters] if disable_filters else []
    cmd += ['--thread', str(os.cpu_count()), '--html', html, '--json', json]

    run_command(cmd)

    return out_fq1, out_fq2, out_fq12


def save_final_fastq_by_mode(merge_reads: bool, mock_merge_reads: bool, in_fq2_path: str,
                             processed_fq1_path: str, processed_fq2_path: str, processed_fq12_path: str,
                             out_fq1_path: str, out_fq2_path: str, out_fq12_path: str):
    """Saves final FASTQ files by provided modes"""
    if mock_merge_reads:
        move_file(processed_fq12_path, out_fq12_path)
        return

    move_file(processed_fq1_path, out_fq1_path)

    if in_fq2_path:
        move_file(processed_fq2_path, out_fq2_path)

    if merge_reads:
        move_file(processed_fq12_path, out_fq12_path)


def run(in_fq1: str, in_fq2: Optional[str], disable_filters: Optional[list[str]], merge_reads: Optional[bool],
        mock_merge_reads: Optional[bool], reads_chunk_size: Optional[int], inner_distance_size: Optional[int],
        out_fq1: Optional[str], out_fq2: Optional[str], out_fq12: Optional[str], html: str, json: str) -> None:

    merge_reads = True if mock_merge_reads else merge_reads

    tmp_fq1, tmp_fq2, tmp_fq12 = run_fastp(in_fq1, in_fq2, disable_filters, merge_reads, json, html)

    if mock_merge_reads:
        mock_merge_by_chunks(tmp_fq1, tmp_fq2, inner_distance_size, reads_chunk_size, tmp_fq12)

    save_final_fastq_by_mode(merge_reads, mock_merge_reads, in_fq2, tmp_fq1, tmp_fq2,
                             tmp_fq12, out_fq1, out_fq2, out_fq12)


if __name__ == "__main__":
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {args}")

    run(**args)
    logger.info("Run is completed successfully.")
