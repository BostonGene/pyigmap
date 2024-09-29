import argparse
import shutil
import subprocess
import os
import sys
import tempfile
from typing import Optional, List

from mock_merge import mock_merge_by_chunks

from logger import set_logger

logger = set_logger(name=__file__)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if args.merge_reads and args.mock_merge_reads:
        msg_list += ["--merge-reads and --mock-merge-reads cannot be used at the same time"]
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
    parser.add_argument('--out-fq1', help='Output fastq file, SE or PE pair 1', type=str)
    parser.add_argument('--out-fq2', help='Output fastq file, PE pair 2', type=str)
    parser.add_argument('--out-fq12', help='Output merged fastq file, PE pairs 1 and 2', type=str)
    parser.add_argument('--html', help='Output html file', required=True, type=str)
    parser.add_argument('--json', help='Output json file', required=True, type=str)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error("\n".join(error_message_list))

    return vars(args)


def move_file(src_file: str, dst_file: str) -> None:
    """Move file from source to the destination"""
    shutil.move(src_file, dst_file)
    logger.info(f'{src_file} moved to {dst_file}')


def check_if_exist(file: str) -> None:
    """Checks if file exists"""
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    if os.path.getsize(file) == 0:
        logger.critical(f'Output file {file} is an empty, exiting..')
        sys.exit(1)
    logger.info(f"Expected file {file} found.")


def print_error_message(error_message: Optional[str]) -> None:
    if not error_message:
        error_message = '< EMPTY >'
    logger.warning(f'stderr: {error_message}')


def run_and_check_with_message(
    cmd: List[str],
    fail_message: str,
    exit_on_error: bool = True,
    return_proc: bool = False,
    **subprocess_args
) -> Optional[subprocess.CompletedProcess[str]]:
    logger.info(f"Running command {' '.join(cmd)}")
    if 'stderr' not in subprocess_args:
        subprocess_args['stderr'] = subprocess.PIPE
    try:
        proc = subprocess.run(cmd, text=True, check=True, **subprocess_args)
        print_error_message(proc.stderr)
        if return_proc:
            return proc
    except subprocess.CalledProcessError as e:
        logger.critical(f'{fail_message} failed with code {e.returncode}')
        print_error_message(e.stderr)
        if exit_on_error:
            logger.critical(f'{exit_on_error=}, now exiting.')
            sys.exit(1)


def run_fastp(in_fq1: str, in_fq2: str, disable_filters: list[str], merge: bool,
              json: str, html: str) -> tuple[str, str, str]:
    out_fq1, out_fq2, out_fq12 = (tempfile.NamedTemporaryFile(suffix=".fastq.gz").name for _ in range(3))

    cmd = ['fastp', '-i', in_fq1]
    cmd += ['-o', out_fq1]

    if in_fq2:
        cmd += ['-I', in_fq2, '-O', out_fq2, '--detect_adapter_for_pe']

    if merge:
        cmd += ['--merge', '--merged_out', out_fq12]

    cmd += [f'--disable_{mode}' for mode in disable_filters] if disable_filters else []
    cmd += ['--thread', str(os.cpu_count()), '--html', html, '--json', json]

    run_and_check_with_message(cmd, 'fastp')

    return out_fq1, out_fq2, out_fq12


def save_final_fastq_by_mode(merge_reads: bool, mock_merge_reads: bool, in_fq2_path: str,
                             processed_fq1_path: str, processed_fq2_path: str, processed_fq12_path: str,
                             out_fq1_path: str, out_fq2_path: str, out_fq12_path: str) -> None:
    """Saves final FASTQ files by provided modes"""

    logger.info("Going to move final FASTQ file(s) from /tmp to /output dir...")

    if mock_merge_reads:
        move_file(processed_fq12_path, out_fq12_path)
        logger.info("All files has been successfully moved from /tmp to /output dir.")
        return

    move_file(processed_fq1_path, out_fq1_path)

    if in_fq2_path:
        move_file(processed_fq2_path, out_fq2_path)

    if merge_reads:
        move_file(processed_fq12_path, out_fq12_path)

    logger.info("All files has been successfully moved from /tmp to /output dir.")


def main() -> None:
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {args}")

    merge_reads = True if args.mock_merge_reads else args.merge_reads

    tmp_fq1, tmp_fq2, tmp_fq12 = run_fastp(args.in_fq1, args.in_fq2, args.disable_filters, merge_reads,
                                           args.json, args.html)

    if args.mock_merge_reads:
        mock_merge_by_chunks(tmp_fq1, tmp_fq2, args.inner_distance_size, args.reads_chunk_size, tmp_fq12)

    save_final_fastq_by_mode(merge_reads, args.mock_merge_reads, args.in_fq2, tmp_fq1, tmp_fq2,
                             tmp_fq12, args.out_fq1, args.out_fq2, args.out_fq12)

    logger.info("Run is completed successfully.")


if __name__ == "__main__":
    main()
