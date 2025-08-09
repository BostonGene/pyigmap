import argparse
import subprocess
import os
import sys
import tempfile
from typing import Optional, List
from logger import set_logger

from mock_merge_rs import mock_merge_by_chunks

logger = set_logger(name=__file__)

CPU_COUNT = os.cpu_count()
TEMPDIR_NAME = os.environ.get('TEMPDIR', '/tmp/')
CHUNK_SIZE = 20 * 1024 * 1024  # 20MB


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if args.merge_reads and args.mock_merge_reads:
        msg_list += ["--merge-reads and --mock-merge-reads cannot be used at the same time"]
    return msg_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input fastq.gz file, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2', help='Input fastq.gz file, PE pair 2')
    parser.add_argument('--disable-filters', help='Mode(s) to disable',
                        choices=["length_filtering", "quality_filtering", "adapter_trimming"], nargs='+', default=[],
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

    return args


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


def generate_temp_fastq_files() -> tuple[str, str, str]:
    """Generate three temporary FASTQ file paths in the temp directory."""
    return tuple(
        os.path.join(TEMPDIR_NAME, os.path.basename(tempfile.NamedTemporaryFile(suffix=".fastq").name))
        for _ in range(3)
    )

def run_fastp(in_fq1: str, in_fq2: Optional[str], disable_filters: list[str], merge: bool,
              json: str, html: str) -> tuple[str, str, str]:
    out_fq1, out_fq2, out_fq12 = generate_temp_fastq_files()

    cmd = ['fastp', '-i', in_fq1, '-o', out_fq1]

    if in_fq2:
        cmd += ['-I', in_fq2, '--detect_adapter_for_pe', '-O', out_fq2]

    if merge and in_fq2:
        cmd += ['--merge', '--merged_out', out_fq12]

    cmd += [f'--disable_{mode}' for mode in disable_filters]
    cmd += ['--thread', str(CPU_COUNT), '--html', html, '--json', json]

    run_and_check_with_message(cmd, 'fastp')
    return out_fq1, out_fq2, out_fq12


def compress_to_gzip(src_file: str, dst_file: str):
    with open(dst_file, "w") as f_out:
        run_and_check_with_message(
            ["pigz", "-c", src_file],
            "pigz",
            stdout=f_out
        )
    logger.info(f"Compressed '{src_file}' to '{dst_file}' successfully.")


def save_final_fastq(merge_reads: bool, mock_merge_reads: bool, paired_end: bool,
                     fq1_src: str, fq2_src: str, fq12_src: str,
                     fq1_dest: str, fq2_dest: str, fq12_dest: str) -> None:
    logger.info("Compressing and moving final FASTQ file(s) to output dir...")

    if paired_end:
        if merge_reads or mock_merge_reads:
            compress_to_gzip(fq12_src, fq12_dest)
        if merge_reads:
            compress_to_gzip(fq1_src, fq1_dest)
            compress_to_gzip(fq2_src, fq2_dest)
    else:
        if merge_reads or mock_merge_reads:
            compress_to_gzip(fq1_src, fq12_dest)
        else:
            compress_to_gzip(fq1_src, fq1_dest)

    logger.info("All files have been successfully compressed and moved.")


def main() -> None:
    args = parse_args()
    logger.info(f"Starting with arguments: {args}")

    paired_end = bool(args.in_fq2)
    run_reads_merging = (args.mock_merge_reads or args.merge_reads) and paired_end

    tmp_fq1, tmp_fq2, tmp_fq12 = run_fastp(args.in_fq1, args.in_fq2, args.disable_filters, run_reads_merging,
                                           args.json, args.html)

    if args.mock_merge_reads and paired_end:
        mock_merge_by_chunks(tmp_fq1, tmp_fq2, args.inner_distance_size, args.reads_chunk_size, tmp_fq12)

    save_final_fastq(args.merge_reads, args.mock_merge_reads, paired_end,
                     tmp_fq1, tmp_fq2, tmp_fq12,
                     args.out_fq1, args.out_fq2, args.out_fq12)

    logger.info("Run completed successfully.")


if __name__ == "__main__":
    main()
