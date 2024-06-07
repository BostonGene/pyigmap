import argparse
import gzip
import shutil
import subprocess
import os
import sys
import tempfile

from mock_merge import mock_merge_reads

from logger import set_logger

logger = set_logger(name=__file__)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if args.insert_size and not args.mock_merge:
        msg_list += ["--insert-size cannot be provided without --mock-merge"]
    if args.mock_merge and (args.out_fq1 or args.out_fq2):
        msg_list += ["--out-fq1 or --out-fq2 cannot be used with --mock-merge"]
    return msg_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input fastq.gz file, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2', help='Input fastq.gz file, PE pair 2')
    parser.add_argument('--trimq',
                        help='Reads trimming: minimal base quality for 5'' & 3'' read ends to keep them', type=int)
    parser.add_argument('--disable', help='Mode(s) to disable',
                        choices=["length_filtering", "quality_filtering", "adapter_trimming"], nargs='+', action='extend')
    parser.add_argument('--merge', help='Merge reads', action='store_true')
    parser.add_argument('--mock-merge', help='Enable mock merging of reads', action='store_true')
    parser.add_argument('--insert-size', help='Insert size for mock merging', type=int)
    parser.add_argument('--out-fq1', help='Output fastq file, SE or PE pair 1')
    parser.add_argument('--out-fq2', help='Output fastq file, PE pair 2')
    parser.add_argument('--out-fq12', help='Output merged fastq file, PE pairs 1 and 2')
    parser.add_argument('--out-html', help='Output html file', required=True)
    parser.add_argument('--out-json', help='Output json file', required=True)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error("\n".join(error_message_list))

    return args


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


def run_fastp(in_fq1: str, in_fq2: str, trimq: int, disable: str, merge: bool, out_json: str, out_html: str,
              out_fq1: str, out_fq2: str, out_fq12: str):
    cmd = ['fastp', '-i', in_fq1]
    cmd += ['-o', out_fq1]

    if in_fq2:
        cmd += ['-I', in_fq2, '-O', out_fq2, '--detect_adapter_for_pe']

    if merge:
        cmd += ['--merge', '--merged_out', out_fq12]

    cmd += ['--cut_tail', '--cut_front', '--cut_window_size', str(1), '--cut_mean_quality',
            str(trimq)] if trimq else []
    cmd += [f'--disable_{mode}' for mode in disable] if disable else []
    cmd += ['--thread', str(os.cpu_count()), '--html', out_html, '--json', out_json]

    run_command(cmd)


def concat_gz_files(gz_files: list[str]) -> str:
    """Concatenates gz files into one file"""
    out_file_path = tempfile.NamedTemporaryFile(suffix=".gz").name
    with gzip.open(out_file_path, "ab") as concat_file:
        for gz_file in gz_files:
            with gzip.open(gz_file, "rb") as f:
                concat_file.write(f.read())
                logger.info(f"{gz_file} appended to {out_file_path}")
    check_if_exist(out_file_path)
    logger.info(f"{out_file_path} concatenation has been done.")
    return out_file_path


def run(args: argparse.Namespace) -> None:

    merge = True if args.mock_merge else args.merge

    out_fq1 = tempfile.NamedTemporaryFile(suffix=".gz").name
    out_fq2 = tempfile.NamedTemporaryFile(suffix=".gz").name
    out_fq12 = tempfile.NamedTemporaryFile(suffix=".gz").name

    run_fastp(args.in_fq1, args.in_fq2, args.trimq, args.disable, merge, args.out_json, args.out_html,
              out_fq1, out_fq2, out_fq12)

    if args.in_fq2 and args.mock_merge:
        fq12_mock = mock_merge_reads(out_fq1, out_fq2, args.insert_size)
        out_fq12 = concat_gz_files([out_fq12, fq12_mock])

    if args.out_fq1:
        replace_file(out_fq1, args.out_fq1)

    if args.out_fq2:
        replace_file(out_fq2, args.out_fq2)

    if args.out_fq12:
        replace_file(out_fq12, args.out_fq12)


if __name__ == "__main__":
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)
    logger.info("Run is completed successfully.")
