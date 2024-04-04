import argparse
import logging
import subprocess
import os
import sys


LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger(name=__file__)


def set_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input fastq.gz file, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2', help='Input fastq.gz file, PE pair 2')
    parser.add_argument('--trimq',
                        help='Reads trimming: minimal base quality for 5'' & 3'' read ends to keep them', type=int)
    parser.add_argument('--disable', help='Mode(s) to disable', nargs='+', action='extend')
    parser.add_argument('--merge', help='Merge reads', action='store_true')
    parser.add_argument('--out-fq1', help='Output fastq file, SE or PE pair 1')
    parser.add_argument('--out-fq2', help='Output fastq file, PE pair 2')
    parser.add_argument('--out-fq12', help='Output merged fastq file, PE pairs 1 and 2')
    parser.add_argument('--out-html', help='Output html file', required=True)
    parser.add_argument('--out-json', help='Output json file', required=True)

    return parser.parse_args()


def run_command(command: list[str], stdin=None):
    logger.info(f'Running {command}...')

    try:
        command_process = subprocess.run(command, text=True, capture_output=True, input=stdin, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(f"Failed to run '{' '.join(command)}'")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"Undefined error: {e}")
        sys.exit(1)

    logger.info(command_process.stdout)
    logger.info(command_process.stderr)


def run(args: argparse.Namespace) -> None:

    cmd = ['fastp', '-i', args.in_fq1]
    cmd += ['-o', args.out_fq1]

    if args.in_fq2:
        cmd += ['-I', args.in_fq2, '-O', args.out_fq2, '--detect_adapter_for_pe']

    if args.merge:
        cmd += ['--merge', '--merged_out', args.out_fq12]

    cmd += ['--cut_tail', '--cut_front', '--cut_window_size', str(1), '--cut_mean_quality', str(args.trimq)] if args.trimq else []
    cmd += [f'--disable_{mode}' for mode in args.disable] if args.disable else []
    cmd += ['--thread', str(os.cpu_count()), '--html', args.out_html, '--json', args.out_json]

    run_command(cmd)


if __name__ == "__main__":
    set_logger(LOGGER_FORMAT)

    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)
    logger.info("Run is completed successfully.")
