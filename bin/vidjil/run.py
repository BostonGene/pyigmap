import argparse
import glob
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import typing
from typing import Optional, Union, List

LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger()

TEMPDIR_NAME = tempfile.gettempdir()

ORGANISM_GLOSSARY = {
    'human': 'homo-sapiens',
    'rat': 'rattus-norvegicus',
    'mouse': 'mus-musculus'
}

CPU_COUNT = os.cpu_count()

REF_DIR = tempfile.TemporaryDirectory().name


def configure_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args() -> argparse.Namespace:
    """Parses arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in-fastq',
        help='Input FASTQ'
    )
    parser.add_argument(
        '--out-fasta',
        help='Output FASTA with detected CDR3',
        required=True,
        type=str
    )
    parser.add_argument(
        '--vdj-ref',
        help='V(D)J reference',
        required=True
    )
    parser.add_argument(
        '--reads-chunk-size',
        type=int,
        default=5_000_000
    )
    parser.add_argument(
        '--organism',
        help="Organism name: 'human', 'rat' or 'mouse'",
        choices=["human", "rat", "mouse"],
        default="human"
    )
    parser.add_argument(
        '--logs',
        help='Output logs file',
        type=str
    )
    parser.add_argument(
        '--debug',
        help='Enables saving logs',
        action="store_true"
    )

    return parser.parse_args()


def exit_with_error(message: Optional[str]) -> None:
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


def read_fastq_file_chunk(file_obj: typing.TextIO, reads_chunk_size: int) -> str:
    """Reads a FASTQ file chunk and returns a concatenated string of reads."""
    reads_list = []
    while len(reads_list) < reads_chunk_size and (read_header := file_obj.readline()):
        reads_list.append(read_header + file_obj.readline() + file_obj.readline() + file_obj.readline())
    return ''.join(reads_list)


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
        logger.critical(f"{fail_message} failed with code {e.returncode}")
        print_error_message(e.stderr)
        if exit_on_error:
            logger.critical(f"{exit_on_error=}, now exiting.")
            sys.exit(1)


def check_if_exist_and_not_empty(file: str) -> None:
    """Checks if file exists"""
    if not os.path.exists(file):
        exit_with_error(f"Expected file {file} not found, exiting...")
    if os.path.getsize(file) == 0:
        exit_with_error(f"Output file {file} is an empty, exiting...")
    logger.info(f"Expected file {file} successfully found.")


def move_file(src_file: str, dst_file: str) -> None:
    """Moves file from source to the destination"""
    shutil.move(src_file, dst_file)
    logger.info(f"{src_file} moved to {dst_file}")


def get_file_obj(file_path: str, mode: str) -> Union[typing.TextIO, gzip.GzipFile]:
    """Returns reader or writer file object"""
    if 'b' in mode:
        return gzip.open(file_path, mode)
    return open(file_path, mode)


def concat_files(output_file_path: str, files_to_concatenate: list[str], file_type: str) -> None:
    """Concatenates files into one file"""
    with get_file_obj(output_file_path, "a" + file_type) as write_obj:
        for file in files_to_concatenate:
            with get_file_obj(file, "r" + file_type) as read_obj:
                write_obj.write(read_obj.read())
                logger.info(f"{file} appended to {output_file_path}!")
    check_if_exist_and_not_empty(output_file_path)
    logger.info(f"Output {output_file_path} file successfully wrote.")


def concat_and_move_file(out_file: str, files_to_concatenate: list[str], file_type: str = "") -> None:
    """Prepares output file(s)"""
    temp_out_file = os.path.join(TEMPDIR_NAME, os.path.basename(out_file))
    concat_files(temp_out_file, files_to_concatenate, file_type)
    move_file(temp_out_file, out_file)


def save_vidjil_results(out_fasta: str, logs: str, debug: bool) -> None:
    """Save Vidjil results"""
    fasta_files_paths = glob.glob(os.path.join(TEMPDIR_NAME, '*detected.vdj.fa.gz'))
    concat_and_move_file(out_fasta, fasta_files_paths, file_type="b")
    if debug:
        log_files_paths = glob.glob(os.path.join(TEMPDIR_NAME, '*.vidjil.log'))
        concat_and_move_file(logs, log_files_paths)


def run_vidjil_in_parallel(command: list, stdin: str) -> None:
    """Parallels vidjil calculation"""
    parallel_cmd = ['parallel', '-j', str(CPU_COUNT), '--pipe', '-L', '4', '--round-robin'] + command
    vidjil_log_file = tempfile.NamedTemporaryFile(suffix=".vidjil.log").name
    with open(vidjil_log_file, 'w') as file_obj:
        run_and_check_with_message(parallel_cmd, 'vidjil', input=stdin, stderr=None, stdout=file_obj)


def detect_cdr3_in_reads(reads_chunk: str, output_basename: str, organism: str) -> None:
    """Detects CDR3 in FASTQ file reads"""
    logger.info("Starting detect reads with CDR3...")
    vidjil_cmd = [
        'vidjil-algo',
        '--germline', os.path.join(REF_DIR, f'{organism}.g'),
        '--out-detected',
        '--dir', TEMPDIR_NAME,
        '--gz',
        '--clean-memory',
        '--base', output_basename,
        '-c', 'detect',
        '-'
    ]

    run_vidjil_in_parallel(vidjil_cmd, reads_chunk)
    logger.info("All reads with CDR3 successfully detected.")


def run_vidjil_by_fastq_chunk(fastq_file: str, reads_chunk_size: int, organism: str) -> None:
    with gzip.open(fastq_file, "rt") as fq_file_obj:
        chunk_number = 1
        while reads_chunk := read_fastq_file_chunk(fq_file_obj, reads_chunk_size):
            output_basename = f"{os.path.basename(fastq_file)}.chunk{chunk_number}.part{{#}}"
            detect_cdr3_in_reads(reads_chunk, output_basename, organism)
            chunk_number += 1


def unpack_reference(archive: str) -> None:
    """Unpack V(D)J reference into selected directory"""
    logger.info(f'Unpacking {archive} into {REF_DIR}...')
    os.mkdir(REF_DIR)
    tar_cmd = ['tar', '-xvf', archive, '-C', REF_DIR]
    run_and_check_with_message(tar_cmd, 'tar')
    logger.info(f'V(D)J reference {archive} successfully unpacked.')


def main() -> None:
    configure_logger()

    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    unpack_reference(args.vdj_ref)

    organism = ORGANISM_GLOSSARY.get(args.organism)
    run_vidjil_by_fastq_chunk(args.in_fastq, args.reads_chunk_size, organism)

    save_vidjil_results(args.out_fasta, args.logs, args.debug)

    logger.info("Run is completed successfully.")


if __name__ == "__main__":
    main()
