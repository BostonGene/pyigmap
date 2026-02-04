import argparse
import glob
import shutil
import subprocess
import sys
import tempfile
import os
from typing import Optional, List
import logging

LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger()

IGBLAST_DIR = os.environ.get('IGBLAST_DIR')

IMGT_URL = "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"

SPECIE = 'Homo_sapiens'

ORGANISMS_BLACKLIST = ["mouse", "rabbit", "rat", "rhesus_monkey"]

CHAINS_GLOSSARY = {
    'Ig.V': ['IGHV', 'IGKV', 'IGLV'],
    'Ig.J': ['IGHJ', 'IGKJ', 'IGLJ'],
    'Ig.D': ['IGHD'],
    'TCR.V': ['TRAV', 'TRBV', 'TRDV', 'TRGV'],
    'TCR.J': ['TRAJ', 'TRBJ', 'TRDJ', 'TRGJ'],
    'TCR.D': ['TRBD', 'TRDD']
}

REFERENCE_DIR = os.path.join(IGBLAST_DIR, 'database')


def configure_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--all-alleles', action='store_true',
                        help='Will use all alleles provided in the antigen receptor segment database '
                             '(*01, *02, etc. according to IMGT).')
    parser.add_argument('-o', '--out-archive', help='Output IgBLAST reference archive path', required=True)

    return parser.parse_args()


def exit_with_error(message: Optional[str]) -> None:
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


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


def check_if_exist(file: str) -> None:
    """Checks if file exists"""
    if not os.path.exists(file):
        exit_with_error(f"Expected file {file} not found, exiting...")
    logger.info(f"Expected file {file} found.")


def download_fasta_from_imgt_database(chains_list: list[str], specie: str) -> list[str]:
    """Downloads VDJ fasta reference from https://www.imgt.org/"""
    fasta_local_paths = []
    for chain in chains_list:
        fasta_link = f"{IMGT_URL}{SPECIE}/{chain[:2]}/{chain}.fasta"
        fasta_local_path = tempfile.NamedTemporaryFile().name
        cmd = ['wget', fasta_link, '-q', '-O', fasta_local_path]
        run_and_check_with_message(cmd, "wget")
        fasta_local_paths.append(fasta_local_path)
    return fasta_local_paths


def filter_minor_alleles(fasta_path: str) -> str:
    """Filters out minor alleles: *02, *03, etc."""
    filtered_fasta_path = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'grep', fasta_path, '-r', '-p', '\\*01', '-o', filtered_fasta_path]
    run_and_check_with_message(cmd, "seqkit")
    return filtered_fasta_path


def concat_fasta_files(fasta_paths: list[str]) -> str:
    """Concatenates fasta files into one file"""
    concatenated_fasta_path = tempfile.NamedTemporaryFile().name
    with open(concatenated_fasta_path, 'a') as out_f:
        for file in fasta_paths:
            with open(file, 'r') as f:
                out_f.write(f.read())
                logger.info(f'{file} appended to {concatenated_fasta_path}!')
    check_if_exist(concatenated_fasta_path)
    logger.info(f'Output {concatenated_fasta_path} file has been successfully wrote.')
    return concatenated_fasta_path


def clean_imgt_fasta(fasta_path: str) -> str:
    """Cleans the header of the fasta downloaded from IMGT"""
    cmd = [os.path.join(IGBLAST_DIR, 'bin', 'edit_imgt_file.pl'), fasta_path]
    edit_imgt_file_process = run_and_check_with_message(cmd, "edit_imgt_file.pl", return_proc=True,
                                                        capture_output=True, stderr=None)
    clean_fasta_path = tempfile.NamedTemporaryFile().name
    with open(clean_fasta_path, 'w') as f:
        f.write(edit_imgt_file_process.stdout)
    return clean_fasta_path


def make_blast_db(fasta_path: str, output_basename: str) -> None:
    """Makes the blast database from cleaned IMGT fasta"""
    makeblastdb_bin = os.path.join(IGBLAST_DIR, 'bin', 'makeblastdb')
    cmd = [makeblastdb_bin, '-parse_seqids', '-dbtype', 'nucl', '-in', fasta_path, '-out', output_basename]
    run_and_check_with_message(cmd, 'makeblastdb')


def archive_reference_as_tar_gz(archive_path: str) -> str:
    """Makes tar.gz archive with IgBLAST final reference"""

    cmd = ['tar', '-czf', archive_path, '-C', IGBLAST_DIR, 'database', 'internal_data', 'optional_file']
    run_and_check_with_message(cmd, "tar")
    check_if_exist(archive_path)
    logger.info(f"'{archive_path}' has been successfully created.")

    return archive_path


def remove_duplicates_by_sequence_id(fasta_path: str) -> str:
    """Removes duplicated sequences by sequence id (header).
    If you don't do this, makeblastdb tool will fall
    """
    output_fasta = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'rmdup', '--by-name', fasta_path, '-o', output_fasta]
    run_and_check_with_message(cmd, "seqkit")
    return output_fasta


def remove_path(path: str) -> None:
    try:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)
    except Exception as e:
        logger.error(f"Error removing {path}: {e}")


def remove_unnecessary_organisms() -> None:
    for organism in ORGANISMS_BLACKLIST:
        paths_to_remove = glob.glob(os.path.join(IGBLAST_DIR, "**", f"{organism}*"), recursive=True)
        for path in paths_to_remove:
            if os.path.exists(path):
                remove_path(path)


def main() -> None:
    configure_logger()

    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    for library, chains_list in CHAINS_GLOSSARY.items():

        fasta_paths = download_fasta_from_imgt_database(chains_list, SPECIE)

        if not args.all_alleles:
            fasta_paths = [filter_minor_alleles(fasta_path) for fasta_path in fasta_paths]

        clean_fasta_paths = [clean_imgt_fasta(fasta_path) for fasta_path in fasta_paths]

        concatenated_fasta_path = concat_fasta_files(clean_fasta_paths)

        fasta_without_duplicates_path = remove_duplicates_by_sequence_id(concatenated_fasta_path)

        output_basename = os.path.join(REFERENCE_DIR, f'{SPECIE}.{library}')
        make_blast_db(fasta_without_duplicates_path, output_basename)

    remove_unnecessary_organisms()
    archive_reference_as_tar_gz(args.out_archive)

    logger.info("Run is completed successfully.")


if __name__ == "__main__":
    main()
