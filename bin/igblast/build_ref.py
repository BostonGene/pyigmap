import argparse
import tempfile
import os

from utils import IGBLAST_DIR, run_command, check_if_exist
from logger import set_logger

logger = set_logger(name=__file__)

TMP_DIR = '/tmp'

SPECIES_GLOSSARY = {'human': 'Homo_sapiens',
                    'mouse': 'Mus_musculus'}

CHAINS_GLOSSARY = {'Ig.V': ['IGHV', 'IGKV', 'IGLV'],
                   'Ig.J': ['IGHJ', 'IGKJ', 'IGLJ'],
                   'Ig.D': ['IGHD'],
                   'TCR.V': ['TRAV', 'TRBV', 'TRDV', 'TRGV'],
                   'TCR.J': ['TRAJ', 'TRBJ', 'TRDJ', 'TRGJ'],
                   'TCR.D': ['TRBD', 'TRDD']}

REFERENCE_DIR = os.path.join(IGBLAST_DIR, 'database')


def parse_args() -> argparse.Namespace:
    """Parses arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--all-alleles', action='store_true',
                        help='Will use all alleles provided in the antigen receptor segment database '
                             '(*01, *02, etc. according to IMGT).')
    parser.add_argument('-o', '--out-archive', help='Output IgBLAST reference archive path', required=True)

    return parser.parse_args()


def download_fasta_from_imgt_database(chains_list: list[str], specie: str) -> list[str]:
    """Downloads VDJ fasta reference from https://www.imgt.org/"""
    fasta_local_paths = []
    for chain in chains_list:
        specie_in_imgt_format = SPECIES_GLOSSARY[specie]
        fasta_link = (f"https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
                      f"{specie_in_imgt_format}/{chain[:2]}/{chain}.fasta")
        fasta_local_path = tempfile.NamedTemporaryFile().name
        cmd = ['wget', fasta_link, '-q', '-O', fasta_local_path]
        _ = run_command(cmd)
        fasta_local_paths.append(fasta_local_path)
    return fasta_local_paths


def filter_minor_alleles(fasta_path: str) -> str:
    """Filters out minor alleles: *02, *03, etc."""
    filtered_fasta_path = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'grep', fasta_path, '-r', '-p', '\\*01', '-o', filtered_fasta_path]
    cmd_process = run_command(cmd)
    logger.info(cmd_process.stderr)
    return filtered_fasta_path


def concat_fasta_files(fasta_paths: list[str]) -> str:
    """Concatenates fasta files into one file"""
    concatenated_fasta_path = tempfile.NamedTemporaryFile().name
    with open(concatenated_fasta_path, 'a') as out_f:
        for i, file in enumerate(fasta_paths):
            with open(file, 'r') as f:
                out_f.write(f.read())
                logger.info(f'{file} appended to {concatenated_fasta_path}!')
    check_if_exist(concatenated_fasta_path)
    logger.info(f'Output {concatenated_fasta_path} file has been successfully wrote.')
    return concatenated_fasta_path


def clean_imgt_fasta(fasta_path: str) -> str:
    """Cleans the header of the fasta downloaded from IMGT"""
    cmd = [os.path.join(IGBLAST_DIR, 'bin', 'edit_imgt_file.pl'), fasta_path]
    clean_fasta = run_command(cmd)
    clean_fasta_path = tempfile.NamedTemporaryFile().name
    with open(clean_fasta_path, 'w') as f:
        f.write(clean_fasta.stdout)
    return clean_fasta_path


def make_blast_db(fasta_path: str, output_basename: str):
    """Makes the blast database from cleaned IMGT fasta"""
    cmd = [os.path.join(IGBLAST_DIR, 'bin', 'makeblastdb'), '-parse_seqids', '-dbtype', 'nucl', '-in', fasta_path,
           '-out', output_basename]
    cmd_process = run_command(cmd)
    logger.info(cmd_process.stdout)


def archive_reference_as_tar_gz(archive_path: str) -> str:
    """Makes tar.gz archive with IgBLAST final reference"""

    cmd = ['tar', '-czf', archive_path, '-C', IGBLAST_DIR, 'database', 'internal_data', 'optional_file']
    _ = run_command(cmd)
    check_if_exist(archive_path)
    logger.info(f"'{archive_path}' has been successfully created.")

    return archive_path


def remove_duplicates_by_sequence_id(fasta_path: str) -> str:
    """Removes duplicated sequences by sequence id (header).
    If you don't do this, makeblastdb tool will fall
    """
    output_fasta = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'rmdup', '--by-name', fasta_path, '-o', output_fasta]
    cmd_process = run_command(cmd)
    logger.info(cmd_process.stderr)
    return output_fasta


def run(args: argparse.Namespace):
    for specie in SPECIES_GLOSSARY:
        for library, chains_list in CHAINS_GLOSSARY.items():

            fasta_paths = download_fasta_from_imgt_database(chains_list, specie)

            if not args.all_alleles:
                fasta_paths = [filter_minor_alleles(fasta_path) for fasta_path in fasta_paths]

            clean_fasta_paths = [clean_imgt_fasta(fasta_path) for fasta_path in fasta_paths]

            concatenated_fasta_path = concat_fasta_files(clean_fasta_paths)

            fasta_without_duplicates_path = remove_duplicates_by_sequence_id(concatenated_fasta_path)

            output_basename = os.path.join(REFERENCE_DIR, f'{specie}.{library}')
            make_blast_db(fasta_without_duplicates_path, output_basename)

    archive_reference_as_tar_gz(args.out_archive)


if __name__ == "__main__":
    args = parse_args()

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
