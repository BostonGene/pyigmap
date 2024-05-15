import argparse
import tempfile
import os

from utils import IGBLAST_DIR, run_command, check_if_exist
from logger import set_logger

logger = set_logger(name=__file__)

TMP_DIR = '/tmp'

SPECIES_GLOSSARY = {'human': 'Homo_sapiens',
                    'mouse': 'Mus_musculus'}

LOCI_GLOSSARY = {'Ig.V': ['IGHV', 'IGKV', 'IGLV'],
                 'Ig.J': ['IGHJ', 'IGKJ', 'IGLJ'],
                 'Ig.D': ['IGHD'],
                 'TCR.V': ['TRAV', 'TRBV', 'TRDV', 'TRGV'],
                 'TCR.J': ['TRAJ', 'TRBJ', 'TRDJ', 'TRGJ'],
                 'TCR.D': ['TRBD', 'TRDD']}

REFERENCE_DIR = os.path.join(IGBLAST_DIR, 'database')


def parse_args() -> argparse.Namespace:
    """
    Parses run.py script arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allow-minor-alleles', action='store_true', help='Allow only minor alleles (*02, *03, etc.)')
    parser.add_argument('-o', '--out-archive', help='Output IgBLAST reference archive basename', required=True)

    return parser.parse_args()


def download_fasta_from_imgt(gene: str, specie: str) -> list[str]:
    fasta_local_paths = []
    loci_list = LOCI_GLOSSARY[gene]
    for locus in loci_list:
        specie_imgt = SPECIES_GLOSSARY[specie]
        fasta_link = f"https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{specie_imgt}/{locus[:2]}/{locus}.fasta"
        fasta_local_path = tempfile.NamedTemporaryFile().name
        cmd = ['wget', fasta_link, '-O', fasta_local_path]
        _ = run_command(cmd)
        fasta_local_paths.append(fasta_local_path)
    return fasta_local_paths


def filter_minor_alleles(fasta_path: str) -> str:
    filtered_fasta_path = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'grep', fasta_path, '-r', '-p', '\\*01', '-o', filtered_fasta_path]

    seqkit_logs = run_command(cmd)
    logger.info(seqkit_logs.stdout)

    return filtered_fasta_path


def concat_fasta_files(fasta_paths: list[str]) -> str:
    """Concatenates fasta files into one file"""
    concatenated_fasta_path = tempfile.NamedTemporaryFile().name
    writer = open(concatenated_fasta_path, 'a')
    with writer as out_f:
        for i, file in enumerate(fasta_paths):
            reader = open(file, 'r')
            with reader as f:
                out_f.write(f.read())
                logger.info(f'{file} appended to {concatenated_fasta_path}!')
    check_if_exist(concatenated_fasta_path)
    logger.info(f'Output {concatenated_fasta_path} file has been successfully wrote.')
    return concatenated_fasta_path


def clean_imgt_fasta(fasta_path: str) -> str:
    cmd = [os.path.join(IGBLAST_DIR, 'bin', 'edit_imgt_file.pl'), fasta_path]
    edited_fasta_stdout = run_command(cmd).stdout
    edited_fasta_path = tempfile.NamedTemporaryFile().name
    with open(edited_fasta_path, 'w') as f:
        f.write(edited_fasta_stdout)
    return edited_fasta_path


def make_blast_db(fasta_path: str, output_basename: str):
    cmd = [os.path.join(IGBLAST_DIR, 'bin', 'makeblastdb'), '-parse_seqids', '-dbtype', 'nucl', '-in', fasta_path,
           '-out', output_basename]
    makeblastdb_logs = run_command(cmd).stdout
    logger.info(makeblastdb_logs)


def archive_reference_as_tar_gz(archive_path: str) -> str:
    """Makes tar.gz archive"""

    cmd = ['tar', '-czf', archive_path, '-C', IGBLAST_DIR, os.path.join(IGBLAST_DIR, 'database'),
           os.path.join(IGBLAST_DIR, 'internal_data'), os.path.join(IGBLAST_DIR, 'optional_file')]
    _ = run_command(cmd)
    check_if_exist(archive_path)
    logger.info(f"'{archive_path}' has been successfully created.")

    return archive_path


def remove_duplicates_by_id(fasta_path: str) -> str:
    output_fasta = tempfile.NamedTemporaryFile().name
    cmd = ['seqkit', 'rmdup', '--by-name', fasta_path, '-o', output_fasta]
    seqkit_logs = run_command(cmd).stdout
    logger.info(seqkit_logs)
    return output_fasta


def run(args: argparse.Namespace):
    for specie in SPECIES_GLOSSARY:
        for gene in LOCI_GLOSSARY:
            fasta_paths = download_fasta_from_imgt(gene, specie)
            if not args.allow_minor_alleles:
                fasta_paths = [filter_minor_alleles(fasta_path) for fasta_path in fasta_paths]
            edited_fasta_paths = [clean_imgt_fasta(fasta_path) for fasta_path in fasta_paths]
            concatenated_fasta_path = concat_fasta_files(edited_fasta_paths)
            fasta_without_duplicates_path = remove_duplicates_by_id(concatenated_fasta_path)
            output_basename = os.path.join(REFERENCE_DIR, f'{specie}.{gene}')
            make_blast_db(fasta_without_duplicates_path, output_basename)
    archive_reference_as_tar_gz(os.path.join(TMP_DIR, args.out_archive))


if __name__ == "__main__":
    args = parse_args()

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
