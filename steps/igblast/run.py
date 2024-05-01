import argparse
import gzip
import logging
import os
from tqdm import tqdm
import subprocess
import sys
import shutil

from logger import set_logger, TqdmToLogger

logger = set_logger(name=__file__)
tqdm_out = TqdmToLogger(logger, level=logging.INFO)

IGBLAST_DIR = os.environ.get('IGBLAST_DIR')
TMP_DIR = '/tmp'

RECEPTOR_GLOSSARY = {'BCR': 'Ig',
                     'TCR': 'TCR',
                     'all': 'all'}
ORGANISM_GLOSSARY = {'human': 'human',
                     'mouse': 'mouse'}

FASTA_CHUNKS_NUM = 50000
FASTA_CHUNKS_DIR = os.path.join(TMP_DIR, "fasta_chunks")


def parse_args() -> argparse.Namespace:
    """
    Parses run.py script arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input with the forward fastq')
    parser.add_argument('--in-fq2', help='Input with the reverse fastq')
    parser.add_argument('--in-fq12', help='Input with the merged (forward + reverse) fastq')
    parser.add_argument('--in-fasta', help='Input fasta file with detected V(D)J segments')
    parser.add_argument('--receptor', help="Receptor type: 'BCR' or 'TCR'", required=True)
    parser.add_argument('--organism', help="Organism name: 'human' or 'mouse'", default='human')
    parser.add_argument('--in-ref', help='FASTA reference with V(D)J segments', required=True)
    parser.add_argument('--out-annotation', help='Output BCR annotation table', required=True)

    return parser.parse_args()


def run_command(command: list[str], stdin=None) -> subprocess.CompletedProcess[str]:
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

    return command_process


def decompress(archive: str) -> None:
    """
    Decompresses archive into selected directory.

    :param archive: path to the zip archive
    """
    logger.info(f'Unzipping {archive} into {IGBLAST_DIR}...')
    tar_cmd = ['tar', '-I', 'pigz', '-xvf', archive, '-C', IGBLAST_DIR]
    _ = run_command(tar_cmd)
    logger.info(f'V(D)J reference {archive} is successfully unzipped!')


def read_gz(file_path: str) -> str:
    """
    Reads .gz archive contents into stdout.

    :param file_path: path to the .gz file

    :return: decompressed file stdout
    """
    logger.info(f"Reading {file_path} file...")
    pigz_cmd = ['pigz', '-dc', file_path]
    pigz_stdout = run_command(pigz_cmd).stdout
    logger.info(f"{file_path} file is successfully read.")

    return pigz_stdout


def convert_fastq_to_fasta(fastq_stdin: str) -> str:
    """
    Converts fastq to fasta.

    param: fastq_stdin: fastq stdout

    :return: fasta stdout
    """
    logger.info("Converting fastq -> fasta...")
    seqtk_cmd = ['seqtk', 'seq', '-a']
    sed_stdout = run_command(seqtk_cmd, fastq_stdin).stdout
    logger.info("Converting fastq -> fasta has been done.")

    return sed_stdout


def generate_annotation(fasta_files: list[str], receptor: str, organism: str) -> list[str]:
    """
    Generates IgBLAST annotation

    :param fasta_files: fasta paths
    :param receptor: receptor name: TCR or Ig
    :param organism: organism

    :return: annotation stdout
    """
    output_annotations = []
    for fasta_file in tqdm(fasta_files, file=tqdm_out, mininterval=1):
        logger.info(f'Generating annotation for {fasta_file}...')

        output_annotation = os.path.join(TMP_DIR, os.path.basename(fasta_file) + f'.{receptor}.tsv')

        igblast_cmd = ['bin/igblastn',
                       '-query', fasta_file,
                       '-germline_db_V', os.path.join(IGBLAST_DIR, 'database', f'human.{receptor}.V'),
                       '-germline_db_D', os.path.join(IGBLAST_DIR, 'database', f'human.{receptor}.D'),
                       '-germline_db_J', os.path.join(IGBLAST_DIR, 'database', f'human.{receptor}.J'),
                       '-organism', organism,
                       '-auxiliary_data', os.path.join(IGBLAST_DIR, 'optional_file', 'human_gl.aux'),
                       '-ig_seqtype', receptor,
                       '-show_translation',
                       '-outfmt', str(19),
                       '-num_threads', str(os.cpu_count()),
                       '-out', output_annotation]
        igblast_cmd += ['-c_region_db',
                        os.path.join(IGBLAST_DIR, 'database', 'ncbi_human_c_genes')] if organism == 'human' else []
        _ = run_command(igblast_cmd)
        output_annotations.append(output_annotation)
        logger.info(f'{output_annotation} annotation has been generated.')

    return output_annotations


def check_if_exist(file: str) -> None:
    """
    Checks if file exists

    :param file: path to the file
    """
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    logger.info(f"Expected file {file} found.")


def replace_file(src_file: str, dst_file: str) -> None:
    """
    Replaces file from source to the destination

    :param src_file: path to the file, that we need to move
    :param dst_file: path to destination file
    """
    shutil.move(src_file, dst_file)
    logger.info(f'{src_file} replaced to {dst_file}.')


def concat_annotation(annotation_files: list[str], out_annotation_path: str) -> None:
    """
    Concatenates temporary annotations into one file

    :param out_annotation_path: path to the output annotation
    :param annotation_files: list of calculated annotations paths
    """
    with gzip.open(out_annotation_path, 'ab') as out_annotation_obj:
        for i, annotation_file in enumerate(annotation_files):
            with open(annotation_file, 'r') as annotation_obj:
                if i > 0:
                    annotation_obj.readline()  # remove duplicated header line
                out_annotation_obj.write(annotation_obj.read().encode('utf-8'))
                logger.info(f'{annotation_file} appended to {out_annotation_path}')
    check_if_exist(out_annotation_path)
    logger.info(f'{out_annotation_path} concatenation has been done.')


def save_annotation(out_annotation_path: str, annotations_paths: list[str]):
    """
    Prepares and saves result annotation

    :param out_annotation_path: path to the output annotation
    :param annotations_paths: list of calculated annotations paths
    """
    temp_out_annotation_path = os.path.join(TMP_DIR, os.path.basename(out_annotation_path))
    logger.info(f'Defined {temp_out_annotation_path} as an output annotation file.')
    concat_annotation(annotations_paths, temp_out_annotation_path)
    replace_file(temp_out_annotation_path, out_annotation_path)


def split_by_chunks(fasta_stdin: str) -> list[str]:
    """
    Splits fasta by chunks
    """
    logger.info(f'Splitting fasta into chunks by {FASTA_CHUNKS_NUM} reads...')
    seqkit_cmd = ['seqkit', 'split2', '--by-size', str(FASTA_CHUNKS_NUM), '--out-dir', FASTA_CHUNKS_DIR]
    _ = run_command(seqkit_cmd, fasta_stdin)
    logger.info('Splitting has been done.')
    check_if_entries_are_exist_and_not_empty(FASTA_CHUNKS_DIR)

    return [entry.path for entry in os.scandir(FASTA_CHUNKS_DIR)]


def check_if_entries_are_exist_and_not_empty(directory: str):
    entries = [entry.path for entry in os.scandir(directory) if entry.is_file()]
    if not entries:
        logger.critical(f'{directory} is an empty, exiting...')
        sys.exit(-1)
    if not all([os.path.getsize(entry) != 0 for entry in entries]):
        logger.critical(f'Some of files in {directory} are an empty, exiting...')


def process_data(*files: str, receptor: str, organism: str, is_fastq=False):
    annotations = []
    for file in files:
        if file:
            if is_fastq:
                fastq_stdout = read_gz(file)
                fasta_stdout = convert_fastq_to_fasta(fastq_stdout)
            else:
                fasta_stdout = read_gz(file)
            fasta_chunk_paths = split_by_chunks(fasta_stdout)

            if receptor == 'all':
                annotations += generate_annotation(fasta_chunk_paths, 'TCR', organism)
                annotations += generate_annotation(fasta_chunk_paths, 'Ig', organism)
            else:
                annotations += generate_annotation(fasta_chunk_paths, receptor, organism)

    return annotations


def check_args(args: argparse.Namespace):
    if not (args.in_fq1 or args.in_fq2 or args.in_fq12 or args.in_fasta):
        logger.critical('One of the arguments --in-fq1, --in-fq2, --in-fq12, --in-fasta is required.')
        sys.exit(1)
    if args.receptor not in RECEPTOR_GLOSSARY:
        logger.critical(f"Invalid receptor name: {args.receptor}. Should be: 'BCR', 'TCR' or 'all', exiting...")
        sys.exit(1)
    if args.organism not in ORGANISM_GLOSSARY:
        logger.critical(f"Invalid organism name: {args.organism}. Should be: 'human' or 'mouse', exiting...")
        sys.exit(1)


def run(args: argparse.Namespace):
    os.chdir(IGBLAST_DIR)
    logger.info(f'Moved inside: {IGBLAST_DIR}')

    decompress(args.in_ref)

    receptor = RECEPTOR_GLOSSARY.get(args.receptor)
    organism = ORGANISM_GLOSSARY.get(args.organism)

    fasta_annotations = process_data(args.in_fasta, receptor=receptor, organism=organism)
    fastq_annotations = process_data(args.in_fq1, args.in_fq2, args.in_fq12,
                                     receptor=receptor, organism=organism, is_fastq=True)

    annotations = fasta_annotations + fastq_annotations
    save_annotation(args.out_annotation, annotations)


if __name__ == "__main__":
    args = parse_args()
    check_args(args)

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
