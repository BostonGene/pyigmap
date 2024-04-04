import argparse
import glob
import gzip
import json
import os
import subprocess
import shutil
import sys

from logger import set_logger

logger = set_logger(name=__file__)

TMP_DIR = '/tmp'

REF_DIR = '/reference'
os.mkdir(REF_DIR)

RECEPTOR_GLOSSARY = {'TCR': 'TR',
                     'BCR': 'IG',
                     'all': 'all'}

ORGANISM_GLOSSARY = {'human': 'homo-sapiens',
                     'rat': 'rattus-norvegicus',
                     'mouse': 'mus-musculus'}


def parse_args() -> argparse.Namespace:
    """
    Parses run.py script arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input with the forward fastq')
    parser.add_argument('--in-fq2', help='Input with the reverse fastq')
    parser.add_argument('--in-fq12', help='Input with the merged (forward + reverse) fastq')
    parser.add_argument('--out-fasta', help='Output fasta file with detected V(D)J segments')
    parser.add_argument('--out-annotation', help='Output annotation file')
    parser.add_argument('--ref', help='Reference vidjil germline sequences', required=True)
    parser.add_argument('--mode',
                        help='Execution mode: "detect" - returns .fasta file with detected V(D)J segments, '
                             '"annotate" - returns .tsv table with annotated V(D)J segments '
                             'or "all" - returns .fasta and .tsv', required=True)
    parser.add_argument('--receptor', help="Receptor type: 'BCR', 'TCR' or 'all'", required=True)
    parser.add_argument('--organism', help="Organism name: 'human', 'rat' or 'mouse'", default='human')
    parser.add_argument('--logs', help='Output logs file', required=True)

    return parser.parse_args()


def read_gz(file_path: str) -> str:
    """
    Reads .gz archive contents into stdout.

    :param file_path: path to the .gz file

    :return: decompressed file stdout
    """
    logger.info(f"Reading {file_path} file...")
    pigz_cmd = ['pigz', '-dc', file_path]
    pigz_stdout = run_command(pigz_cmd).stdout
    logger.info(f"{file_path} file has been successfully read.")
    return pigz_stdout


def parallels(command: list, stdin: str) -> subprocess.CompletedProcess[str]:
    """
    Parallels calculation

    :param stdin: stdin stream
    :param command: command, which we need to parallelize

    :return: vidjil logs stdout
    """
    parallel_cmd = ['parallel', '-j', str(os.cpu_count()), '--pipe', '-L', '4', '--round-robin'] + command
    parallel_stdout = run_command(parallel_cmd, stdin)
    return parallel_stdout


def run_command(command: list[str], stdin=None) -> subprocess.CompletedProcess[str]:
    logger.info(f'Running {command}...')
    try:
        command_process = subprocess.run(command, text=True, capture_output=True, input=stdin, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(f"Failed to run '{' '.join(command)}', exiting...")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"Undefined error: {e}, exiting...")
        sys.exit(1)

    return command_process


def decompress(archive: str) -> None:
    """
    Decompresses archive into selected directory.

    :param archive: path to the zip archive
    """
    logger.info(f'Decompressing {archive} into {REF_DIR}...')
    tar_cmd = ['tar', '-xvf', archive, '-C', REF_DIR]
    _ = run_command(tar_cmd)
    logger.info(f'Reference {archive} has been successfully decompressed.')


def detect_vdj(fastq_stdin: str, output_basename: str, mode: str, organism: str, receptor_loci: list[str]):
    """
    Detects V(D)J segments

    :param fastq_stdin: fastq standard input
    :param mode: "fasta", "tsv" or "all"
    :param output_basename: output basename
    :param organism: organism name: 'homo-sapiens', 'rattus-norvegicus' or 'mus-musculus'
    :param receptor_loci: receptor loci: ['TRA', 'TRB', 'TRB+', ...]

    return: vidjil standard output (logs)
    """
    logger.info(f'Detecting V(D)J...')
    vidjil_organism_loci_preset = f'{organism}.g:' + ','.join(receptor_loci) if receptor_loci else f'{organism}.g' # example: homo-sapiens.g:IGH,IGK,IGL
    vidjil_cmd = ['vidjil-algo',
                  '--germline', os.path.join(REF_DIR, vidjil_organism_loci_preset),
                  '--out-detected',
                  '--dir', TMP_DIR,
                  '--gz',
                  '--clean-memory',
                  '--base', output_basename,
                  '-']
    if mode == 'detect':
        vidjil_cmd += ['-c', 'detect']
    else:
        vidjil_cmd += ['-c', 'clones', '--all']

    vidjil_logs = parallels(vidjil_cmd, fastq_stdin).stdout
    logger.info(f'V(D)J segments has been successfully detected.')

    logs_file = os.path.join(TMP_DIR, output_basename + '.log')
    with open(logs_file, 'w') as f:
        f.write(vidjil_logs)


def check_if_exist(file: str) -> None:
    """
    Checks if file exists

    :param file: path to the file
    """
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    if os.path.getsize(file) == 0:
        logger.critical(f'Output file {file} is an empty, exiting....')
        sys.exit(1)
    logger.info(f"Expected file {file} successfully found.")


def replace_file(src_file: str, dst_file: str) -> None:
    """
    Replaces file from source to the destination

    :param src_file: path to the file, that we need to move
    :param dst_file: path to destination file
    """
    check_if_exist(src_file)
    shutil.move(src_file, dst_file)
    logger.info(f'{src_file} replaced to {dst_file}')


def concat_files(output_file_path: str, files_to_concatenate: list[str], is_compressed: bool) -> None:
    """
    Concatenates files into one file

    :param output_file_path: path to the output file
    :param files_to_concatenate: files, that we need to concatenate
    :param is_compressed: if file is compressed (archived) or not
    """
    writer = gzip.open(output_file_path, 'ab') if is_compressed else open(output_file_path, 'a')
    with writer as out_f:
        for i, file in enumerate(files_to_concatenate):
            reader = gzip.open(file, 'rb') if is_compressed else open(file, 'r')
            with reader as f:
                if file.endswith('tsv.gz') and i > 0:
                    f.readline()  # remove duplicated header line in tsv.gz
                out_f.write(f.read())
                logger.info(f'{file} appended to {output_file_path}!')
    check_if_exist(output_file_path)
    logger.info(f'Output {output_file_path} file has been successfully wrote.')


def prepare_output_files(out_file: str, files_to_concatenate: list[str], is_compressed=True):
    """
    Prepares output file(s)

    :param out_file: path to the output file
    :param files_to_concatenate: files, that we need to concatenate into out_file
    :param is_compressed: if file is compressed (archived) or not
    """
    temp_out_file = os.path.join(TMP_DIR, os.path.basename(out_file))
    logger.info(f'Defined {temp_out_file} as a result file.')
    concat_files(temp_out_file, files_to_concatenate, is_compressed)
    replace_file(temp_out_file, out_file)


def save_output_by_mode(mode: str, out_fasta: str, out_annotation: str, out_logs: str):
    """
    Save outputs files by selected mode
    """
    fasta_files = glob.glob(os.path.join(TMP_DIR, '*detected.vdj.fa.gz'))
    tsv_files = glob.glob(os.path.join(TMP_DIR, '*tsv.gz'))
    log_files = glob.glob(os.path.join(TMP_DIR, '*log'))
    if mode == 'detect' or mode == 'all':
        prepare_output_files(out_fasta, fasta_files)
    if mode == 'annotate' or mode == 'all':
        prepare_output_files(out_annotation, tsv_files)
    prepare_output_files(out_logs, log_files, is_compressed=False)


def get_receptor_loci_for_organism(organism: str, receptor: str) -> list[str]:
    if receptor == 'all':
        return []
    vidjil_preset = os.path.join(REF_DIR, organism + '.g')
    with open(vidjil_preset) as f:
        vidjil_config_dict = json.load(f)
    return [locus for locus in vidjil_config_dict['systems'] if locus.upper().startswith(receptor)]


def check_args(args: argparse.Namespace):
    if not (args.out_fasta or args.out_annotation):
        logger.critical('One of the arguments --out-fasta, --out-annotation is required.')
        sys.exit(1)
    if args.receptor not in RECEPTOR_GLOSSARY:
        logger.critical(f"Invalid receptor name: {args.receptor}. Supported receptors: 'BCR', 'TCR' or 'all', exiting...")
        sys.exit(1)
    if args.organism not in ORGANISM_GLOSSARY:
        logger.critical(f"Invalid organism name: {args.organism}. Supported organisms: 'human', 'mouse' or 'rat', exiting...")
        sys.exit(1)


def run(args: argparse.Namespace) -> None:
    decompress(args.ref)

    receptor = RECEPTOR_GLOSSARY.get(args.receptor)
    organism = ORGANISM_GLOSSARY.get(args.organism)
    receptor_loci = get_receptor_loci_for_organism(organism, receptor)

    for fastq_file in [args.in_fq1, args.in_fq2, args.in_fq12]:
        if fastq_file:
            fastq_stdout = read_gz(fastq_file)
            output_basename = os.path.basename(fastq_file) + '.part{#}'
            detect_vdj(fastq_stdout, output_basename, args.mode, organism, receptor_loci)

    save_output_by_mode(args.mode, args.out_fasta, args.out_annotation, args.logs)


if __name__ == "__main__":
    args = parse_args()
    check_args(args)

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
