import argparse
import gzip
import logging
import os
import tempfile
from typing import Optional, List
import shutil
import subprocess
import sys
from utils import TEMPDIR_NAME

LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger()

IGBLAST_DIR = os.environ.get('IGBLAST_DIR')

FASTA_CHUNKS_DIR = tempfile.TemporaryDirectory().name

RECEPTOR_GLOSSARY = {
    'BCR': ['Ig'],
    'TCR': ['TCR'],
    'all': ['Ig', 'TCR']
}


def configure_logger(logger_format: str = LOGGER_FORMAT) -> logging.Logger:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger


def parse_args() -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fastq', help='Input FASTQ (for amplicon)', nargs='+',
                        action='extend')
    parser.add_argument('--in-fasta', help='Input FASTA (Vidjil output; for RNASeq-bulk)')
    parser.add_argument('--organism', help="Organism name: 'human' or 'mouse'",
                        choices=["human", "mouse"], default='human')
    parser.add_argument('--receptor', help="Receptor type: 'BCR', 'TCR' or 'all'",
                        choices=["BCR", "TCR", "all"], type=str, default="all")
    parser.add_argument('--ref', help='IgBLAST reference archive', required=True)
    parser.add_argument('--reads-chunk-size', help='Count of sequences processed in one run of IgBLAST',
                        type=int, default=50_000)
    parser.add_argument('--out-annotation', help='Output AIRR-formatted annotation table', required=True,
                        type=str)

    args = parser.parse_args()

    if args.in_fasta and args.in_fastq:
        parser.error("--in-fasta and --in-fastq cannot be used at the same time.")

    return args


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


def read_seq_gz_file(file_path: str) -> str:
    """Reads .gz file to stdout"""
    logger.info(f"Reading {file_path} file...")
    pigz_cmd = ['pigz', '-dc', file_path]
    pigz_process = run_and_check_with_message(pigz_cmd, "pigz", return_proc=True,
                                              capture_output=True, stderr=None)
    logger.info(f"{file_path} file is successfully read.")

    return pigz_process.stdout


def convert_fastq_to_fasta(fastq_stdin: str) -> str:
    """Converts FASTQ to FASTA"""
    logger.info("Converting FASTQ -> FASTA...")
    seqtk_cmd = ['seqtk', 'seq', '-a']
    seqtk_process = run_and_check_with_message(seqtk_cmd, "seqtk", return_proc=True,
                                               capture_output=True, stderr=None, input=fastq_stdin)
    logger.info("Converting FASTQ -> FASTA successfully done.")

    return seqtk_process.stdout


def run_igblast(seq_file: str, receptor: str, organism: str) -> str:
    logger.info("Going to run IgBLAST...")

    output_annotation = tempfile.NamedTemporaryFile(suffix=".tsv").name

    igblast_cmd = [
        'bin/igblastn',
        '-query', seq_file,
        '-germline_db_V', os.path.join(IGBLAST_DIR, 'database', f'{organism}.{receptor}.V'),
        '-germline_db_D', os.path.join(IGBLAST_DIR, 'database', f'{organism}.{receptor}.D'),
        '-germline_db_J', os.path.join(IGBLAST_DIR, 'database', f'{organism}.{receptor}.J'),
        '-organism', organism,
        '-auxiliary_data', os.path.join(IGBLAST_DIR, 'optional_file', f'{organism}_gl.aux'),
        '-ig_seqtype', receptor,
        '-show_translation',
        '-outfmt', str(19),
        '-num_threads', str(os.cpu_count()),
        '-out', output_annotation
    ]

    if organism == "human":
        igblast_cmd += ['-c_region_db', os.path.join(IGBLAST_DIR, "database", "ncbi_human_c_genes")]

    run_and_check_with_message(igblast_cmd, "igblast")

    logger.info("IgBLAST calculations successfully finished.")

    return output_annotation


def is_file_empty(file: str) -> bool:
    return os.path.getsize(file) == 0


def check_if_entries_exist_and_not_empty(directory: str) -> None:
    entries = [entry.path for entry in os.scandir(directory) if entry.is_file()]
    if not entries:
        exit_with_error(f'{directory} is an empty, exiting...')
    if any(is_file_empty(entry) for entry in entries):
        exit_with_error(f'Some file in {directory} is an empty, exiting...')


def generate_annotations(fasta_chunks: list[str], receptor_type: str, organism: str) -> list[str]:
    """Generate IgBLAST annotation"""
    annotation_paths = []
    for fasta_chunk in fasta_chunks:
        for receptor in RECEPTOR_GLOSSARY[receptor_type]:
            annotation = run_igblast(fasta_chunk, receptor, organism)
            if not is_file_empty(annotation):
                annotation_paths.append(annotation)

    return annotation_paths


def check_if_exist(file: str) -> None:
    """Checks if file exists"""
    if not os.path.exists(file):
        exit_with_error(f"Expected file {file} not found, exiting...")
    logger.info(f"Expected file {file} found.")


def move_file(src_file: str, dst_file: str) -> None:
    """Replaces file from source to the destination"""
    shutil.move(src_file, dst_file)
    logger.info(f'{src_file} moved to {dst_file}.')


def concat_annotations(annotation_files: list[str]) -> str:
    """Concatenates temporary annotations into one file"""
    out_annotation_path = os.path.join(
        TEMPDIR_NAME,
        os.path.basename(tempfile.NamedTemporaryFile(suffix=".tsv").name)
    )
    with gzip.open(out_annotation_path, 'ab') as out_annotation_obj:
        for i, annotation_file in enumerate(annotation_files):
            with open(annotation_file, 'r') as annotation_obj:
                if i > 0:
                    annotation_obj.readline()  # remove duplicated header line
                out_annotation_obj.write(annotation_obj.read().encode('utf-8'))
                logger.info(f'{annotation_file} appended to {out_annotation_path}')
    check_if_exist(out_annotation_path)
    logger.info(f'{out_annotation_path} concatenation has been done.')

    return out_annotation_path


def save_igblast_result(annotations_paths: list[str], out_annotation_path: str) -> None:
    """Concatenate and save IgBLAST outputs"""
    annotation_path = concat_annotations(annotations_paths)
    move_file(annotation_path, out_annotation_path)


def split_fasta_by_chunks(fasta_stdin: str, reads_chunk_size: int) -> list[str]:
    """Splits fasta by chunks"""
    logger.info(f"Splitting fasta into chunks by {reads_chunk_size} reads...")
    seqkit_cmd = ['seqkit', 'split2', '--by-size', str(reads_chunk_size), '--out-dir', FASTA_CHUNKS_DIR]
    run_and_check_with_message(seqkit_cmd, "seqkit", input=fasta_stdin)
    logger.info('Splitting has been done.')
    check_if_entries_exist_and_not_empty(FASTA_CHUNKS_DIR)
    return [entry.path for entry in os.scandir(FASTA_CHUNKS_DIR)]


def get_seq_chunks(seq_file: str, chunk_size: int, is_fastq: bool = False) -> list[str]:
    seq_stdout = read_seq_gz_file(seq_file)
    fasta_stdout = convert_fastq_to_fasta(seq_stdout) if is_fastq else seq_stdout
    return split_fasta_by_chunks(fasta_stdout, chunk_size)


def unpack_reference(archive: str) -> None:
    """Unpack V(D)J reference into selected directory"""
    logger.info(f'Going to unpack {archive}...')

    os.chdir(IGBLAST_DIR)  # It is necessary to correctly run IgBLAST
    logger.info(f'Moved inside: {IGBLAST_DIR}')

    tar_cmd = ['tar', '-I', 'pigz', '-xvf', archive]
    run_and_check_with_message(tar_cmd, "tar")
    logger.info(f'V(D)J reference {archive} successfully unpacked.')


def main() -> None:
    configure_logger()

    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    unpack_reference(args.ref)

    all_annotation_paths = []
    for seq_file in args.in_fastq or []:
        fasta_chunks_list = get_seq_chunks(seq_file, args.reads_chunk_size, is_fastq=True)
        all_annotation_paths += generate_annotations(fasta_chunks_list, args.receptor, args.organism)

    if args.in_fasta:
        fasta_chunks_list = get_seq_chunks(args.in_fasta, args.reads_chunk_size)
        all_annotation_paths += generate_annotations(fasta_chunks_list, args.receptor, args.organism)

    save_igblast_result(all_annotation_paths, args.out_annotation)

    logger.info("Run is completed successfully.")


if __name__ == "__main__":
    main()
