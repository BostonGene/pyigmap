import argparse
import os
import sys
import json
import subprocess
import pandas as pd
import tarfile

import airr
import filter
from correct import ClonotypeCorrector
from pgen import PgenModel

from logger import set_logger
logger = set_logger(name=__file__)

OLGA_MODELS_DIR = '/tmp/olga_models'
os.makedirs(OLGA_MODELS_DIR, exist_ok=True)

pd.options.mode.copy_on_write = True


def parse_args() -> argparse.Namespace:
    """
    Parses arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-tcr-annotation', help='Input raw TCR annotation', nargs='+',
                        action='extend', required=True, type=str)
    parser.add_argument('--in-bcr-annotation', help='Input raw BCR annotation', nargs='+',
                        action='extend', required=True, type=str)
    parser.add_argument('--in-json', help='Input json(s) with total reads', nargs='+',
                        action='extend', type=str)
    parser.add_argument('--pgen-threshold', type=float, default=0, help='Pgen (generation probability) threshold value')
    parser.add_argument('--clonotype-collapse-factor', type=float, default=0.05,
                        help='Factor value, that involved in collapsing of clonotype duplicates')
    parser.add_argument('--only-productive', help='Filter out non-productive clonotypes', action='store_true')
    parser.add_argument('--remove-chimeras', action='store_true',
                        help='Remove chimeras clonotypes, that have different locus in v-/j-genes')
    parser.add_argument('--only-functional', help='Filter out non-functional clonotypes', action='store_true')
    parser.add_argument('--olga-models', type=str, help='Archive with OLGA models', required=True)
    parser.add_argument('--out-corrected-annotation', type=str, help='Output corrected annotation', required=True)
    parser.add_argument('--out-json', type=str, help='Output json with metrics')
    parser.add_argument('--out-archive', type=str, help='Output archive with all results', required=True)
    return parser.parse_args()


def run_command(command: list[str]) -> str:
    logger.info(f'Running {command}...')
    try:
        command_process = subprocess.run(command, text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        logger.critical(f"Failed to run '{' '.join(command)}', exiting...")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"Undefined error: {e}, exiting...")
        sys.exit(1)

    return command_process.stdout


def decompress(archive: str) -> None:
    """
    Unzips archive into selected directory.

    :param archive: Path to the zip archive
    """
    logger.info(f'Unzipping {archive} into {OLGA_MODELS_DIR}...')
    tar_cmd = ['tar', '-xvf', archive, '-C', OLGA_MODELS_DIR]
    _ = run_command(tar_cmd)
    logger.info('OLGA model files have been successfully decompressed.')


def check_if_exist_and_not_empty(file: str) -> None:
    """
    Checks if file exists

    :param file: path to the file
    """
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    if os.path.getsize(file) == 0:
        logger.critical(f'Output files {file} is an empty, exiting..')
        sys.exit(1)

    logger.info(f"Expected file {file} successfully found and not empty.")


def make_archive(*files: str, archive_path: str):
    """
    Makes tar.gz archive

    :param files: list of file path, which we need to archive
    :param archive_path: path to the archive
    """
    try:
        with tarfile.open(archive_path, 'w:gz') as tar:
            for file in files:
                tar.add(file, arcname=os.path.basename(file))
                logger.info(f'{file} has been added to {archive_path}')
    except Exception as e:
        logger.critical(e)
        logger.critical(f'Failed to create archive: {archive_path}, exiting...')
        sys.exit(1)

    check_if_exist_and_not_empty(archive_path)

    logger.info(f'Archiving has been successfully completed: {archive_path}.')


def parse_total_reads(json_files: list[str]) -> dict:
    total_reads_count = 0
    for json_file in json_files:
        with open(json_file, 'r') as f:
            json_content = json.load(f)
        if not json_content.get('summary', False):
            logger.warning(f'Total reads count could not be found in {json_file}')
        else:
            total_reads_count += int(json_content.get('summary', 0).get('before_filtering', 0).get('total_reads', 0))
    return {"total_reads": total_reads_count}


def save_metrics(*metrics: dict, output_json: str):
    all_metrics = dict()

    for metric in metrics:
        all_metrics.update(metric)

    json_content = json.dumps(all_metrics)

    with open(output_json, 'w') as f:
        f.write(json_content)

    check_if_exist_and_not_empty(output_json)


def save_corrected_annotation(annotation: pd.DataFrame, annotation_path: str):
    annotation.to_csv(annotation_path, sep='\t', index=False)
    check_if_exist_and_not_empty(annotation_path)


def run(args: argparse.Namespace) -> None:
    decompress(args.olga_models)

    annotation, loci_count = airr.read_annotation(*args.in_tcr_annotation, *args.in_bcr_annotation,
                                                  only_functional=args.only_functional,
                                                  remove_chimeras=args.remove_chimeras)

    if len(annotation):
        corrected_annotations = []
        for annotation_by_locus, locus in zip(*airr.split_by_loci(annotation)):
            logger.info(f'Processing {locus} locus...')

            corrector = ClonotypeCorrector(args.clonotype_collapse_factor)
            corrected_annotation = corrector.correct_full(annotation_by_locus)

            pgen_model = PgenModel(OLGA_MODELS_DIR, locus)
            corrected_annotation['pgen'] = pgen_model.get_pgen(corrected_annotation['junction_aa'],
                                                               corrected_annotation['junction'])

            corrected_annotations.append(corrected_annotation)
            logger.info(f'{locus} locus has been processed.')

        concatenated_annotation = pd.concat(corrected_annotations)

        filtered_annotation = filter.run_filtration(concatenated_annotation, args.only_productive, args.pgen_threshold)
    else:
        filtered_annotation = annotation

    save_corrected_annotation(filtered_annotation, args.out_corrected_annotation)

    total_reads_count = parse_total_reads(args.in_json)
    save_metrics(total_reads_count, loci_count, output_json=args.out_json)

    make_archive(args.out_json, args.out_corrected_annotation, archive_path=args.out_archive)


if __name__ == "__main__":
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
