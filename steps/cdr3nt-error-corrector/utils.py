import sys
import json
import subprocess
import tarfile
import os
import pandas as pd

from logger import set_logger
logger = set_logger(name=__file__)

OLGA_MODELS_DIR = '/tmp/olga_models'

pd.options.mode.copy_on_write = True


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
            if 'fastp_version' in json_content['summary']:
                total_reads_count //= 2
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
