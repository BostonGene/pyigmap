import argparse
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from typing import Optional, List

import pandas as pd

import filter
import airr
from correct import ClonotypeCorrector
from pgen import PgenModel
from logger import set_logger

pd.options.mode.copy_on_write = True

logger = set_logger(name=__file__)

OLGA_MODELS_DIR = tempfile.TemporaryDirectory().name
os.makedirs(OLGA_MODELS_DIR, exist_ok=True)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if not args.olga_models and not args.skip_pgen_calculation:
        msg_list += ["Pgen calculation is enabled, but no OLGA model is provided"]
    if args.skip_pgen_calculation and (args.filter_pgen_all is not None or args.filter_pgen_singletons is not None):
        msg_list += ["Pgen calculation is disabled, but '--filter-pgen-all' or '--filter-pgen-singletons' is provided"]
    if args.filter_pgen_all is not None and args.filter_pgen_singletons is not None:
        msg_list += ["Flags '--filter-pgen-all' and '--filter-pgen-singletons' cannot be provided at the same time"]
    return msg_list


def exit_with_error(message: Optional[str]) -> None:
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


def parse_args() -> argparse.Namespace:
    """Parses arguments"""
    parser = argparse.ArgumentParser()

    inputs_group = parser.add_argument_group('input files')
    inputs_group.add_argument(
        '--in-annotation',
        help='Input raw annotation',
        nargs='+',
        action='extend',
        required=True,
        type=str
    )
    inputs_group.add_argument(
        '--in-json',
        help='Input json(s) with total reads',
        nargs='+',
        action='extend',
        type=str
    )
    inputs_group.add_argument(
        '--olga-models',
        help='Archive with OLGA models',
        type=str
    )

    # Arguments for filtering clonotypes
    filters_group = parser.add_argument_group('filters')
    filters_group.add_argument(
        '--only-productive',
        help='Filter out non-productive clonotypes',
        action='store_true'
    )
    filters_group.add_argument(
        '--remove-chimeras',
        help='Remove chimeras clonotypes, that have different locus in v-/j-genes',
        action='store_true'
    )
    filters_group.add_argument(
        '--only-functional',
        help='Filter out non-functional clonotypes',
        action='store_true'
    )
    filters_group.add_argument(
        '--only-canonical',
        help='Filter out non-canonical clonotypes',
        action='store_true'
    )
    filters_group.add_argument(
        "--only-best-alignment",
        help="Store the best aligned V, D, J and C genes call",
        action="store_true"
    )
    filters_group.add_argument(
        "--discard-junctions-with-n",
        help="Discard clonotypes with undefined nucleotide or amino acid in CDR3 sequence",
        action="store_true"
    )
    filters_group.add_argument(
        '--skip-pgen-calculation',
        help='Skip pgen calculation via OLGA tool',
        action='store_true'
    )
    filters_group.add_argument(
        '--error-rate',
        help='Error rate, that involved in collapsing of clonotype duplicates',
        type=float,
        default=0.001
    )
    filters_group.add_argument(
        "--filter-pgen-all",
        help="All clonotypes with 'pgen <= pgen_threshold' will be removed",
        type=float
    )
    filters_group.add_argument(
        "--filter-pgen-singletons",
        help="All clonotypes with 'duplicate_count == 1 && pgen <= pgen_threshold' will be removed",
        type=float
    )

    # Arguments, that involved in clonotypes aggregation
    aggregation_group = parser.add_argument_group('aggregation')
    aggregation_group.add_argument(
        '--top-c-call',
        help='Returns clonotypes with the most frequent c_call',
        action='store_true'
    )
    aggregation_group.add_argument(
        '--top-v-alignment-call',
        help='Returns clonotypes with the most frequent v alignment',
        action='store_true'
    )

    outputs_group = parser.add_argument_group('outputs')
    outputs_group.add_argument(
        '--out-corrected-annotation',
        help='Output corrected annotation',
        required=True
    )
    outputs_group.add_argument(
        '--out-json',
        help='Output json with metrics',
        required=True
    )
    outputs_group.add_argument(
        '--out-archive',
        help='Output archive with all results',
        required=True
    )

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        exit_with_error("\n".join(error_message_list))

    return args


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
        logger.critical(f'{fail_message} failed with code {e.returncode}')
        print_error_message(e.stderr)
        if exit_on_error:
            logger.critical(f'{exit_on_error=}, now exiting.')
            sys.exit(1)


def unpack_olga_models(archive: str) -> None:
    """Unzips archive into selected directory"""
    logger.info(f'Unzipping {archive} into {OLGA_MODELS_DIR}...')
    tar_cmd = ['tar', '-xvf', archive, '-C', OLGA_MODELS_DIR]
    run_and_check_with_message(tar_cmd, "tar")
    logger.info('OLGA model files have been successfully decompressed.')


def check_if_exist_and_not_empty(file: str) -> None:
    """Checks if file exists"""
    if not os.path.exists(file):
        logger.critical(f"Expected file {file} not found, exiting...")
        sys.exit(1)
    if os.path.getsize(file) == 0:
        logger.critical(f'Output files {file} is an empty, exiting..')
        sys.exit(1)

    logger.info(f"Expected file {file} successfully found and not empty.")


def create_final_archive(*files: str, archive_path: str) -> None:
    """Creates final tar.gz archive"""
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


def save_metrics(*metrics: dict, output_json: str) -> None:
    all_metrics = {}

    for metric in metrics:
        all_metrics.update(metric)

    json_content = json.dumps(all_metrics)

    with open(output_json, 'w') as f:
        f.write(json_content)

    check_if_exist_and_not_empty(output_json)


def save_corrected_annotation(annotation: pd.DataFrame, annotation_path: str) -> None:
    annotation.to_csv(annotation_path, sep='\t', index=False)
    check_if_exist_and_not_empty(annotation_path)


def get_pgen_threshold_value(filter_pgen_singletons: Optional[float],
                             filter_pgen_all: Optional[float]) -> tuple[float, bool]:
    pgen_threshold = filter_pgen_singletons if filter_pgen_singletons is not None else filter_pgen_all
    return pgen_threshold, filter_pgen_singletons is not None


def get_filtered_annotation(annotation: pd.DataFrame, top_c_call: bool, top_v_alignment_call: bool,
                            error_rate: float, skip_pgen_calculation: bool, only_productive: bool,
                            pgen_threshold: float, filter_only_pgen_singletons: bool) -> pd.DataFrame:
    """Filters annotation and calculates Pgen value for each clonotype"""
    corrected_annotations = []
    for annotation_by_locus, locus in zip(*airr.split_by_loci(annotation)):
        logger.info(f'Processing {locus} locus...')

        corrector = ClonotypeCorrector(top_c_call, top_v_alignment_call, error_rate)
        corrected_annotation = corrector.correct_full(annotation_by_locus)

        if not skip_pgen_calculation:
            pgen_model = PgenModel(OLGA_MODELS_DIR, locus)
            corrected_annotation['pgen'] = pgen_model.get_pgen(corrected_annotation['junction_aa'])

        corrected_annotations.append(corrected_annotation)
        logger.info(f'{locus} locus has been processed.')

    concatenated_annotation = pd.concat(corrected_annotations)

    return filter.run_filtration(
        concatenated_annotation, only_productive, pgen_threshold, filter_only_pgen_singletons
    )


def main() -> None:
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {args}")

    if not args.skip_pgen_calculation:
        unpack_olga_models(args.olga_models)

    annotation, metrics_dict = airr.read_annotation(*args.in_annotation,
                                                    only_functional=args.only_functional,
                                                    only_canonical=args.only_canonical,
                                                    remove_chimeras=args.remove_chimeras,
                                                    only_best_alignment=args.only_best_alignment,
                                                    discard_junctions_with_n=args.discard_junctions_with_n)

    pgen_threshold, filter_only_pgen_singletons = get_pgen_threshold_value(args.filter_pgen_singletons,
                                                                           args.filter_pgen_all)

    if len(annotation):
        filtered_annotation = get_filtered_annotation(annotation, args.top_c_call, args.top_v_alignment_call,
                                                      args.error_rate, args.skip_pgen_calculation,
                                                      args.only_productive, pgen_threshold, filter_only_pgen_singletons)
    else:
        filtered_annotation = annotation

    save_corrected_annotation(filtered_annotation, args.out_corrected_annotation)

    total_reads_count = parse_total_reads(args.in_json)
    metrics_dict.update(total_reads_count)

    save_metrics(metrics_dict, output_json=args.out_json)

    create_final_archive(args.out_json, args.out_corrected_annotation, archive_path=args.out_archive)

    logger.info("Run is completed successfully.")


if __name__ == "__main__":
    main()
