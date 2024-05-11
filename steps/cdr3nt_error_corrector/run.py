import argparse
import os
import pandas as pd

import airr
import filter
from correct import ClonotypeCorrector
from pgen import PgenModel

from utils import OLGA_MODELS_DIR, decompress, save_corrected_annotation, parse_total_reads, save_metrics, make_archive
from logger import set_logger

logger = set_logger(name=__file__)

os.makedirs(OLGA_MODELS_DIR, exist_ok=True)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if not args.olga_models and not args.keep_pgen_calculation:
        msg_list += ["Pgen calculation is enabled, but no OLGA model is provided"]
    if args.keep_pgen_calculation and (args.filter_pgen_all or args.filter_pgen_singletons):
        msg_list += ["Pgen calculation is disabled, but --filter-pgen-all / --filter-pgen-singletons is provided"]
    if args.filter_pgen_all and args.filter_pgen_singletons:
        msg_list += ["Flags--filter-pgen-all and --filter-pgen-singletons cannot be provided at the same time"]
    return msg_list


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
    parser.add_argument('--keep-pgen-calculation', action='store_false', help='Keep pgen calculation via OLGA tool')
    parser.add_argument('--clonotype-collapse-factor', type=float, default=0.05,
                        help='Factor value, that involved in collapsing of clonotype duplicates')
    parser.add_argument('--only-productive', help='Filter out non-productive clonotypes', action='store_true')
    parser.add_argument('--remove-chimeras', action='store_true',
                        help='Remove chimeras clonotypes, that have different locus in v-/j-genes')
    parser.add_argument('--only-functional', help='Filter out non-functional clonotypes', action='store_true')
    parser.add_argument("--filter-pgen-all", type=float,
                        help="All clonotypes with 'pgen <= pgen_threshold' will be removed")
    parser.add_argument("--filter-pgen-singletons", type=float,
                        help="All clonotypes with 'duplicate_cound == 1 && pgen <= pgen_threshold' will be removed")
    parser.add_argument('--olga-models', type=str, help='Archive with OLGA models')
    parser.add_argument('--out-corrected-annotation', type=str, help='Output corrected annotation', required=True)
    parser.add_argument('--out-json', type=str, help='Output json with metrics')
    parser.add_argument('--out-archive', type=str, help='Output archive with all results', required=True)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error("\n".join(error_message_list))

    return args


def run(args: argparse.Namespace) -> None:
    if not args.keep_pgen_calculation:
        decompress(args.olga_models)

    annotation, metrics_dict = airr.read_annotation(*args.in_tcr_annotation, *args.in_bcr_annotation,
                                                    only_functional=args.only_functional,
                                                    remove_chimeras=args.remove_chimeras)

    pgen_threshold = args.filter_pgen_all or args.filter_pgen_singletons

    if len(annotation):
        corrected_annotations = []
        for annotation_by_locus, locus in zip(*airr.split_by_loci(annotation)):
            logger.info(f'Processing {locus} locus...')

            corrector = ClonotypeCorrector(args.clonotype_collapse_factor)
            corrected_annotation = corrector.correct_full(annotation_by_locus)

            if not args.keep_pgen_calculation:
                pgen_model = PgenModel(OLGA_MODELS_DIR, locus)
                corrected_annotation['pgen'] = pgen_model.get_pgen(corrected_annotation['junction_aa'])

            corrected_annotations.append(corrected_annotation)
            logger.info(f'{locus} locus has been processed.')

        concatenated_annotation = pd.concat(corrected_annotations)

        filtered_annotation = filter.run_filtration(concatenated_annotation, args.only_productive,
                                                    pgen_threshold, args.filter_pgen_singletons)
    else:
        filtered_annotation = annotation

    save_corrected_annotation(filtered_annotation, args.out_corrected_annotation)

    total_reads_count = parse_total_reads(args.in_json)
    metrics_dict.update(total_reads_count)

    save_metrics(metrics_dict, output_json=args.out_json)

    make_archive(args.out_json, args.out_corrected_annotation, archive_path=args.out_archive)


if __name__ == "__main__":
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
