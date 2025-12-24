import argparse

import pandas as pd
from logger import set_logger

from bin.reporter.utils import (
    parse_UMI_data,
    parse_vidjil_fasta,
    parse_igblast_tsv,
    load_final_outputs
)
from bin.reporter.viz.report import create_report

logger = set_logger(name=__file__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1-pyumi',
                        help='Input fastq.gz file after pyumi step, SE or PE pair 1. Needed if UMIs are specified',
                        required=False)
    parser.add_argument('--in-fq2-pyumi',
                        help='Input fastq.gz file after pyumi step, PE pair 2. Needed if UMIs are specified')

    parser.add_argument('--in-fq1-calib',
                        help='Input fastq.gz file after calib step, SE or PE pair 2. Needed if UMIs are specified',
                        required=False)
    parser.add_argument('--in-fq2-calib',
                        help='Input fastq.gz file after calib step, PE pair 2. Needed if UMIs are specified')

    parser.add_argument('--umi-reverse', action='store_true')

    parser.add_argument('--vidjil-fasta',
                        help='Input file with vidjil algo results. Needed if this step was used')

    parser.add_argument('--igblast-tsv',
                        help='File with IgBlast results',
                        required=True)

    parser.add_argument('--final-stat',
                        help='File with final results (json stat file)',
                        required=True)

    parser.add_argument('--final-clones',
                        help='File with IgBlast results (clone data file)',
                        required=True)

    parser.add_argument('--report-file')

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {args}")

    res_umi, sequences_umi = parse_UMI_data(args)
    vidjil_df = parse_vidjil_fasta(args.vidjil_fasta)
    igblast_df = parse_igblast_tsv(args.igblast_tsv)
    stats, final_clones = load_final_outputs(args.final_stat, args.final_clones)

    logger.info('Started report creation...')
    create_report(vidjil_df=vidjil_df, igblast_df=igblast_df,
                  umi_res=res_umi, umi_sequences=sequences_umi,
                  final_clones_df=final_clones,
                  report_file_name=args.report_file)
    logger.info(f'Report saved successfully in {args.report_file}.')


if __name__ == "__main__":
    main()
