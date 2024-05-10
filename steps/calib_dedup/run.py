import argparse
import sys
import logging

from utils import (get_prepared_pattern, check_error_tolerance_size, keep_only_paired_reads, split_by_chunks,
                   extract_umi, cluster_umi, generate_consensus, save_total_read_count, save_results)
from logger import set_logger, TqdmToLogger

logger = set_logger(name=__file__)
tqdm_out = TqdmToLogger(logger, level=logging.INFO)

TMP_DIR = '/tmp'
FASTQ_CHUNK_SIZE = 2_000_000  # reads count in one fastq chunk


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input forward fastq(.gz)', required=True)
    parser.add_argument('--in-fq2', help='Input reverse fastq(.gz)', required=True)
    parser.add_argument('--kmer-size', type=int, default=8)
    parser.add_argument('--minimizer-count', type=int, default=7)
    parser.add_argument('--fq1-barcode-pattern', type=str)
    parser.add_argument('--fq2-barcode-pattern', type=str)
    parser.add_argument('--find-umi-in-reverse-complement', action='store_true')
    parser.add_argument('--pattern-max-error-budget', type=int, default=10)
    parser.add_argument('--minimizer-threshold', type=int, default=7,
                        help='Error threshold between reads (without UMI)')
    parser.add_argument('--error-tolerance', help='Hamming distance between UMIs', type=int, default=2)
    parser.add_argument('--min-reads-per-cluster', type=int, default=1)
    parser.add_argument('--max-reads-per-cluster', type=int, default=1000)
    parser.add_argument('--out-fq1', help='Output forward fastq without PCR duplicates', required=True)
    parser.add_argument('--out-fq2', help='Output reverse fastq without PCR duplicates', required=True)
    parser.add_argument('--out-json', help='Output json with some metrics', required=True)
    return parser.parse_args()


def run(args: argparse.Namespace):
    fq1_prepared_pattern, fq1_umi_length = get_prepared_pattern(args.fq1_barcode_pattern,
                                                                args.pattern_max_error_budget)
    fq2_prepared_pattern, fq2_umi_length = get_prepared_pattern(args.fq2_barcode_pattern,
                                                                args.pattern_max_error_budget)

    check_error_tolerance_size(fq1_umi_length, args.error_tolerance)
    check_error_tolerance_size(fq2_umi_length, args.error_tolerance)

    fq1_filtered, fq2_filtered = keep_only_paired_reads(args.in_fq1, args.in_fq2)

    fq12_chunks = split_by_chunks(fq1_filtered, fq2_filtered)

    processed_fq1, processed_fq2, total_reads_count = extract_umi(fq12_chunks,
                                                                  fq1_prepared_pattern,
                                                                  fq2_prepared_pattern,
                                                                  args.find_umi_in_reverse_complement)

    fq1_filtered, fq2_filtered = keep_only_paired_reads(processed_fq1, processed_fq2, clear=True)

    cluster_file = cluster_umi(fq1_filtered, fq2_filtered, fq1_umi_length, fq2_umi_length, args.kmer_size,
                               args.minimizer_count, args.minimizer_threshold, args.error_tolerance)

    fq1_cons, fq2_cons = generate_consensus(fq1_filtered, fq2_filtered, cluster_file,
                                            args.min_reads_per_cluster,
                                            args.max_reads_per_cluster)

    save_total_read_count(total_reads_count, args.out_json)
    save_results(fq1_cons, fq2_cons, args.out_fq1, args.out_fq2)


def check_args(args: argparse.Namespace):
    if not args.fq1_barcode_pattern and not args.fq2_barcode_pattern:
        logger.critical('One of the arguments --fq1-barcode-pattern or --fq2-barcode-pattern is required.')
        sys.exit(1)


if __name__ == '__main__':
    args = parse_args()
    logger.info(f"Starting program with the following arguments: {vars(args)}")

    check_args(args)

    run(args)
    logger.info("Run is completed successfully.")
