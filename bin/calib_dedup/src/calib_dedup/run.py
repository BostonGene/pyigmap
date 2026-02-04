import argparse

from calib_dedup.logger import set_logger
from calib_dedup.utils import check_error_tolerance_size, cluster_umi, decompress, generate_consensus, save_results

logger = set_logger(name=__file__)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if args.fq1_umi_length is None and args.fq2_umi_length is None:
        msg_list += ['One of the arguments --fq1-umi-length or --fq2-umi-length must be provided.']
    return msg_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input first FASTQ', required=True)
    parser.add_argument('--in-fq2', help='Input second FASTQ', required=True)
    parser.add_argument('--kmer-size', type=int, default=4)
    parser.add_argument('--minimizer-count', type=int, default=7)
    parser.add_argument('--fq1-umi-length', type=int)
    parser.add_argument('--fq2-umi-length', type=int)
    parser.add_argument(
        '--minimizer-threshold', type=int, default=3, help='Error threshold between reads (without UMI)'
    )
    parser.add_argument('--error-tolerance', help='Hamming distance between UMIs', type=int, default=2)
    parser.add_argument('--min-reads-per-cluster', type=int, default=1)
    parser.add_argument('--max-reads-per-cluster', type=int, default=1000)
    parser.add_argument('--out-fq1', help='Output first deduplicated FASTQ', required=True)
    parser.add_argument('--out-fq2', help='Output second deduplicated FASTQ', required=True)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error('\n'.join(error_message_list))

    return parser.parse_args()


def run(args: argparse.Namespace):
    check_error_tolerance_size(args.fq1_umi_length, args.error_tolerance)
    check_error_tolerance_size(args.fq2_umi_length, args.error_tolerance)

    in_fq1 = decompress(args.in_fq1)
    in_fq2 = decompress(args.in_fq2)

    cluster_file = cluster_umi(
        in_fq1,
        in_fq2,
        args.fq1_umi_length,
        args.fq2_umi_length,
        args.kmer_size,
        args.minimizer_count,
        args.minimizer_threshold,
        args.error_tolerance,
    )

    fq1_cons, fq2_cons = generate_consensus(
        in_fq1, in_fq2, cluster_file, args.min_reads_per_cluster, args.max_reads_per_cluster
    )

    save_results(fq1_cons, fq2_cons, args.out_fq1, args.out_fq2)


def main():
    args = parse_args()
    logger.info(f'Starting program with the following arguments: {vars(args)}')

    run(args)
    logger.info('Run is completed successfully.')


if __name__ == '__main__':
    main()
