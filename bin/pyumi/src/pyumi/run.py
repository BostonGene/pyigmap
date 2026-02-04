import argparse

from pyumi.logger import set_logger
from pyumi.pattern import get_prepared_pattern_and_umi_len
from pyumi.utils import extract_umi, keep_only_paired_reads, save_metrics, save_results, split_by_chunks

logger = set_logger(name=__file__)


def check_argument_consistency(args: argparse.Namespace) -> list[str]:
    msg_list = []
    if not args.fq1_pattern and not args.fq2_pattern:
        msg_list += ['One of the arguments --fq1-pattern or --fq2-pattern is required.']
    return msg_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input forward FASTQ', required=True)
    parser.add_argument('--in-fq2', help='Input reverse FASTQ', required=True)
    parser.add_argument('--fq1-pattern', help='Barcode pattern of forward FASTQ', type=str)
    parser.add_argument('--fq2-pattern', help='Barcode pattern of reverse FASTQ', type=str)
    parser.add_argument('--find-in-reverse-complement', action='store_true')
    parser.add_argument('--max-error', help='Max error (mismatch) per ten nucleotides', type=int, default=1)
    parser.add_argument('--out-fq1', help='Output forward fastq without PCR duplicates', required=True)
    parser.add_argument('--out-fq2', help='Output reverse fastq without PCR duplicates', required=True)
    parser.add_argument('--out-json', help='Output json with umi metrics', required=True)

    args = parser.parse_args()

    error_message_list = check_argument_consistency(args)
    if error_message_list:
        parser.error('\n'.join(error_message_list))

    return parser.parse_args()


def run(args: argparse.Namespace):
    fq1_pattern, fq1_umi_length = get_prepared_pattern_and_umi_len(args.fq1_pattern, max_error=args.max_error)
    fq2_pattern, fq2_umi_length = get_prepared_pattern_and_umi_len(args.fq2_pattern, max_error=args.max_error)

    fq1_filtered, fq2_filtered = keep_only_paired_reads(args.in_fq1, args.in_fq2)

    fq12_chunks = split_by_chunks(fq1_filtered, fq2_filtered)

    processed_fq1, processed_fq2, total_reads_count = extract_umi(fq12_chunks,
                                                                  fq1_pattern,
                                                                  fq2_pattern,
                                                                  args.find_in_reverse_complement)

    fq1_filtered, fq2_filtered = keep_only_paired_reads(processed_fq1, processed_fq2, clear=True)

    save_metrics({'summary': {'before_filtering': {'total_reads': int(total_reads_count)}}},
                 {'fq1_umi_length': fq1_umi_length},
                 {'fq2_umi_length': fq2_umi_length},
                 output_json=args.out_json)

    save_results(fq1_filtered, fq2_filtered, args.out_fq1, args.out_fq2)


def main():
    args = parse_args()
    logger.info(f'Starting program with the following arguments: {vars(args)}')

    run(args)
    logger.info('Run is completed successfully.')


if __name__ == '__main__':
    main()
