import argparse

import pandas as pd

from reporter.logger import set_logger
from reporter.utils import get_consensus_group_size_per_read, get_read_to_umi_mapping
from reporter.viz import create_report

logger = set_logger(name=__file__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1-pyumi', help='Input fastq.gz file after pyumi step, SE or PE pair 1', required=True)
    parser.add_argument('--in-fq2-pyumi', help='Input fastq.gz file after pyumi step, PE pair 2')

    parser.add_argument('--in-fq1-calib', help='Input fastq.gz file after calib step, SE or PE pair 2', required=True)
    parser.add_argument('--in-fq2-calib', help='Input fastq.gz file after calib step, PE pair 2')

    parser.add_argument('--umi-reverse', action='store_true')

    parser.add_argument('--report-file')

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logger.info(f'Starting program with the following arguments: {args}')

    pyumi_data_path = args.in_fq1_pyumi if not args.umi_reverse else args.in_fq2_pyumi
    calib_data_path = args.in_fq1_calib if not args.umi_reverse else args.in_fq2_calib

    logger.info(f'Started reading initial data from {pyumi_data_path}')
    umi_to_count_mapping_pre, id_to_umi, sequences = get_read_to_umi_mapping(pyumi_data_path)

    logger.info(f'Started reading calib data from {calib_data_path}')
    umi_to_count_mapping_post = get_consensus_group_size_per_read(calib_data_path, id_to_umi)

    logger.info('Created dataset')
    res = (
        pd.DataFrame({'umi': umi_to_count_mapping_pre.keys(), 'read_num_pre_dedup': umi_to_count_mapping_pre.values()})
        .merge(
            pd.DataFrame(
                {'umi': umi_to_count_mapping_post.keys(), 'read_num_post_dedup': umi_to_count_mapping_post.values()}
            ),
            how='outer',
        )
        .fillna(0)
    )

    res['umi'] = res['umi'].apply(lambda x: x.decode('ascii'))

    logger.info('Started report creation...')
    create_report(res, sequences, report_file_name=args.report_file)
    logger.info(f'Report saved successfully in {args.report_file}.')


if __name__ == '__main__':
    main()
