import argparse
import logging
import os
import sys

from utils import decompress, process_data, save_annotation, RECEPTOR_GLOSSARY, ORGANISM_GLOSSARY, IGBLAST_DIR
from logger import set_logger, TqdmToLogger

logger = set_logger(name=__file__)
tqdm_out = TqdmToLogger(logger, level=logging.INFO)


def parse_args() -> argparse.Namespace:
    """
    Parses run.py script arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input with the forward fastq')
    parser.add_argument('--in-fq2', help='Input with the reverse fastq')
    parser.add_argument('--in-fq12', help='Input with the merged (forward + reverse) fastq')
    parser.add_argument('--in-fasta', help='Input fasta file with detected V(D)J segments')
    parser.add_argument('--receptor', help="Receptor type: 'BCR' or 'TCR'",
                        choices=["BCR", "TCR", "all"], required=True)
    parser.add_argument('--organism', help="Organism name: 'human' or 'mouse'",
                        choices=["human", "mouse"], default='human')
    parser.add_argument('--in-ref', help='FASTA reference with V(D)J segments', required=True)
    parser.add_argument('--out-annotation', help='Output BCR annotation table', required=True)

    return parser.parse_args()


def check_args(args: argparse.Namespace):
    if not (args.in_fq1 or args.in_fq2 or args.in_fq12 or args.in_fasta):
        logger.critical('One of the arguments --in-fq1, --in-fq2, --in-fq12, --in-fasta is required.')
        sys.exit(1)


def run(args: argparse.Namespace):
    os.chdir(IGBLAST_DIR)
    logger.info(f'Moved inside: {IGBLAST_DIR}')

    decompress(args.in_ref)

    receptor = RECEPTOR_GLOSSARY.get(args.receptor)
    organism = ORGANISM_GLOSSARY.get(args.organism)

    fasta_annotations = process_data(args.in_fasta, receptor=receptor, organism=organism)
    fastq_annotations = process_data(args.in_fq1, args.in_fq2, args.in_fq12,
                                     receptor=receptor, organism=organism, is_fastq=True)

    annotations = fasta_annotations + fastq_annotations
    save_annotation(args.out_annotation, annotations)


if __name__ == "__main__":
    args = parse_args()
    check_args(args)

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
