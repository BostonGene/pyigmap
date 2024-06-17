import argparse
import os
import sys

from utils import (RECEPTOR_GLOSSARY, ORGANISM_GLOSSARY, REF_DIR, decompress, get_receptor_loci_for_organism,
                   read_gz, detect_vdj, save_output_by_mode)

from logger import set_logger

logger = set_logger(name=__file__)

os.mkdir(REF_DIR)


def parse_args() -> argparse.Namespace:
    """
    Parses run.py script arguments.

    :return: argparse.Namespace object with parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fq1', help='Input with the forward fastq')
    parser.add_argument('--in-fq2', help='Input with the reverse fastq')
    parser.add_argument('--in-fq12', help='Input with the merged (forward + reverse) fastq')
    parser.add_argument('--out-fasta', help='Output fasta file with detected V(D)J segments')
    parser.add_argument('--out-annotation', help='Output annotation file')
    parser.add_argument('--ref', help='Reference vidjil germline sequences', required=True)
    parser.add_argument('--mode',
                        help='Execution mode: "detect" - returns .fasta file with detected V(D)J segments, '
                             '"annotate" - returns .tsv table with annotated V(D)J segments '
                             'or "all" - returns .fasta and .tsv', choices=["detect", "annotate", "all"], default="detect")
    parser.add_argument('--receptor', help="Receptor type: 'BCR', 'TCR' or 'all'",
                        choices=["BCR", "TCR", "all"], default="all")
    parser.add_argument('--organism', help="Organism name: 'human', 'rat' or 'mouse'",
                        choices=["human", "rat", "mouse"], default="human")
    parser.add_argument('--logs', help='Output logs file', required=True)

    return parser.parse_args()


def check_args(args: argparse.Namespace):
    if not (args.out_fasta or args.out_annotation):
        logger.critical('One of the arguments --out-fasta, --out-annotation is required.')
        sys.exit(1)


def run(args: argparse.Namespace) -> None:
    decompress(args.ref)

    receptor = RECEPTOR_GLOSSARY.get(args.receptor)
    organism = ORGANISM_GLOSSARY.get(args.organism)
    receptor_loci = get_receptor_loci_for_organism(organism, receptor)

    for fastq_file in [args.in_fq1, args.in_fq2, args.in_fq12]:
        if fastq_file:
            fastq_stdout = read_gz(fastq_file)
            output_basename = os.path.basename(fastq_file) + '.part{#}'
            detect_vdj(fastq_stdout, output_basename, args.mode, organism, receptor_loci)

    save_output_by_mode(args.mode, args.out_fasta, args.out_annotation, args.logs)


if __name__ == "__main__":
    args = parse_args()
    check_args(args)

    logger.info(f"Starting program with the following arguments: {vars(args)}")

    run(args)

    logger.info("Run is completed successfully.")
