import gzip
from collections import defaultdict, Counter
import pandas as pd
import re
from logger import set_logger
import numpy as np
logger = set_logger(name=__file__)



def get_consensus_group_size_per_read(fastq_file: str, id_to_umi: dict[int, str]) -> dict[str, int]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    umi_to_group_size_per_read = defaultdict(int)
    num_reads = 0
    with gzip.open(fastq_file, 'rb') as f:
        while header := f.readline().strip().decode("utf-8"):  # read header string
            if not header:
                break
            _ = f.readline()
            _ = f.readline()
            _ = f.readline()

            group_ids = header.split('\t')[1].split(';')[:-1]
            consensus_group_size = len(group_ids)

            umi = Counter([id_to_umi[int(x)] for x in group_ids]).most_common(1)[0][0]
            umi_to_group_size_per_read[umi] += consensus_group_size
            num_reads += 1

    return umi_to_group_size_per_read


def get_read_to_umi_mapping(fastq_file: str, umi_len: int = 12) -> dict[str, int]:
    umi_to_count = defaultdict(int)
    id_to_umi = {}
    read_id = 0
    sequences = []
    with gzip.open(fastq_file) as f:
        while header := f.readline().strip().decode("utf-8"):  # read header string
            if not header:
                break
            sequence = f.readline()
            sequences.append(sequence)
            _ = f.readline()
            _ = f.readline()

            umi = sequence[:umi_len]

            umi_to_count[umi] += 1
            id_to_umi[read_id] = umi
            read_id += 1

    return umi_to_count, id_to_umi, sequences

def parse_vidjil_fasta(path: str) -> pd.DataFrame:
    if path is None:
        return None
    logger.info(f"Started reading vidjil data from {path}")
    headers = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                sequence_id = header.split()[0]
                match = re.search(r"(\w+) SEG_", header)
                locus_from_seg = match.group(1) if match else None
                headers.append({
                    "header": header,
                    "sequence_id": sequence_id,
                    "locus": locus_from_seg
                })
    return pd.DataFrame(headers)


def parse_UMI_data(args):
    if args.in_fq1_pyumi is None:
        return None, None
    pyumi_data_path = args.in_fq1_pyumi if not args.umi_reverse else args.in_fq2_pyumi
    calib_data_path = args.in_fq1_calib if not args.umi_reverse else args.in_fq2_calib

    logger.info(f"Started reading initial data from {pyumi_data_path}")
    umi_to_count_mapping_pre, id_to_umi, sequences = get_read_to_umi_mapping(pyumi_data_path)

    logger.info(f"Started reading calib data from {calib_data_path}")
    umi_to_count_mapping_post = get_consensus_group_size_per_read(calib_data_path, id_to_umi)

    logger.info(f"Created dataset for UMI processing")
    res = pd.DataFrame({'umi': umi_to_count_mapping_pre.keys(),
                        'read_num_pre_dedup': umi_to_count_mapping_pre.values()}).merge(
        pd.DataFrame({'umi': umi_to_count_mapping_post.keys(),
                      'read_num_post_dedup': umi_to_count_mapping_post.values()}), how='outer'
    ).fillna(0)
    return res, sequences

def parse_igblast_tsv(path: str) -> pd.DataFrame:
    logger.info(f"Started reading igblast data from {path}")
    return pd.read_csv(path, sep="\t", comment="#", dtype=str)


def load_final_outputs(stat_path, clone_path):
    import json
    logger.info(f"Started reading final pyigmap data from {stat_path}, {clone_path}")

    with open(stat_path) as f:
        stats = json.load(f)
    clone_df = pd.read_csv(clone_path, sep="\t")
    return stats, clone_df

def shannon_entropy(probs, base=np.e):
    probs = np.array(probs)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log(probs)) / np.log(base)
