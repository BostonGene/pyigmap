import gzip
from collections import Counter, defaultdict


def get_consensus_group_size_per_read(fastq_file: str, id_to_umi: dict[int, bytes]) -> dict[bytes, int]:
    """Reads a FASTQ file chunk and returns a list of reads."""
    umi_to_group_size_per_read: defaultdict[bytes, int] = defaultdict(int)
    num_reads = 0
    with gzip.open(fastq_file, 'rb') as f:
        while header := f.readline().strip().decode('utf-8'):  # read header string
            if not header:
                break
            f.readline()
            f.readline()
            f.readline()

            group_ids = header.split('\t')[1].split(';')[:-1]
            consensus_group_size = len(group_ids)  # @0  145607;265853;563279;

            umi = Counter([id_to_umi[int(x)] for x in group_ids]).most_common(1)[0][0]
            umi_to_group_size_per_read[umi] += consensus_group_size
            num_reads += 1

    return umi_to_group_size_per_read


def get_read_to_umi_mapping(
    fastq_file: str, umi_len: int = 12
) -> tuple[defaultdict[bytes, int], dict[int, bytes], list[bytes]]:
    umi_to_count: defaultdict[bytes, int] = defaultdict(int)
    id_to_umi: dict[int, bytes] = {}
    read_id = 0
    sequences: list[bytes] = []
    with gzip.open(fastq_file) as f:
        while header := f.readline().strip().decode('utf-8'):  # read header string
            if not header:
                break
            sequence = f.readline()
            sequences.append(sequence)
            f.readline()
            f.readline()

            umi = sequence[:umi_len]

            umi_to_count[umi] += 1
            id_to_umi[read_id] = umi
            read_id += 1

    return umi_to_count, id_to_umi, sequences
