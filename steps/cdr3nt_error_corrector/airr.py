import pandas as pd
import re
from collections import namedtuple

import filter
from logger import set_logger

logger = set_logger(name=__file__)

Cdr3Markup = namedtuple('Cdr3Markup', 'junction cdr3_sequence_start cdr3_sequence_end')
CYS_CODON = re.compile('TG[TC]')
STOP_CODON = re.compile('T(?:AA|AG|GA)')
FGXG_CODON = re.compile('T(?:GG|TT|TC)GG....GG')
FGXG_SHORT_CODON = re.compile('T(?:GG|TT|TC)GG')
CODONS = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
          'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
          'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
          'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
          'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
          'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
          'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
          'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
          'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
          'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
          'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
          'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
          'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
          'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
          'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
          'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}


def split_by_loci(annotation: pd.DataFrame) -> tuple[list[pd.DataFrame], list[str]]:
    annotations_by_loci = []
    loci_list = []
    for locus in annotation['locus'].unique().tolist():
        annotations_by_loci.append(annotation[annotation['locus'] == locus])
        loci_list.append(locus)
    return annotations_by_loci, loci_list


def get_loci_count(annotation: pd.DataFrame, suffix='_aligned_reads'):
    loci_dict = annotation.groupby(['locus']).size().to_dict()
    return {locus + suffix: count for locus, count in loci_dict.items()}


def get_no_call_count(annotation: pd.DataFrame) -> dict[str, int]:
    """Returns count of no call in v and j genes"""
    return {
        "no_v_call": len(annotation[annotation['v_call'].isna()]),
        "no_j_call": len(annotation[annotation['j_call'].isna()]),
    }


def _concat_annotations(*annotation_paths: str) -> pd.DataFrame:
    annotations = []
    for annotation_path in annotation_paths:
        if annotation_path:
            annotation = pd.read_csv(annotation_path, sep='\t', low_memory=False)
            annotations.append(annotation)
    concatenated_annotation = pd.concat(annotations)
    concatenated_annotation = concatenated_annotation.loc[:, ~concatenated_annotation.columns.str.startswith('Unnamed:')]

    return concatenated_annotation


def read_annotation(*annotation_paths: str, only_functional: bool, only_canonical: bool, remove_chimeras: bool,
                    discard_junctions_with_N: bool) -> tuple[pd.DataFrame, dict]:
    logger.info('Reading annotation...')
    metrics_dict = {}
    annotation = _concat_annotations(*annotation_paths)

    if not len(annotation):
        logger.warning('Annotation is an empty.')
        return annotation, {}

    no_call_count = get_no_call_count(annotation)
    metrics_dict.update(no_call_count)

    annotation = _prepare_vj_columns(annotation)

    if "duplicate_count" in annotation.columns:
        # run CDR3 processing on Vidjil's annotation
        annotation = _process_cdr3_sequences(annotation)

    annotation, no_junction_count = filter.remove_no_junction(annotation)
    metrics_dict.update(no_junction_count)

    annotation = filter.discard_junctions_with_n(annotation)

    annotation = _prepare_duplicate_count_column(annotation)

    annotation = filter.remove_chimeras(annotation) if remove_chimeras else annotation
    loci_count = get_loci_count(annotation)
    metrics_dict.update(loci_count)
    annotation = filter.remove_non_canonical(annotation) if only_canonical else annotation
    annotation = filter.remove_non_functional(annotation) if only_functional else annotation
    annotation = filter.drop_duplicates_in_different_loci(annotation)

    logger.info('Annotation has been read.')

    return annotation, metrics_dict


def _prepare_duplicate_count_column(annotation: pd.DataFrame):
    if 'duplicate_count' not in annotation.columns:
        annotation['duplicate_count'] = 1
    duplicate_count_column = annotation.pop("duplicate_count")
    annotation.insert(0, duplicate_count_column.name, duplicate_count_column)
    annotation.sort_values('duplicate_count', inplace=True, ascending=False)
    return annotation


def _prepare_vj_columns(annotation: pd.DataFrame) -> pd.DataFrame:
    annotation.dropna(subset=['v_call', 'j_call'], inplace=True)
    annotation['v_sequence_end'] = annotation['v_sequence_end'].fillna(-1).astype(int)
    annotation['j_sequence_start'] = annotation['j_sequence_start'].fillna(-1).astype(int)
    return annotation


def _process_cdr3_sequences(annotation: pd.DataFrame) -> pd.DataFrame:
    cdr3markup_list = [_find_cdr3nt_simple(seq, v_end, j_start) for seq, v_end, j_start in
                       zip(annotation['sequence'].values,
                           annotation['v_sequence_end'].values,
                           annotation['j_sequence_start'].values)]
    annotation['junction'] = [cdr3markup.junction for cdr3markup in cdr3markup_list]
    annotation['cdr3_sequence_start'] = [cdr3markup.cdr3_sequence_start for cdr3markup in cdr3markup_list]
    annotation['cdr3_sequence_end'] = [cdr3markup.cdr3_sequence_end for cdr3markup in cdr3markup_list]
    annotation['junction_aa'] = [_translate_cdr3(junction) for junction in annotation['junction'].values]
    annotation['cdr3aa'] = [junction_aa[1:-1] if junction_aa else '' for junction_aa in
                            annotation['junction_aa'].values]

    return annotation


def _find_inframe_patterns(sequence: str, pattern: re.Pattern) -> list[int]:
    bad_frame_positions = {codon.start() % 3 for codon in STOP_CODON.finditer(sequence)}
    positions = [codon.start() for codon in pattern.finditer(sequence)]
    return [position for position in positions if position % 3 not in bad_frame_positions]


def _find_cdr3nt_simple(sequence: str, v_seq_end=-1, j_seq_start=-1) -> Cdr3Markup:
    v_seq_end = len(sequence) if v_seq_end < 0 else v_seq_end
    j_seq_start = 1 if j_seq_start <= 0 else j_seq_start
    cys_pos = max(_find_inframe_patterns(sequence[:v_seq_end], CYS_CODON), default=-1)
    phe_pos = max(_find_inframe_patterns(sequence[(j_seq_start-1):], FGXG_CODON), default=-1) + j_seq_start
    if cys_pos < 0 or phe_pos <= cys_pos or phe_pos < j_seq_start:
        return Cdr3Markup('', -1, -1)
    return Cdr3Markup(sequence[cys_pos:(phe_pos+2)], cys_pos+4, phe_pos-1)


def _translate_cdr3(sequence: str):
    sequence_length = len(sequence)
    middle = sequence_length // 2
    shift = sequence_length % 3
    pad = '' if shift == 0 else '_' * (3 - shift)
    sequence = sequence[:middle] + pad + sequence[middle:]
    return ''.join([CODONS.get(sequence[i:(i + 3)], '_') for i in range(0, len(sequence), 3)])
