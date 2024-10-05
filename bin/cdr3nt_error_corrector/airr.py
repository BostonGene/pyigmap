import re
from collections import namedtuple
import pandas as pd

from filter import remove_clones_without_junction, discard_clones_with_n_in_junction, remove_chimeras_clones, \
    remove_non_canonical_clones, remove_non_functional_clones, drop_clones_with_duplicates_in_different_loci
from logger import set_logger

logger = set_logger(name=__file__)

Cdr3Markup = namedtuple('Cdr3Markup', ['junction', 'cdr3_sequence_start', 'cdr3_sequence_end'])
CYS_CODON = re.compile('TG[TC]')
STOP_CODON = re.compile('T(?:AA|AG|GA)')
FGXG_CODON = re.compile('T(?:GG|TT|TC)GG....GG')
FGXG_SHORT_CODON = re.compile('T(?:GG|TT|TC)GG')
CODONS = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
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
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
}


def split_by_loci(annotation: pd.DataFrame) -> tuple[list[pd.DataFrame], list[str]]:
    annotations_by_loci = []
    loci_list = []
    for locus in annotation['locus'].unique().tolist():
        annotations_by_loci.append(annotation[annotation['locus'] == locus])
        loci_list.append(locus)
    return annotations_by_loci, loci_list


def get_loci_count(annotation: pd.DataFrame, suffix: str = '_aligned_reads') -> dict:
    loci_dict = annotation.groupby(['locus']).size().to_dict()
    return {locus + suffix: count for locus, count in loci_dict.items()}


def get_no_call_count(annotation: pd.DataFrame) -> dict[str, int]:
    """Returns the number of uncalled V, D, J and C genes"""
    return {f'no_{column}': int(annotation[column].isna().sum()) for column in ['v_call', 'j_call', 'd_call', 'c_call']}


def concat_annotations(*annotation_paths: str) -> pd.DataFrame:
    annotations = []
    for annotation_path in annotation_paths:
        if annotation_path:
            annotation = pd.read_csv(annotation_path, sep='\t', low_memory=False)
            annotations.append(annotation)
    concatenated_annotation = pd.concat(annotations)
    concatenated_annotation = concatenated_annotation.loc[
        :, ~concatenated_annotation.columns.str.startswith('Unnamed:')
    ]

    return concatenated_annotation


def read_annotation(*annotation_paths: str, only_functional: bool, only_canonical: bool, remove_chimeras: bool,
                    only_best_alignment: bool, discard_junctions_with_n: bool) -> tuple[pd.DataFrame, dict]:
    logger.info('Reading annotation...')
    metrics_dict = {}
    annotation = concat_annotations(*annotation_paths)

    if not len(annotation):
        logger.warning('Annotation is an empty.')
        return annotation, {}

    no_call_count = get_no_call_count(annotation)
    metrics_dict.update(no_call_count)

    annotation = prepare_vdjc_genes_columns(annotation, only_best_alignment)

    annotation, no_junction_count = remove_clones_without_junction(annotation, "junction")
    metrics_dict.update(no_junction_count)

    annotation, no_junction_aa_count = remove_clones_without_junction(annotation, "junction_aa")
    metrics_dict.update(no_junction_aa_count)

    if discard_junctions_with_n:
        annotation = discard_clones_with_n_in_junction(annotation)

    annotation = prepare_duplicate_count_column(annotation)

    annotation = remove_chimeras_clones(annotation) if remove_chimeras else annotation
    loci_count = get_loci_count(annotation)
    metrics_dict.update(loci_count)
    annotation = remove_non_canonical_clones(annotation) if only_canonical else annotation
    annotation = remove_non_functional_clones(annotation) if only_functional else annotation
    annotation = drop_clones_with_duplicates_in_different_loci(annotation)

    logger.info('Annotation has been read.')

    return annotation, metrics_dict


def prepare_duplicate_count_column(annotation: pd.DataFrame) -> pd.DataFrame:
    if 'duplicate_count' not in annotation.columns:
        annotation['duplicate_count'] = 1
    duplicate_count_column = annotation.pop("duplicate_count")
    annotation.insert(0, duplicate_count_column.name, duplicate_count_column)
    annotation.sort_values('duplicate_count', inplace=True, ascending=False)
    return annotation


def correct_c_call(c_call: str, j_call: str) -> str:
    """Corrects TRxC gene names"""
    if j_call.startswith("TR"):
        c_call = j_call[:3] + "C"
        if c_call == "TRBC":
            c_call = c_call + j_call[4]
    return c_call


def prepare_vdjc_genes_columns(annotation: pd.DataFrame, only_best_alignment: bool) -> pd.DataFrame:
    """Filter and prepare V, D, J and C genes columns according to AIRR standards objects"""
    annotation.dropna(subset=['v_call', 'j_call'], inplace=True)
    annotation.fillna({'d_call': '', 'c_call': ''}, inplace=True)

    if only_best_alignment:
        for gene in ['v_call', 'j_call', 'd_call', 'c_call']:
            annotation[gene] = annotation[gene].str.split(',').str[0]

    for column in [
        'v_sequence_end', 'j_sequence_start', 'd_sequence_end',
        'd_sequence_start', 'c_sequence_end', 'c_sequence_start'
    ]:
        annotation[column] = annotation[column].fillna(-1).astype(int)

    annotation["c_call"] = [
        correct_c_call(c_call, v_call) for c_call, v_call in
        zip(annotation["c_call"].values, annotation["j_call"].values)
    ]

    return annotation
