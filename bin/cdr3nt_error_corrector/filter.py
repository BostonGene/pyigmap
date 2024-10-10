import pandas as pd

from logger import set_logger

logger = set_logger(name=__file__)

ALLOWED_LOCUS_CHIMERAS = {'TRA', 'TRD'}

STOP_CODON_SIGN = "\\*"  # this sign in AIRR notation means stop codon
FRAME_SHIFT_SIGN = "_"  # this sign in AIRR notation means frame shift

# Regex pattern of immunoreceptor tyrosine-based inhibition motif
# For more information: https://www.pnas.org/doi/10.1073/pnas.121101598
IMMUNORECEPTOR_MOTIF = r"[FW]G.G"

UNKNOWN_NUCLEOTIDE = "N"
UNKNOWN_AMINO_ACID = "X"

CYS = "C"  # CYSTEINE amino acid
TRP = "W"  # TRYPTOPHAN amino acid
PHE = "F"  # PHENYLALANINE amino acid


def run_filtration(annotation: pd.DataFrame, only_productive: bool, pgen_threshold: float,
                   filter_pgen_singletons: bool) -> pd.DataFrame:
    if pgen_threshold is not None:
        annotation = filter_clones_by_pgen(annotation, pgen_threshold, filter_pgen_singletons)
    if only_productive:
        annotation = remove_non_productive_clones(annotation)
    return annotation


def filter_duplicates_by_vj_score(annotation: pd.DataFrame) -> pd.DataFrame:
    """Filters out duplicate sequences in a DataFrame
    by selecting the entries with the highest combined V and J gene scores."""
    vj_score_annotation = annotation[['sequence_id', 'locus', 'v_score', 'j_score']].copy()

    vj_score_annotation['v_score'] = vj_score_annotation['v_score'].fillna(0)
    vj_score_annotation['j_score'] = vj_score_annotation['j_score'].fillna(0)
    vj_score_annotation['score'] = vj_score_annotation['v_score'] + vj_score_annotation['j_score']

    max_vj_score_annotation = vj_score_annotation \
        .groupby('sequence_id', as_index=False)['score'].max() \
        .merge(vj_score_annotation, on=['sequence_id', 'score'])

    return annotation.merge(max_vj_score_annotation[['sequence_id', 'locus']], on=['sequence_id', 'locus'])


def get_duplicates_in_different_loci(annotation: pd.DataFrame) -> list[pd.DataFrame]:
    """Returns list of cdr3 duplicates in different loci"""
    duplicates_with_different_loci = (annotation[annotation.duplicated(subset=['junction'], keep=False)]
                                      .groupby('junction'))
    return [group for _, group in duplicates_with_different_loci]


def drop_clones_with_duplicates_in_different_loci(annotation: pd.DataFrame) -> pd.DataFrame:
    """Drops cdr3 duplicates in different loci"""
    annotation = annotation.reset_index(drop=True)
    loci_duplicates = get_duplicates_in_different_loci(annotation)

    corrected_loci_duplicates = [
        filter_cdr3_duplicates_by_metrics(loci_duplicate) for loci_duplicate in loci_duplicates
    ]

    annotation_without_duplicates = annotation.drop(pd.concat(loci_duplicates).index) if loci_duplicates else annotation
    corrected_annotation = pd.concat([annotation_without_duplicates] + corrected_loci_duplicates)

    corrected_annotation = corrected_annotation.sort_values(by=['locus', 'duplicate_count'], ascending=[True, False])

    logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} duplicates in different loci.')

    return corrected_annotation


def filter_cdr3_duplicates_by_metrics(annotation: pd.DataFrame) -> pd.DataFrame:
    """Filters cdr3 duplicates by j_support and v_support values"""
    locus = ''
    if not annotation['j_support'].isna().all():
        locus = annotation[annotation['j_support'] == annotation['j_support'].min()]['locus'].iloc[0]
    elif not annotation['v_support'].isna().all():
        locus = annotation[annotation['v_support'] == annotation['v_support'].min()]['locus'].iloc[0]
    if locus:
        return annotation[annotation['locus'] == locus]
    return annotation


def remove_chimeras_by_segment(annotation: pd.DataFrame, segment: str) -> pd.DataFrame:
    """Removes chimeras: sequences that have different locus"""
    if annotation.empty:
        return annotation

    locus_values = annotation['locus'].values
    segment_calls = annotation[f'{segment}_call'].str.split(',').values

    # Example,
    # call = IGHV2-5*01 => call[:3] => IGH
    # locus = IGL
    # call != locus => this is chimera, and we need to filter it out.
    not_chimeric_clonotypes = [
        all(not receptor_chain_segment_allele
            or receptor_chain_segment_allele[:3].upper() == locus
            or receptor_chain_segment_allele[:3].upper() in ALLOWED_LOCUS_CHIMERAS
            for receptor_chain_segment_allele in calls)
        for locus, calls in zip(locus_values, segment_calls)
    ]

    filtered_annotation = annotation[not_chimeric_clonotypes]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_annotation.shape[0]} '
                f'chimeras in {segment.upper()} segment.')

    return filtered_annotation


def remove_clones_without_junction(annotation: pd.DataFrame, junction_column: str) -> tuple[pd.DataFrame, dict]:
    """Removes clonotypes without junction (or junction_aa) region nucleotide sequence"""
    filtered_annotation = annotation[~annotation[junction_column].isna()]
    no_junction_count = annotation.shape[0] - filtered_annotation.shape[0]
    logger.info(f'Filtered out {no_junction_count} clones with {junction_column} == None.')
    return filtered_annotation, {f'no_{junction_column}': no_junction_count}


def remove_chimeras_clones(annotation: pd.DataFrame) -> pd.DataFrame:
    """Removes chimeras in the V, J, and C genes.
    D-genes are skipped due to their short length and lack of relevance for downstream analysis"""
    annotation = remove_chimeras_by_segment(annotation, 'v')
    annotation = remove_chimeras_by_segment(annotation, 'j')
    annotation = remove_chimeras_by_segment(annotation, 'c')
    return annotation


def filter_clones_by_pgen(annotation: pd.DataFrame, pgen_threshold: float,
                          filter_pgen_singletons: bool) -> pd.DataFrame:
    """Filters out clones using pgen threshold"""
    condition = (annotation['pgen'].isna()) | (annotation['pgen'] > pgen_threshold)
    if filter_pgen_singletons:
        condition |= annotation['duplicate_count'] != 1
    filtered_annotation = annotation[condition]
    diff_count = annotation.shape[0] - filtered_annotation.shape[0]
    logger.info(f'Filtered out {diff_count} clones with pgen <= {pgen_threshold}.')
    return filtered_annotation


def remove_non_productive_clones(annotation: pd.DataFrame) -> pd.DataFrame:
    """Filters out non-productive (non-canonical) clones based on AIRR fields"""
    filtered_productive = annotation[
        (annotation['stop_codon'] == 'F')
        & (annotation['vj_in_frame'] == 'T')
        & (annotation['v_frameshift'] == 'F')
        & (annotation['productive'] == 'T')
    ]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_productive.shape[0]} non productive clones.')
    return filtered_productive


def remove_non_canonical_clones(annotation: pd.DataFrame) -> pd.DataFrame:
    """Filters out non-canonical clones with junction sequences:
    * not start with a Cys codon
    * and not end with a Phe/Trp codon"""
    filtered_canonical = annotation[
        annotation['j_sequence_alignment_aa'].str.contains(IMMUNORECEPTOR_MOTIF, na=False, regex=True)
        & annotation['junction_aa'].str.startswith(CYS, na=False)
        & (
            annotation['junction_aa'].str.endswith(PHE, na=False)
            | annotation['junction_aa'].str.endswith(TRP, na=False)
        )
    ]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_canonical.shape[0]} non canonical clones.')
    return filtered_canonical


def remove_non_functional_clones(annotation: pd.DataFrame) -> pd.DataFrame:
    """Filters out non-functional clones: clones with stop codon or with frame shift"""
    filtered_functional = annotation[
        ~annotation['junction'].isna() &
        ~annotation['junction_aa'].str.contains(FRAME_SHIFT_SIGN, na=False)
        & ~annotation['junction_aa'].str.contains(STOP_CODON_SIGN, na=False)
    ]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_functional.shape[0]} non functional clones.')
    return filtered_functional


def discard_clones_with_n_in_junction(annotation: pd.DataFrame) -> pd.DataFrame:
    """Discards clones with undefined nucleotides or amino acids in junctions"""
    filtered_annotation = annotation[
        ~annotation['junction'].astype(str).str.contains(UNKNOWN_NUCLEOTIDE, na=False)
        & ~annotation['junction_aa'].astype(str).str.contains(UNKNOWN_AMINO_ACID, na=False)
    ]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_annotation.shape[0]} clones with undefined nucleotide '
                f'and amino acid in CDR3.')
    return filtered_annotation
