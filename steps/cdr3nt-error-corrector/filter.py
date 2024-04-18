import pandas as pd
from logger import set_logger

logger = set_logger(name=__file__)

ALLOWED_LOCUS_CHIMERAS = {'TRAV', 'TRDV', 'TRAJ', 'TRDJ'}


def run_filtration(annotation: pd.DataFrame, only_productive: bool, pgen_threshold: float) -> pd.DataFrame:
    annotation = filter_pgen(annotation, pgen_threshold) if pgen_threshold is not None else annotation
    annotation = remove_non_productive(annotation) if only_productive else annotation
    annotation = remove_out_of_frame(annotation)
    annotation = drop_duplicates_in_different_loci(annotation, use_pgen=True)
    return annotation


def _get_duplicates_in_different_loci(annotation: pd.DataFrame) -> list[pd.DataFrame]:
    """Returns list of cdr3 duplicates in different loci"""
    duplicates_with_different_loci = (annotation[annotation.duplicated(subset=['junction'], keep=False)]
                                      .groupby('junction'))
    return [group for _, group in duplicates_with_different_loci]


def drop_duplicates_in_different_loci(annotation: pd.DataFrame, use_pgen=False) -> pd.DataFrame:
    """Drops cdr3 duplicates in different loci"""
    loci_duplicates = _get_duplicates_in_different_loci(annotation)

    corrected_loci_duplicates = [_filter_cdr3_duplicates_by_metrics(loci_duplicate, use_pgen)
                                 for loci_duplicate in loci_duplicates]

    annotation_without_duplicates = annotation.drop(pd.concat(loci_duplicates).index) if loci_duplicates else annotation
    corrected_annotation = pd.concat([annotation_without_duplicates] + corrected_loci_duplicates)

    corrected_annotation = corrected_annotation.sort_values(by=['locus', 'duplicate_count'], ascending=[True, False])

    logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} duplicates in different loci.')

    return corrected_annotation


def _filter_cdr3_duplicates_by_metrics(annotation: pd.DataFrame, use_pgen: bool) -> pd.DataFrame:
    """Filters cdr3 duplicates by pgen, j_support, and v_support values"""
    locus = ''
    if use_pgen and not annotation['pgen'].isna().all():
        locus = annotation.loc[annotation['pgen'].idxmax(), 'locus']
    elif not annotation['j_support'].isna().all():
        locus = annotation.loc[annotation['j_support'].idxmin(), 'locus']
    elif not annotation['v_support'].isna().all():
        locus = annotation.loc[annotation['v_support'].idxmin(), 'locus']
    if locus:
        return annotation[annotation['locus'] == locus]
    return annotation


def _remove_chimeras_by_segment(annotation: pd.DataFrame, segment: str):
    not_chimera_mask = []
    for locus, segment_call in zip(annotation['locus'].values, annotation[f'{segment}_call'].values):
        if all(call[:3].upper() == locus or call[:3].upper() in ALLOWED_LOCUS_CHIMERAS for call in segment_call.split()):
            not_chimera_mask.append(True)
        else:
            not_chimera_mask.append(False)
    filtered_annotation = annotation[not_chimera_mask]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_annotation.shape[0]} chimeras in {segment.upper()} segment.')
    return filtered_annotation


def remove_chimeras(annotation: pd.DataFrame) -> pd.DataFrame:
    annotation = _remove_chimeras_by_segment(annotation, 'v')
    annotation = _remove_chimeras_by_segment(annotation, 'j')
    return annotation


def filter_pgen(annotation: pd.DataFrame, pgen_threshold) -> pd.DataFrame:
    filtered_annotation = annotation[(annotation['pgen'].isna()) | (annotation['pgen'].values > pgen_threshold)]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_annotation.shape[0]} clones with pgen <= {pgen_threshold}.')
    return filtered_annotation


def remove_out_of_frame(annotation: pd.DataFrame) -> pd.DataFrame:
    filtered_annotation = annotation[annotation['junction'].str.len() % 3 == 0]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_annotation.shape[0]} clones out of frame.')
    return filtered_annotation


def remove_non_productive(annotation: pd.DataFrame) -> pd.DataFrame:
    filtered_productive = annotation[(annotation['stop_codon'] == 'F') &
                                     (annotation['vj_in_frame'] == 'T') &
                                     (annotation['v_frameshift'] == 'F') &
                                     (annotation['productive'] == 'T')]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_productive.shape[0]} non productive clones.')
    return filtered_productive


def remove_non_functional(annotation: pd.DataFrame) -> pd.DataFrame:
    filtered_functional = annotation[~annotation['junction'].isna() &
                                     ~annotation['junction_aa'].str.contains('_', na=False) &
                                     ~annotation['junction_aa'].str.contains('\*', na=False) &
                                     annotation['junction_aa'].str.startswith('C', na=False) &
                                     (annotation['junction_aa'].str.endswith('F', na=False) |
                                      annotation['junction_aa'].str.endswith('W', na=False))]
    logger.info(f'Filtered out {annotation.shape[0] - filtered_functional.shape[0]} non functional clones.')
    return filtered_functional
