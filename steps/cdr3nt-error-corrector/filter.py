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
    duplicated_loci_groups = annotation[annotation['junction'].duplicated(keep=False)].groupby(['junction'])
    loci_duplicates = []

    for _, duplicate_loci_df in duplicated_loci_groups:
        loci_duplicates.append(duplicate_loci_df)
    return loci_duplicates


def drop_duplicates_in_different_loci(annotation: pd.DataFrame, use_pgen=False) -> pd.DataFrame:
    loci_duplicates = _get_duplicates_in_different_loci(annotation)

    corrected_loci_duplicates = []

    for loci_duplicate in loci_duplicates:
        corrected_loci_duplicate = _filter_cdr3_by_metrics(loci_duplicate, use_pgen)
        corrected_loci_duplicates.append(corrected_loci_duplicate)

    annotation_without_duplicates = annotation.drop_duplicates(subset=['junction'], keep=False)
    corrected_annotation = pd.concat([annotation_without_duplicates] + corrected_loci_duplicates)

    corrected_annotation = corrected_annotation.sort_values(by=['locus', 'duplicate_count'],
                                                            ascending=[True, False])

    logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} duplicates in different loci.')

    return corrected_annotation


def _filter_cdr3_by_metrics(annotation: pd.DataFrame, use_pgen: bool) -> pd.DataFrame:
    if use_pgen and not any(annotation['pgen'].isna()):
        return annotation[annotation['pgen'] == annotation['pgen'].max()]
    elif not all(annotation['j_support'].isna()):
        return annotation[annotation['j_support'] == annotation['j_support'].min()].drop_duplicates(subset=['junction'],
                                                                                                    keep='first')
    elif not all(annotation['v_support'].isna()):
        return annotation[annotation['v_support'] == annotation['v_support'].min()].drop_duplicates(subset=['junction'],
                                                                                                    keep='first')
    return annotation.drop_duplicates(subset=['junction'], keep='first')


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
