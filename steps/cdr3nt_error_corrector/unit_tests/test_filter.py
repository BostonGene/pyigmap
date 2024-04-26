from pytest import fixture
import pandas as pd

import filter


@fixture(scope='module')
def annotation_with_duplicates_in_different_loci() -> pd.DataFrame:
    return pd.DataFrame(data={'duplicate_count': [1, 1, 1, 1, 1],
                              'junction': ['AAA', 'AAA', 'AAA', 'AAA', 'AAA'],
                              'locus': ['TRA', 'TRA', 'IGH', 'IGH', 'IGL'],
                              'pgen': [0.0001, 0.0001, 0.6, 0.6, 0.5],
                              'j_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002],
                              'v_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002]},
                        index=[0, 1, 0, 1, 2])


@fixture(scope='module')
def annotation_pgen() -> pd.DataFrame:
    return pd.DataFrame(data={'duplicate_count': [1, 1, 2, 100, 2, 100, 2],
                              'pgen': [0, 0, None, 0, 0, 0.9, 0.9]},
                        index=[0, 1, 2, 3, 4, 5, 6])


@fixture(scope='module')
def annotation_non_productive() -> pd.DataFrame:
    return pd.DataFrame(data={'stop_codon': ['F', 'T', 'F', 'F', 'F'],
                              'vj_in_frame': ['T', 'F', 'T', 'T', 'T'],
                              'v_frameshift': ['F', 'F', 'F', 'T', 'F'],
                              'productive': ['T', 'T', 'T', 'T', 'F'], })


@fixture(scope='module')
def annotation_non_functional() -> pd.DataFrame:
    return pd.DataFrame(data={'junction_aa': ['CAAAAW', 'CAA*W', 'CAAAAAF', 'CAA_W', 'AAAAAAA'],
                              'junction': [None, 'AAAAAA', 'AAAAAAA', 'AAAAA', None]
                              })


@fixture(scope='module')
def annotation_out_of_frame() -> pd.DataFrame:
    return pd.DataFrame(data={'junction': ['TTTTTTT', 'GGGG', 'GGG', 'AAAAAA'], })


@fixture(scope='module')
def annotation_with_v_chimeras() -> pd.DataFrame:
    return pd.DataFrame(data={'locus': ['IGH', 'IGL', 'TRA'],
                              'v_call': ['IGLV', 'IGHV', 'TRAV']})


@fixture(scope='module')
def annotation_with_j_chimeras() -> pd.DataFrame:
    return pd.DataFrame(data={'locus': ['IGH', 'IGL', 'TRA'],
                              'j_call': ['IGLJ', 'IGHJ', 'TRAJ']})


def test_remove_non_functional(annotation_non_functional):
    filtered_annotation = filter.remove_non_functional(annotation_non_functional)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'junction_aa': ['CAAAAAF'],
                           'junction': ['AAAAAAA']},
                     index=[2])
    )


def test_remove_non_productive(annotation_non_productive):
    filtered_annotation = filter.remove_non_productive(annotation_non_productive)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'stop_codon': ['F', 'F'],
                           'vj_in_frame': ['T', 'T'],
                           'v_frameshift': ['F', 'F'],
                           'productive': ['T', 'T'], },
                     index=[0, 2])
    )


def test_get_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci):
    assert filter._get_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci)[0].equals(
        annotation_with_duplicates_in_different_loci
    )


def test_drop_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci):
    filtered_annotation = filter.drop_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [1, 1],
                           'junction': ['AAA', 'AAA'],
                           'locus': ['IGH', 'IGH'],
                           'pgen': [0.6, 0.6],
                           'j_support': [0.0001, 0.0001],
                           'v_support': [0.0001, 0.0001]},
                     index=[0, 1])
    )


def test_remove_out_of_frame(annotation_out_of_frame):
    filtered_annotation = filter.remove_out_of_frame(annotation_out_of_frame)
    assert len(filtered_annotation[filtered_annotation['junction'].str.len() % 3 != 0]) == 0


def test_remove_v_chimeras(annotation_with_v_chimeras):
    filtered_annotation = filter._remove_chimeras_by_segment(annotation_with_v_chimeras, 'v')
    assert filtered_annotation.equals(
        pd.DataFrame(data={'locus': ['TRA'],
                           'v_call': ['TRAV']},
                     index=[2])
    )


def test_remove_j_chimeras(annotation_with_j_chimeras):
    filtered_annotation = filter._remove_chimeras_by_segment(annotation_with_j_chimeras, 'j')
    assert filtered_annotation.equals(
        pd.DataFrame(data={'locus': ['TRA'],
                           'j_call': ['TRAJ']},
                     index=[2])
    )


def test_filter_pgen_singletons(annotation_pgen):
    filtered_annotation = filter.filter_pgen(annotation_pgen, pgen_threshold=0, filter_pgen_singletons=True)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [2, 100, 2, 100, 2],
                           'pgen': [None, 0, 0, 0.9, 0.9]},
                     index=[2, 3, 4, 5, 6])
    )


def test_filter_pgen_default(annotation_pgen):
    filtered_annotation = filter.filter_pgen(annotation_pgen, pgen_threshold=0, filter_pgen_singletons=False)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [2, 100, 2],
                           'pgen': [None, 0.9, 0.9]},
                     index=[2, 5, 6])
    )
