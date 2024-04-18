from pytest import fixture
import pandas as pd
import os

import airr
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
def annotation_with_chimeras() -> pd.DataFrame:
    return pd.DataFrame(data={'locus': ['IGH', 'IGL'],
                              'v_call': ['IGHV', 'IGLV'],
                              'j_call': ['IGLJ', 'IGHJ']})


def test_remove_non_functional(annotation_non_functional):
    filtered_annotation = filter.remove_non_functional(annotation_non_functional)
    assert len(filtered_annotation[filtered_annotation['junction'].isna() |
                                   filtered_annotation['junction_aa'].str.contains('_', na=False) |
                                   filtered_annotation['junction_aa'].str.contains('\*', na=False) |
                                   ~filtered_annotation['junction_aa'].str.startswith('C', na=False) |
                                   ~(filtered_annotation['junction_aa'].str.endswith('F', na=False) |
                                     filtered_annotation['junction_aa'].str.endswith('W', na=False))]) == 0


def test_remove_non_productive(annotation_non_productive):
    filtered_annotation = filter.remove_non_productive(annotation_non_productive)
    assert len(filtered_annotation[(filtered_annotation['stop_codon'] == 'T') |
                                   (filtered_annotation['vj_in_frame'] == 'F') |
                                   (filtered_annotation['v_frameshift'] == 'T') |
                                   (filtered_annotation['productive'] == 'F')]) == 0


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


def test_remove_v_chimeras(annotation_with_chimeras):
    filtered_annotation = filter.remove_chimeras(annotation_with_chimeras)
    v_chimeras_count = filtered_annotation.T.apply(lambda x: chimeras_check(x, segment='v')).sum()
    assert v_chimeras_count.empty


def test_remove_j_chimeras(annotation_with_chimeras):
    filtered_annotation = filter.remove_chimeras(annotation_with_chimeras)
    j_chimeras_count = filtered_annotation.T.apply(lambda x: chimeras_check(x, segment='j')).sum()
    assert j_chimeras_count.empty


def chimeras_check(clone, segment):
    return any([locus[:3].upper() != clone.locus and locus[:3].upper() not in filter.ALLOWED_LOCUS_CHIMERAS
                for locus in clone[f'{segment}_call'].split(',')])
