from pytest import fixture
import pandas as pd
import os

import airr
import filter

test_annotation_tcr_paths = os.path.join('unit_tests', 'test_data', 'test_annotation_tcr.tsv.gz')
test_annotation_bcr_paths = os.path.join('unit_tests', 'test_data', 'test_annotation_bcr.tsv.gz')


@fixture(scope='module')
def annotation_object() -> pd.DataFrame:
    return airr.read_annotation(test_annotation_tcr_paths, test_annotation_bcr_paths,
                                only_functional=False,
                                remove_chimeras=False
                                )[0]


@fixture(scope='module')
def annotation_with_duplicates_in_different_loci() -> pd.DataFrame:
    return pd.DataFrame(data={'duplicate_count': [1, 1, 1, 1, 1],
                              'junction': ['AAA', 'AAA', 'AAA', 'AAA', 'AAA'],
                              'locus': ['TRA', 'TRA', 'IGH', 'IGH', 'IGL'],
                              'pgen': [0.0001, 0.0001, 0.6, 0.6, 0.5],
                              'j_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002],
                              'v_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002]})


def test_remove_non_functional(annotation_object):
    filtered_annotation = filter.remove_non_functional(annotation_object)
    assert len(filtered_annotation[filtered_annotation['junction'].isna() |
                                   filtered_annotation['junction_aa'].str.contains('_', na=False) |
                                   filtered_annotation['junction_aa'].str.contains('\*', na=False) |
                                   ~filtered_annotation['junction_aa'].str.startswith('C', na=False) |
                                   ~(filtered_annotation['junction_aa'].str.endswith('F', na=False) |
                                     filtered_annotation['junction_aa'].str.endswith('W', na=False))]) == 0


def test_remove_non_productive(annotation_object):
    filtered_annotation = filter.remove_non_productive(annotation_object)
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
                     index=[2, 3])
    )


def test_remove_out_of_frame(annotation_object):
    filtered_annotation = filter.remove_out_of_frame(annotation_object)
    assert len(filtered_annotation[filtered_annotation['junction'].str.len() % 3 != 0]) == 0


def test_remove_chimeras(annotation_object):
    filtered_annotation = filter.remove_chimeras(annotation_object)
    v_flag = filtered_annotation.T.apply(lambda x: chimeras_check(x, segment='v')).sum()
    j_flag = filtered_annotation.T.apply(lambda x: chimeras_check(x, segment='j')).sum()

    assert not v_flag or j_flag


def chimeras_check(clone, segment):
    return any([locus[:3].upper() != clone.locus and locus[:3].upper() not in filter.ALLOWED_LOCUS_CHIMERAS
                for locus in clone[f'{segment}_call'].split(',')])
