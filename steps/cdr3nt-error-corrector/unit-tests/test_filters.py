from pytest import fixture
import pandas as pd
import os

import airr
import filter

cwd = os.getcwd()
test_annotation_tcr_paths = os.path.join(cwd, 'unit-tests', 'test_data', 'test_annotation_tcr.tsv.gz')
test_annotation_bcr_paths = os.path.join(cwd, 'unit-tests', 'test_data', 'test_annotation_bcr.tsv.gz')

@fixture(scope='module')
def annotation_object() -> pd.DataFrame:
    return airr.read_annotation(test_annotation_tcr_paths, test_annotation_bcr_paths,
                                only_functional=False,
                                remove_chimeras=False
                                )[0]


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


def test_drop_duplicates_in_different_loci(annotation_object):
    filtered_annotation = filter.drop_duplicates_in_different_loci(annotation_object)
    assert len(filtered_annotation[filtered_annotation['junction'].duplicated(keep=False)]) == 0


def test_remove_out_of_frame(annotation_object):
    filtered_annotation = filter.remove_out_of_frame(annotation_object)
    assert len(filtered_annotation[filtered_annotation['junction'].str.len() % 3 != 0]) == 0
