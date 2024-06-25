from pytest import fixture
import pandas as pd

from filter import (remove_non_canonical, remove_non_functional, remove_non_productive, filter_pgen,
                    _get_duplicates_in_different_loci, drop_duplicates_in_different_loci,
                    _remove_chimeras_by_segment, remove_chimeras, remove_no_junction, discard_junctions_with_n)

from logger import set_logger

logger = set_logger(name=__file__)


@fixture(scope='module')
def annotation_with_duplicates_in_different_loci() -> pd.DataFrame:
    return pd.DataFrame(data={'duplicate_count': [1, 1, 1, 1, 1, 100],
                              'junction': ['AAA', 'AAA', 'AAA', 'AAA', 'AAA', 'AAT'],
                              'locus': ['TRA', 'TRA', 'IGH', 'IGH', 'IGL', 'IGK'],
                              'pgen': [0.0001, 0.0001, 0.6, 0.6, 0.5, 0.9],
                              'j_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002, 0.00001],
                              'v_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002, 0.00001]},
                        index=[0, 1, 0, 1, 2, 0])


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
    return pd.DataFrame(data={'junction_aa': ['AAAA', 'AAAA', 'AA*A', 'AA_A'],
                              'junction': [None, 'AAAA', 'AAAA', 'AAAA']
                              })


@fixture(scope='module')
def annotation_junctions_with_n() -> pd.DataFrame:
    return pd.DataFrame(data={"junction_aa": ["CXA", "KC", "CA", "KC"],
                              "junction": ["TGTNGCG", "AAATGT", "TGTNGCG", "AAATGTN"]
                              })


@fixture(scope='module')
def annotation_non_canonical() -> pd.DataFrame:
    return pd.DataFrame(data={'junction_aa': ['CAAW', 'CAAF', 'CAAA', 'AAAF', 'AAAW'],
                              'j_sequence_alignment_aa': ['FGGGG', 'CWGGG', 'FWGGG', None, 'WCCG']})


@fixture(scope='module')
def annotation_out_of_frame() -> pd.DataFrame:
    return pd.DataFrame(data={'junction': ['TTTTTTT', 'GGGG', 'GGG', 'AAAAAA'], })


@fixture(scope='module')
def annotation_with_v_chimeras() -> pd.DataFrame:
    return pd.DataFrame(data={'locus': ['IGH', 'IGL', 'TRA', 'IGK', 'TRA'],
                              'v_call': ['IGLV3-25*03', 'IGHV2-26*03', 'TRAV4-2*01', 'IGLV2-34*01,IGHV4-34*02', 'TRAV7-6*01,TRDV7-6*02']})


@fixture(scope='module')
def annotation_with_j_chimeras() -> pd.DataFrame:
    return pd.DataFrame(data={'locus': ['IGH', 'IGL', 'TRA', 'IGK', 'TRA'],
                              'j_call': ['IGLJ3-25*03', 'IGHJ2-26*03', 'TRAJ4-2*01', 'IGLJ2-34*01,IGHJ4-34*02', 'TRAJ7-6*01,TRDJ7-6*02']})


@fixture(scope='module')
def empty_annotation() -> pd.DataFrame:
    columns = [
        'sequence', 'locus', 'stop_codon', 'vj_in_frame', 'v_frameshift', 'productive', 'v_call', 'j_call', 'junction',
        'junction_aa', 'v_support', 'j_support', 'v_sequence_start', 'v_sequence_end', 'j_sequence_start',
        'j_sequence_end', 'j_sequence_alignment_aa', 'pgen', 'duplicate_count'
    ]
    return pd.DataFrame(columns=columns)


def test_remove_no_junction(annotation_non_functional):
    filtered_annotation, no_junction_count = remove_no_junction(annotation_non_functional)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'junction_aa': ['AAAA', 'AA*A', 'AA_A'],
                           'junction': ['AAAA', 'AAAA', 'AAAA']
                           },
                     index=[1, 2, 3])
    )
    assert no_junction_count == {"no_junction": 1}


def test_remove_no_junction_on_empty_annotation(empty_annotation):
    filtered_annotation, no_junction_count = remove_no_junction(empty_annotation)
    assert filtered_annotation.empty and no_junction_count == {'no_junction': 0}


def test_remove_non_functional(annotation_non_functional):
    filtered_annotation = remove_non_functional(annotation_non_functional)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'junction_aa': ['AAAA'],
                           'junction': ['AAAA']},
                     index=[1])
    )


def test_remove_non_functional_on_empty_annotation(empty_annotation):
    filtered_annotation = remove_non_functional(empty_annotation)
    assert filtered_annotation.empty


def test_remove_non_productive(annotation_non_productive):
    filtered_annotation = remove_non_productive(annotation_non_productive)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'stop_codon': ['F', 'F'],
                           'vj_in_frame': ['T', 'T'],
                           'v_frameshift': ['F', 'F'],
                           'productive': ['T', 'T'], },
                     index=[0, 2])
    )


def test_remove_non_productive_on_empty_annotation(empty_annotation):
    filtered_annotation = remove_non_productive(empty_annotation)
    assert filtered_annotation.empty


def test_remove_non_canonical(annotation_non_canonical):
    filtered_annotation = remove_non_canonical(annotation_non_canonical)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'junction_aa': ['CAAW', 'CAAF'],
                           'j_sequence_alignment_aa': ['FGGGG', 'CWGGG']},
                     index=[0, 1])
    )


def test_remove_non_canonical_on_empty_annotation(empty_annotation):
    filtered_annotation = remove_non_canonical(empty_annotation)
    assert filtered_annotation.empty


def test_get_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci):
    assert _get_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci)[0].equals(
        pd.DataFrame(data={'duplicate_count': [1, 1, 1, 1, 1],
                           'junction': ['AAA', 'AAA', 'AAA', 'AAA', 'AAA'],
                           'locus': ['TRA', 'TRA', 'IGH', 'IGH', 'IGL'],
                           'pgen': [0.0001, 0.0001, 0.6, 0.6, 0.5],
                           'j_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002],
                           'v_support': [0.12, 0.15, 0.0001, 0.0001, 0.0002]},
                     index=[0, 1, 0, 1, 2])
    )


def test_drop_duplicates_in_different_loci_without_pgen(annotation_with_duplicates_in_different_loci):
    filtered_annotation = drop_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [1, 1, 100],
                           'junction': ['AAA', 'AAA', 'AAT'],
                           'locus': ['IGH', 'IGH', 'IGK'],
                           'pgen': [0.6, 0.6, 0.9],
                           'j_support': [0.0001, 0.0001, 0.00001],
                           'v_support': [0.0001, 0.0001, 0.00001]},
                     index=[2, 3, 5])
    )


def test_drop_duplicates_in_different_loci_with_pgen(annotation_with_duplicates_in_different_loci):
    filtered_annotation = drop_duplicates_in_different_loci(annotation_with_duplicates_in_different_loci)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [1, 1, 100],
                           'junction': ['AAA', 'AAA', 'AAT'],
                           'locus': ['IGH', 'IGH', 'IGK'],
                           'pgen': [0.6, 0.6, 0.9],
                           'j_support': [0.0001, 0.0001, 0.00001],
                           'v_support': [0.0001, 0.0001, 0.00001]},
                     index=[2, 3, 5])
    )


def test_remove_v_chimeras(annotation_with_v_chimeras):
    filtered_annotation = _remove_chimeras_by_segment(annotation_with_v_chimeras, 'v')
    logger.info(filtered_annotation)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'locus': ['TRA', 'TRA'],
                           'v_call': ['TRAV4-2*01', 'TRAV7-6*01,TRDV7-6*02']},
                     index=[2, 4])
    )


def test_remove_j_chimeras(annotation_with_j_chimeras):
    filtered_annotation = _remove_chimeras_by_segment(annotation_with_j_chimeras, 'j')
    logger.info(filtered_annotation)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'locus': ['TRA', 'TRA'],
                           'j_call': ['TRAJ4-2*01', 'TRAJ7-6*01,TRDJ7-6*02']},
                     index=[2, 4])
    )


def test_remove_chimeras_on_empty_annotation(empty_annotation):
    filtered_annotation = remove_chimeras(empty_annotation)
    assert filtered_annotation.empty


def test_filter_pgen_singletons(annotation_pgen):
    filtered_annotation = filter_pgen(annotation_pgen, pgen_threshold=0, filter_pgen_singletons=True)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [2, 100, 2, 100, 2],
                           'pgen': [None, 0, 0, 0.9, 0.9]},
                     index=[2, 3, 4, 5, 6])
    )


def test_filter_pgen_singletons_on_empty_annotation(empty_annotation):
    filtered_annotation = filter_pgen(empty_annotation, pgen_threshold=0, filter_pgen_singletons=True)
    assert filtered_annotation.empty


def test_filter_pgen_default(annotation_pgen):
    filtered_annotation = filter_pgen(annotation_pgen, pgen_threshold=0, filter_pgen_singletons=False)
    assert filtered_annotation.equals(
        pd.DataFrame(data={'duplicate_count': [2, 100, 2],
                           'pgen': [None, 0.9, 0.9]},
                     index=[2, 5, 6])
    )


def test_filter_pgen_default_on_empty_annotation(empty_annotation):
    filtered_annotation = filter_pgen(empty_annotation, pgen_threshold=0, filter_pgen_singletons=False)
    assert filtered_annotation.empty


def test_discard_junctions_with_n(annotation_junctions_with_n):
    filtered_annotation = discard_junctions_with_n(annotation_junctions_with_n)
    assert filtered_annotation.equals(
        pd.DataFrame(data={"junction_aa": ["KC"],
                           "junction": ["AAATGT"]},
                     index=[1])
    )


def test_discard_junctions_with_n_on_empty_annotation(empty_annotation):
    filtered_annotation = discard_junctions_with_n(empty_annotation)
    assert filtered_annotation.empty