from pytest import fixture
import pandas as pd

from airr import get_no_call_count, _prepare_vdjc_genes_columns, read_annotation


@fixture(scope="function")
def annotation_with_no_call():
    return pd.DataFrame(
        data={
            "v_call": [None, "IGHV4-59*01", None, "IGHV4-59*01"],
            "j_call": [None, None, "IGHJ4-59*01", "IGHJ4-59*01,IGHJ4-59*02"],
            "d_call": [None, None, None, "IGHD3-22*01,IGHD3-22*02"],
            "c_call": [None, None, None, "IGHM1,IGHM2"],
            "v_sequence_end": [None, None, None, None],
            "j_sequence_start": [None, None, None, None],
            "d_sequence_end": [None, None, None, None],
            "d_sequence_start": [None, None, None, None],
            "c_sequence_end": [None, None, None, None],
            "c_sequence_start": [None, None, None, None],
        },
        index=[0, 1, 2, 3],
    )


def test_get_no_call_count(annotation_with_no_call):
    assert get_no_call_count(annotation_with_no_call) == {"no_v_call": 2,
                                                          "no_j_call": 2,
                                                          "no_d_call": 3,
                                                          "no_c_call": 3}


def test__prepare_vdjc_genes_columns_with_best_alignment(annotation_with_no_call):
    assert _prepare_vdjc_genes_columns(annotation_with_no_call, only_best_alignment=True).equals(
        pd.DataFrame(data={"v_call": ["IGHV4-59*01"],
                           "j_call": ["IGHJ4-59*01"],
                           "d_call": ["IGHD3-22*01"],
                           "c_call": ["IGHM1"],
                           "v_sequence_end": [-1],
                           "j_sequence_start": [-1],
                           "d_sequence_end": [-1],
                           "d_sequence_start": [-1],
                           "c_sequence_end": [-1],
                           "c_sequence_start": [-1]},
                     index=[3])
    )


def test__prepare_vdjc_genes_columns_without_best_alignment(annotation_with_no_call):
    assert _prepare_vdjc_genes_columns(annotation_with_no_call, only_best_alignment=False).equals(
        pd.DataFrame(data={"v_call": ["IGHV4-59*01"],
                           "j_call": ["IGHJ4-59*01,IGHJ4-59*02"],
                           "d_call": ["IGHD3-22*01,IGHD3-22*02"],
                           "c_call": ["IGHM1,IGHM2"],
                           "v_sequence_end": [-1],
                           "j_sequence_start": [-1],
                           "d_sequence_end": [-1],
                           "d_sequence_start": [-1],
                           "c_sequence_end": [-1],
                           "c_sequence_start": [-1]},
                     index=[3])
    )


def test_read_annotation_on_empty_annotation(empty_annotation_file):
    filtered_annotation, metrics_dict = read_annotation(empty_annotation_file,
                                                        only_functional=True,
                                                        only_canonical=True,
                                                        remove_chimeras=True,
                                                        only_best_alignment=True,
                                                        discard_junctions_with_N=True)
    assert filtered_annotation.empty, not metrics_dict
