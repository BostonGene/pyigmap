from pytest import fixture
import pandas as pd

from airr import get_no_call_count, _prepare_vj_columns


@fixture(scope="function")
def annotation_with_no_call():
    return pd.DataFrame(
        data={
            "v_call": [None, "IGHV4-59*01", None, "IGHV4-59*01"],
            "j_call": [None, None, "IGHJ4-59*01", "IGHJ4-59*01,IGHJ4-59*02"],
            "v_sequence_end": [None, None, None, None],
            "j_sequence_start": [None, None, None, None]
        },
        index=[0, 1, 2, 3],
    )


def test_get_no_call_count(annotation_with_no_call):
    assert get_no_call_count(annotation_with_no_call) == {"no_v_call": 2, "no_j_call": 2}


def test__prepare_vj_columns_with_best_alignment(annotation_with_no_call):
    assert _prepare_vj_columns(annotation_with_no_call, only_best_alignment=True).equals(
        pd.DataFrame(data={"v_call": ["IGHV4-59*01"],
                           "j_call": ["IGHJ4-59*01"],
                           "v_sequence_end": [-1],
                           "j_sequence_start": [-1]},
                     index=[3])
    )


def test__prepare_vj_columns_without_best_alignment(annotation_with_no_call):
    assert _prepare_vj_columns(annotation_with_no_call, only_best_alignment=False).equals(
        pd.DataFrame(data={"v_call": ["IGHV4-59*01"],
                           "j_call": ["IGHJ4-59*01,IGHJ4-59*02"],
                           "v_sequence_end": [-1],
                           "j_sequence_start": [-1]},
                     index=[3])
    )