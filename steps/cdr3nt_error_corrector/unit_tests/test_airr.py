from pytest import fixture
import pandas as pd

from airr import get_no_call_count


@fixture(scope="module")
def annotation_with_no_call():
    return pd.DataFrame(
        data={
            "v_call": [None, "IGHV4-59*01", None, "IGHV4-59*01"],
            "j_call": [None, None, "IGHJ4-59*01", "IGHJ4-59*01"],
        },
        index=[0, 1, 2, 3],
    )


def test_get_no_call_count(annotation_with_no_call):
    assert get_no_call_count(annotation_with_no_call) == {"no_v_call": 2, "no_j_call": 2}
