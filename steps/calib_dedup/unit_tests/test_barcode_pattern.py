import pytest
from pytest import fixture

from barcode.pattern import (add_nucleotide_cost, replace_nucleotide_patterns,
                             replace_barcode_type_to_regex_group, add_brackets_around_barcode,
                             validate_pattern, parse_umi_length, add_nucleotide_cost,
                             NORMAL_NUCLEOTIDES, IUPAC_WILDCARDS, ValidationError,
                             get_prepared_pattern_and_umi_len)


@fixture(scope='module')
def pattern1() -> str:
    return "^UMI:N{12}"


@fixture(scope='module')
def pattern2() -> str:
    return "^N{0:2}TGGTATCAACGCAGAGT(UMI:N{14})"


@fixture(scope='module')
def pattern3() -> str:
    return "^N{0:2}TGGTATCAACGCAGAGT(UMI:NNNTNNNTNNN)"


@fixture(scope='module')
def pattern5() -> str:
    return "^N{0:2}TGGTATCAACGCAGAGT(UMI:TNNNTNNNTNNNN)"


@fixture(scope='module')
def pattern6() -> str:
    return "^N{14}"


@fixture(scope='module')
def pattern7() -> str:
    return "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(?P<UMI>[ATGCN]{14})"


@fixture(scope='module')
def pattern8() -> str:
    return "^?P<UMI>[ATGCN]{12}"


def test_parse_barcode_length_when_present(pattern2):
    assert parse_umi_length(pattern2) == 14


def test_parse_barcode_length_when_missing(pattern3, pattern5):
    assert parse_umi_length(pattern3)== 9
    assert parse_umi_length(pattern5) == 10


def test_validate_pattern_when_pass(pattern2):
    validate_pattern(pattern2, 14)


def test_validate_pattern_when_fail(pattern6):
    with pytest.raises(ValidationError) as excinfo:
        validate_pattern(pattern6, 0)
    assert str(excinfo.value) == "UMI placeholder does not found in ^N{14}, exiting..."


def test_add_brackets_around_barcode(pattern1):
    assert add_brackets_around_barcode(pattern1, barcode_type='UMI') == "^(UMI:N{12})"


def test_replace_umi_barcode_to_regex_group(pattern2):
    assert (replace_barcode_type_to_regex_group(pattern2, barcode_type='UMI')
            == "^N{0:2}TGGTATCAACGCAGAGT(?P<UMI>N{14})")


def test_replace_nucleotide_patterns(pattern2, pattern3):
    assert (replace_nucleotide_patterns(pattern2)
            == "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(UMI:[ATGCN]{14})")

    assert (replace_nucleotide_patterns(pattern3)
            == "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(UMI:[ATGCN][ATGCN][ATGCN]T[ATGCN][ATGCN][ATGCN]T[ATGCN][ATGCN][ATGCN])")


def test_add_nucleotide_cost(pattern7):
    assert add_nucleotide_cost(pattern=pattern7) == "^[ATGCN]{0:2}(TGGTATCAACGCAGAGT){2s<=2}(?P<UMI>[ATGCN]{14})"


def test_without_add_nucleotide_cost(pattern8):
    assert add_nucleotide_cost(pattern=pattern8) == pattern8


def test_get_prepared_pattern_and_umi_len(pattern2):
    assert (get_prepared_pattern_and_umi_len(pattern=pattern2)
            == ("^[ATGCN]{0:2}(TGGTATCAACGCAGAGT){2s<=2}(?P<UMI>[ATGCN]{14})", 14))

    # TODO!
    # assert (BarcodePattern(pattern='^N{13}').get_prepared_pattern()
    #         == "^(?P<UMI>[ATGCN]{13})")
