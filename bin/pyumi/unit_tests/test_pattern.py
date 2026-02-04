import pytest
from pytest import fixture
from pyumi.pattern import (
    ValidationError,
    add_brackets_around_barcode,
    add_nucleotide_cost,
    get_prepared_pattern_and_umi_len,
    parse_umi_length,
    replace_barcode_type_to_regex_group,
    validate_pattern,
)


@fixture(scope='module')
def pattern1() -> str:
    return '^UMI:N{12}'


@fixture(scope='module')
def pattern2() -> str:
    return '^N{0:2}TGGTATCAACGCAGAGT(UMI:N{14})'


@fixture(scope='module')
def pattern3() -> str:
    return '^N{0:2}TGGTATCAACGCAGAGT(UMI:NNNTNNNTNNN)'


@fixture(scope='module')
def pattern5() -> str:
    return '^N{0:2}TGGTATCAACGCAGAGT(UMI:TNNNTNNNTNNNT)'


@fixture(scope='module')
def pattern6() -> str:
    return '^N{14}'


@fixture(scope='module')
def pattern7() -> str:
    return '^[ATGCN]{0:2}TGGTATCAACGCAGAGT(?P<UMI>[ATGCN]{14})'


@fixture(scope='module')
def pattern8() -> str:
    return '^?P<UMI>[ATGCN]{12}'


@fixture(scope='module')
def pattern9() -> str:
    return 'UMI{9}'


@fixture(scope='module')
def pattern10() -> str:
    return '(UMI{10})'


@fixture(scope='module')
def pattern11() -> str:
    return '(UMI:N{11})'


@fixture(scope='module')
def pattern12() -> str:
    return '^TGGTATCAACGCAGAGTAC(UMI:N{19})TCTTGGGGG'


@fixture(scope='module')
def pattern13() -> str:
    return 'AAAGACAGTGGTATCAACGCAGAGT(?P<UMI>[ATGCN]{4}T[ATGCN]{4}T[ATGCN]{4}TCTT)'


def test_parse_barcode_length_when_present(pattern2, pattern9, pattern10, pattern11):
    assert parse_umi_length(pattern2) == 14
    assert parse_umi_length(pattern9) == 9
    assert parse_umi_length(pattern10) == 10
    assert parse_umi_length(pattern11) == 11


def test_parse_barcode_length_when_missing(pattern3, pattern5):
    assert parse_umi_length(pattern3) == 11
    assert parse_umi_length(pattern5) == 13


def test_validate_pattern_when_pass(pattern2):
    validate_pattern(pattern2, 14)


def test_validate_pattern_when_fail(pattern6):
    with pytest.raises(ValidationError) as excinfo:
        validate_pattern(pattern6, 0)
    assert str(excinfo.value) == 'UMI placeholder does not found in ^N{14}, exiting...'


def test_add_brackets_around_barcode(pattern1):
    assert add_brackets_around_barcode(pattern1, barcode_type='UMI') == '^(UMI:N{12})'


def test_replace_umi_barcode_to_regex_group(pattern2):
    assert (replace_barcode_type_to_regex_group(pattern2, barcode_type='UMI')
            == '^N{0:2}TGGTATCAACGCAGAGT(?P<UMI>N{14})')


def test_add_nucleotide_cost(pattern7, pattern13):
    assert add_nucleotide_cost(pattern=pattern7) == '^[ATGCN]{0:2}(TGGTATCAACGCAGAGT){s<=4}(?P<UMI>[ATGCN]{14})'
    assert add_nucleotide_cost(pattern=pattern13) == '(AAAGACAGTGGTATCAACGCAGAGT){s<=5}(?P<UMI>[ATGCN]{4}(T){s<=1}[ATGCN]{4}(T){s<=1}[ATGCN]{4}(TCTT){s<=1})'


def test_without_add_nucleotide_cost(pattern8):
    assert add_nucleotide_cost(pattern=pattern8) == pattern8


def test_get_prepared_pattern_and_umi_len(pattern2, pattern12):
    assert (get_prepared_pattern_and_umi_len(pattern=pattern2)
            == ('^[ATGCN]{0:2}(TGGTATCAACGCAGAGT){s<=4}(?P<UMI>[ATGCN]{14})', 14))
    assert (get_prepared_pattern_and_umi_len(pattern=pattern12)
            == ('^(TGGTATCAACGCAGAGTAC){s<=4}(?P<UMI>[ATGCN]{19})(TCTTGGGGG){s<=2}', 19))

    # TODO!
    # assert (BarcodePattern(pattern='^N{13}').get_prepared_pattern()
    #         == "^(?P<UMI>[ATGCN]{13})")
