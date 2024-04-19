from pytest import fixture
from barcode_pattern import BarcodePattern, NORMAL_NUCLEOTIDES, IUPAC_WILDCARDS


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
def pattern4() -> str:
    return "^N{0:2}tggtatcaacgcagagt(UMI:N{14})"


def test_parse_barcode_length_when_present(pattern2):
    assert BarcodePattern(pattern=pattern2)._parse_barcode_length() == 14


def test_parse_barcode_length_when_missing(pattern3):
    pass  # TODO! "^N{0:2}TGGTATCAACGCAGAGT(UMI:NNNTNNNTNNN)" -> UMI len = 9


def test_validate_pattern_when_pass(pattern2):
    BarcodePattern(pattern=pattern2)._validate_pattern(pattern=pattern2, umi_len=14)


def test_add_brackets_around_barcode(pattern1):
    pattern_with_brackets = BarcodePattern(pattern=pattern1)._add_brackets_around_barcode(pattern1, barcode_type='UMI')
    assert pattern_with_brackets == "^(UMI:N{12})"


def test_replace_umi_barcode_to_regex_group(pattern2):
    replaced_umi_barcode = BarcodePattern(pattern=pattern2)._replace_barcode_type_to_regex_group(pattern2, barcode_type='UMI')
    assert replaced_umi_barcode == "^N{0:2}TGGTATCAACGCAGAGT(?P<UMI>N{14})"


def test_replace_nucleotide_patterns(pattern2):
    replaced_pattern = BarcodePattern(pattern=pattern2)._replace_nucleotide_patterns(pattern=pattern2)
    assert replaced_pattern == "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(UMI:[ATGCN]{14})"


def test_add_nucleotide_cost(pattern2):
    pattern_with_cost = BarcodePattern(pattern=pattern2)._add_nucleotide_cost(pattern=pattern2, nucleotides=NORMAL_NUCLEOTIDES, mismatch_cost=6)
    assert pattern_with_cost == "^N{0:2}(TGGTATCAACGCAGAGT){6s<=10}(UMI:N{14})"


def test_without_add_nucleotide_cost(pattern1):
    pattern_with_cost = BarcodePattern(pattern=pattern1)._add_nucleotide_cost(pattern=pattern1, nucleotides=NORMAL_NUCLEOTIDES, mismatch_cost=6)
    assert pattern_with_cost == pattern1
