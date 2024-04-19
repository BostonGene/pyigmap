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
def pattern5() -> str:
    return "^N{0:2}TGGTATCAACGCAGAGT(UMI:TNNNTNNNTNNNN)"


def test_parse_barcode_length_when_present(pattern2):
    assert len(BarcodePattern(pattern=pattern2)) == 14


def test_parse_barcode_length_when_missing(pattern3, pattern5):
    assert len(BarcodePattern(pattern=pattern3)) == 9
    assert len(BarcodePattern(pattern=pattern5)) == 10


def test_validate_pattern_when_pass(pattern2):
    BarcodePattern(pattern=pattern2)._validate_pattern(pattern=pattern2, umi_len=14)


def test_add_brackets_around_barcode(pattern1):
    pattern_with_brackets = BarcodePattern(pattern=pattern1)._add_brackets_around_barcode(pattern1, barcode_type='UMI')
    assert pattern_with_brackets == "^(UMI:N{12})"


def test_replace_umi_barcode_to_regex_group(pattern2):
    assert (BarcodePattern(pattern=pattern2)._replace_barcode_type_to_regex_group(pattern2, barcode_type='UMI')
            == "^N{0:2}TGGTATCAACGCAGAGT(?P<UMI>N{14})")


def test_replace_nucleotide_patterns(pattern2, pattern3):
    assert (BarcodePattern(pattern=pattern2)._replace_nucleotide_patterns(pattern=pattern2)
            == "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(UMI:[ATGCN]{14})")

    assert (BarcodePattern(pattern=pattern3)._replace_nucleotide_patterns(pattern=pattern3)
            == "^[ATGCN]{0:2}TGGTATCAACGCAGAGT(UMI:[ATGCN][ATGCN][ATGCN]T[ATGCN][ATGCN][ATGCN]T[ATGCN][ATGCN][ATGCN])")


def test_add_nucleotide_cost(pattern2):
    pattern_with_cost = BarcodePattern(pattern=pattern2)._add_nucleotide_cost(pattern=pattern2, nucleotides=NORMAL_NUCLEOTIDES, mismatch_cost=2)
    assert pattern_with_cost == "^N{0:2}(TGGTATCAACGCAGAGT){2s<=2}(UMI:N{14})"


def test_without_add_nucleotide_cost(pattern1):
    pattern_with_cost = BarcodePattern(pattern=pattern1)._add_nucleotide_cost(pattern=pattern1, nucleotides=NORMAL_NUCLEOTIDES, mismatch_cost=2)
    assert pattern_with_cost == pattern1
