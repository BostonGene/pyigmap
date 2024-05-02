from pytest import fixture

from extract import (create_new_read, get_reverse_complement, PatternMarkup, BarcodeMarkup, remove_subseq,
                     replace_umi_to_the_seq_start)


@fixture(scope='module')
def read_header() -> str:
    return "M01691:10:000000000-A7F7L:1:1101:13747:1534 2:N:0:1"


@fixture(scope='module')
def read_seq() -> str:
    return "AAAATTTTGGGGCCCC"


@fixture(scope='module')
def read_quality() -> str:
    return "KK?KK.KAKKKK&KKY"


@fixture(scope='module')
def barcode_seq() -> str:
    return "TTTT"


@fixture(scope='module')
def barcode_quality() -> str:
    return "K.KA"


@fixture(scope='module')
def pattern_markup(barcode_seq, barcode_quality) -> str:
    return PatternMarkup(BarcodeMarkup(barcode_seq, barcode_quality), 4, 8)


def test_create_new_read(read_header, read_seq, read_quality):
    assert (create_new_read(read_header, read_seq, read_quality)
            == (f"@{read_header}\n"
                f"{read_seq}\n"
                "+\n"
                f"{read_quality}\n")
            )


def test_get_reverse_complement(read_seq):
    assert get_reverse_complement(read_seq) == 'GGGGCCCCAAAATTTT'


def test_replace_umi_to_the_seq_start(read_seq, read_quality, pattern_markup):
    assert (replace_umi_to_the_seq_start(read_seq, read_quality, pattern_markup)
            == ("TTTTAAAAGGGGCCCC", "K.KAKK?KKKKK&KKY"))


def test_remove_subseq(read_seq):
    assert remove_subseq(read_seq, 0, 8) == 'GGGGCCCC'
