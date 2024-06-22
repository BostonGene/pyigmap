from utils import parse_total_reads


def test_parse_total_reads_from_calib_json(calib_json):
    assert parse_total_reads([calib_json]) == {"total_reads": 10000}


def test_parse_total_reads_from_fastp_json(fastp_json):
    assert parse_total_reads([fastp_json]) == {"total_reads": 10000}
