from pytest import fixture
import tempfile

from mock_merge import mock_merge_reads


@fixture(scope="module")
def fastq1():
    fastq_path = tempfile.NamedTemporaryFile().name
    with open(fastq_path, 'w') as f:
        f.write(
            f"@1 1:N:0:CACTAGACCA+ATAAGGCAGT\n"
            "ATGCRYSWKMBDHVN\n"
            "+\n"
            "FFFFFFFFFFFFFFF"
        )
    return fastq_path


@fixture(scope="module")
def fastq2():
    fastq_path = tempfile.NamedTemporaryFile().name
    with open(fastq_path, 'w') as f:
        f.write(
            f"@1 2:N:0:CACTAGACCA+ATAAGGCAGT\n"
            "ATGCRYSWKMBDHVN\n"
            "+\n"
            "FFFFFFFFFFFFFFF"
        )
    return fastq_path


def test_mock_merge_reads(fastq1, fastq2):
    merged_fastq = mock_merge_reads(fastq1, fastq2, insert_size=1)
    with open(merged_fastq, "rb") as f:
        reads = [line.strip() for line in f.readlines()]
    assert reads == [b'@1 1:N:0:CACTAGACCA+ATAAGGCAGT mock_merged_15_15',
                     b'ATGCRYSWKMBDHVNNAAAAAAAAAAAGCAT',
                     b'+',
                     b'FFFFFFFFFFFFFFF#FFFFFFFFFFFFFFF']
