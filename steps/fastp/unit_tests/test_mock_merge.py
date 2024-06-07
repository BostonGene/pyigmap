import gzip

from mock_merge import mock_merge_reads


def test_mock_merge_reads(fastq1, fastq2):
    merged_fastq = mock_merge_reads(fastq1, fastq2, insert_size=1)
    with gzip.open(merged_fastq, "rb") as f:
        reads = [line.strip() for line in f.readlines()]
    assert reads == [b'@1 1:N:0:CACTAGACCA+ATAAGGCAGT mock_merged_15_15',
                     b'ATGCRYSWKMBDHVNNAAAAAAAAAAAGCAT',
                     b'+',
                     b'FFFFFFFFFFFFFFF#FFFFFFFFFFFFFFF',
                     b'@2 1:N:0:CACTAGACCA+ATAAGGCAGT mock_merged_20_20',
                     b'GGGGAAAATTTTGGGGCCCCNTTTTAAAATTTTGGGGCCCC',
                     b'+',
                     b'FFFFFFFFFFFFFFFFFFFF#FFFFFFFFFFFFFFFFFFFF']