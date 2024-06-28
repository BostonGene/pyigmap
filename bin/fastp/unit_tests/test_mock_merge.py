import gzip
import tempfile

from mock_merge import run_mock_merge_reads


def test_mock_merge_reads(fastq1, fastq2):
    merged_fastq = tempfile.NamedTemporaryFile(prefix=".fastq.gz").name
    run_mock_merge_reads(fastq1, fastq2, inner_distance_size=1, out_fq12=merged_fastq)
    with gzip.open(merged_fastq, "rb") as f:
        reads = [line.strip() for line in f.readlines()]
    assert reads == [b'@1/1 mock_merged_15_15',
                     b'ATGCRYSWKMBDHVNNAAAAAAAAAAAGCAT',
                     b'+',
                     b'FFFFFFFFFFFFFFF#FFFFFFFFFFFFFFF',
                     b'@2/1 mock_merged_20_20',
                     b'GGGGAAAATTTTGGGGCCCCNTTTTAAAATTTTGGGGCCCC',
                     b'+',
                     b'FFFFFFFFFFFFFFFFFFFF#FFFFFFFFFFFFFFFFFFFF',
                     b'@3/1 mock_merged_50_50',
                     b'CCCGAACTCTGCCAGTCTGGAGCCCGGAGCTGAAGTAGGATTAGCCTCCANGGAGCCCGGAGCTGAAGTAGGATTAGCCTCCAATGTGACTTCCAAGTCTG',
                     b'+',
                     b'@C@DFFFFHGHHFJJGBGHIEHEIIEIGE<;FHAHFGGECGGIGIIJGJI#B@@DFEFFGHHFDFHHGHIJIFGGAEGG<<A?4C<<B:??<DH:@GEDA7',
                     b'@4/1 mock_merged_50_50',
                     b'TTGGTGTATATGTTGTAATTGAGATTGCTCGGGGGAATAGGTTATGTGATNGTTGAACATTGTTTGTTGGTGTATATGTTGTAATTGAGATTGCTCGGGGG',
                     b'+',
                     b"+==A+?<A<=CB@AAAB4>ABBBAABABB4A<:A6='8AB>AA>AA=A?A#+++A;D?DFDFHDHGI<B@F>HDB=G@GGGHCA=B7;DF9;B?B8?FG8<",
                     b'@5/1 mock_merged_50_50',
                     b'GTTGGGGCCCTGGCCTTTTCAGCTGCGGATCAGGGTGCCTTTCTTTGTGTNATCAAGAGTACAAATGACTCAGCCCTGGCCTTTTCAGCTGCGGATCAGGG',
                     b'+',
                     b'C@CFFFFFHHHHHJJJJJJJJJJIJJII>@DGIIJ?BCHHGCHHGGHGHC#B@@FFFFFHHHGHJJJJJIJIIIGIJJ@GDHIIBFHGIJJDDHGGG@AE>']