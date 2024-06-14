import gzip
import tempfile
from pytest import fixture


@fixture(scope="module")
def fastq1():
    fastq_path = tempfile.NamedTemporaryFile(suffix=".gz").name
    with gzip.open(fastq_path, "wb") as f:
        f.write(
            b"@1/1\n"
            b"ATGCRYSWKMBDHVN\n"
            b"+\n"
            b"FFFFFFFFFFFFFFF\n"
            b"@2/1\n"
            b"GGGGAAAATTTTGGGGCCCC\n"
            b"+\n"
            b"FFFFFFFFFFFFFFFFFFFF\n"
            b"@3/1\n"
            b"CCCGAACTCTGCCAGTCTGGAGCCCGGAGCTGAAGTAGGATTAGCCTCCA\n"
            b"+\n"
            b"@C@DFFFFHGHHFJJGBGHIEHEIIEIGE<;FHAHFGGECGGIGIIJGJI\n"
            b"@4/1\n"
            b"TTGGTGTATATGTTGTAATTGAGATTGCTCGGGGGAATAGGTTATGTGAT\n"
            b"+\n"
            b"+==A+?<A<=CB@AAAB4>ABBBAABABB4A<:A6='8AB>AA>AA=A?A\n"
            b"@5/1\n"
            b"GTTGGGGCCCTGGCCTTTTCAGCTGCGGATCAGGGTGCCTTTCTTTGTGT\n"
            b"+\n"
            b"C@CFFFFFHHHHHJJJJJJJJJJIJJII>@DGIIJ?BCHHGCHHGGHGHC\n"
        )
    return fastq_path


@fixture(scope="module")
def fastq2():
    fastq_path = tempfile.NamedTemporaryFile(suffix=".gz").name
    with gzip.open(fastq_path, "wb") as f:
        f.write(
            b"@1/2\n"
            b"ATGCRYSWKMBDHVN\n"
            b"+\n"
            b"FFFFFFFFFFFFFFF\n"
            b"@2/2\n"
            b"GGGGCCCCAAAATTTTAAAA\n"
            b"+\n"
            b"FFFFFFFFFFFFFFFFFFFF\n"
            b"@3/2\n"
            b"CAGACTTGGAAGTCACATTGGAGGCTAATCCTACTTCAGCTCCGGGCTCC\n"
            b"+\n"
            b"B@@DFEFFGHHFDFHHGHIJIFGGAEGG<<A?4C<<B:??<DH:@GEDA7\n"
            b"@4/2\n"
            b"CCCCCGAGCAATCTCAATTACAACATATACACCAACAAACAATGTTCAAC\n"
            b"+\n"
            b"+++A;D?DFDFHDHGI<B@F>HDB=G@GGGHCA=B7;DF9;B?B8?FG8<\n"
            b"@5/2\n"
            b"CCCTGATCCGCAGCTGAAAAGGCCAGGGCTGAGTCATTTGTACTCTTGAT\n"
            b"+\n"
            b"B@@FFFFFHHHGHJJJJJIJIIIGIJJ@GDHIIBFHGIJJDDHGGG@AE>\n"
        )
    return fastq_path