import gzip
import tempfile
from pytest import fixture


@fixture(scope="module")
def fastq1():
    fastq_path = tempfile.NamedTemporaryFile().name
    with gzip.open(fastq_path, "wb") as f:
        f.write(
            b"@1 1:N:0:CACTAGACCA+ATAAGGCAGT\n"
            b"ATGCRYSWKMBDHVN\n"
            b"+\n"
            b"FFFFFFFFFFFFFFF\n"
            b"@2 1:N:0:CACTAGACCA+ATAAGGCAGT\n"
            b"GGGGAAAATTTTGGGGCCCC\n"
            b"+\n"
            b"FFFFFFFFFFFFFFFFFFFF"
        )
    return fastq_path


@fixture(scope="module")
def fastq2():
    fastq_path = tempfile.NamedTemporaryFile().name
    with gzip.open(fastq_path, "wb") as f:
        f.write(
            b"@1 2:N:0:CACTAGACCA+ATAAGGCAGT\n"
            b"ATGCRYSWKMBDHVN\n"
            b"+\n"
            b"FFFFFFFFFFFFFFF\n"
            b"@2 2:N:0:CACTAGACCA+ATAAGGCAGT\n"
            b"GGGGCCCCAAAATTTTAAAA\n"
            b"+\n"
            b"FFFFFFFFFFFFFFFFFFFF"
        )
    return fastq_path