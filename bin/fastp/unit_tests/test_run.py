import gzip
import os
import subprocess
import tempfile
import pathlib

from logger import set_logger

logger = set_logger(name=__file__)

STEP_DIR = pathlib.Path(__file__).parents[1]

container_engine = "podman" if os.environ.get("USE_PODMAN", False) else "docker"


def output_file_path(make_gzip=False) -> str:
    suffix = ".gz" if make_gzip else ""
    file_path = tempfile.NamedTemporaryFile(suffix=suffix).name
    with gzip.open(file_path, "wb") if make_gzip else open(file_path, "w") as f:
        f.write(b"" if make_gzip else "")
    return file_path


def run_command(command: list[str]):
    try:
        _ = subprocess.run(command, text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        raise Exception(f"Failed to run {command}")
    except Exception as e:
        raise Exception(f"Undefined error: {e}")


def docker_cmd(in_fq1, in_fq2, out_fq12_path, output_html_path, output_json_path) -> list[str]:
    out_fq12_basename = "R12.fastq.gz"
    out_html_basename = "fastp.html"
    out_json_basename = "fastp.json"
    return [
        container_engine, "run", "--rm",
        "-v", f"{in_fq1}:{in_fq1}",
        "-v", f"{in_fq2}:{in_fq2}",
        "-v", f"{out_fq12_path}:/root/{out_fq12_basename}",
        "-v", f"{output_html_path}:/root/{out_html_basename}",
        "-v", f"{output_json_path}:/root/{out_json_basename}",
        "fastp-tool",
        "--in-fq1", in_fq1,
        "--in-fq2", in_fq2,
        "--out-fq12", f"/root/{out_fq12_basename}",
        "--html", f"/root/{out_html_basename}",
        "--json", f"/root/{out_json_basename}",
        "--disable-filters", "length_filtering", "quality_filtering", "adapter_trimming",
        "--mock-merge-reads", "--inner-distance-size", str(1)
    ]


def test_run(fastq1, fastq2):
    out_fq12_path = output_file_path(make_gzip=True)
    output_html_path = output_file_path(make_gzip=False)
    output_json_path = output_file_path(make_gzip=False)
    cmd = docker_cmd(fastq1, fastq2, out_fq12_path, output_html_path, output_json_path)
    run_command(cmd)
    with gzip.open(out_fq12_path, "rb") as f:
        reads = [line.strip() for line in f.readlines()]
    assert reads == [
        b'@3/1 merged_50_18',
        b'CCCGAACTCTGCCAGTCTGGAGCCCGGAGCTGAAGTAGGATTAGCCTCCAATGTGACTTCCAAGTCTG',
        b'+',
        b'@C@DFFFFHGHHFJJGBGHIEHEIIEIGE<;FHAHFGGECGGIGIIJGJIHGHHFDFHHGFFEFD@@B',
        b'@4/1 merged_35_0',
        b'TTGGTGTATATGTTGTAATTGAGATTGCTCGGGGG',
        b'+',
        b'+==A+?<A<=CB@AAAB4>ABBBAABABB4A<:A6',
        b'@5/1 merged_35_0',
        b'GTTGGGGCCCTGGCCTTTTCAGCTGCGGATCAGGG',
        b'+',
        b'C@CFFFFFHHHHHJJJJJJJJJJIJJII>@DGIIJ',
        b'@1/1 mock_merged_15_15',
        b'ATGCRYSWKMBDHVNNAAAAAAAAAAAGCAT',
        b'+',
        b'FFFFFFFFFFFFFFF#FFFFFFFFFFFFFFF',
        b'@2/1 mock_merged_20_20',
        b'GGGGAAAATTTTGGGGCCCCNTTTTAAAATTTTGGGGCCCC',
        b'+',
        b'FFFFFFFFFFFFFFFFFFFF#FFFFFFFFFFFFFFFFFFFF']
