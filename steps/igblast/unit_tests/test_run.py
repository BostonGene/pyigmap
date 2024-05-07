import os
import subprocess
import tempfile
import gzip
import pathlib

from pytest import fixture
import pandas as pd

from logger import set_logger

logger = set_logger(name=__file__)

IGBLAST_STEP_DIR = pathlib.Path(__file__).parents[1]

@fixture(scope="module")
def fasta_bcr_with_all_alleles():
    fasta_file = tempfile.NamedTemporaryFile().name
    with gzip.open(fasta_file, "wt") as f:
        f.write(
            ">IGLV1-40*01\n"
            "CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACC\n"
            ">IGLV1-40*03\n"
            "GATTATTACTGCCAGTCCTATGACAGCAGCCTGAGTGGT"
        )
    return fasta_file


@fixture(scope="module")
def fasta_bcr_with_only_major_allele():
    fasta_file = tempfile.NamedTemporaryFile().name
    with gzip.open(fasta_file, "wt") as f:
        f.write(">IGLV1-40*01\n"
                "CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACC")
    return fasta_file


@fixture(scope="module")
def ref_with_all_alleles() -> str:
    return os.path.join(IGBLAST_STEP_DIR, "igblast.reference.all_alleles.tar.gz")


@fixture(scope="module")
def ref_with_major_allele() -> str:
    return os.path.join(IGBLAST_STEP_DIR, "igblast.reference.major_allele.tar.gz")


@fixture(scope="function")
def output_annotation_path() -> str:
    annotation_path = tempfile.NamedTemporaryFile().name
    with gzip.open(annotation_path, "wb") as f:
        f.write("".encode("utf-8"))
    return annotation_path


def run_command(command: list[str]):
    try:
        _ = subprocess.run(command, text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        raise Exception(f"Failed to run {command}")
    except Exception as e:
        raise Exception(f"Undefined error: {e}")


def docker_cmd(ref_with_major_allele, fasta_bcr_with_only_major_allele, output_annotation_path) -> list[str]:
    ref_basename = "ref.tar.gz"
    input_fasta_basename = "input.fasta.gz"
    output_annotation_basename = "annotation.tsv.gz"
    return [
        "docker", "run",
        "-v", f"{ref_with_major_allele}:/root/{ref_basename}",
        "-v", f"{fasta_bcr_with_only_major_allele}:/root/{input_fasta_basename}",
        "-v", f"{output_annotation_path}:/root/{output_annotation_basename}",
        "igblast",
        "--in-ref", f"/root/{ref_basename}",
        "--in-fasta", f"/root/{input_fasta_basename}",
        "--receptor", "BCR",
        "--organism", "human",
        "--out-annotation", f"/root/{output_annotation_basename}",
    ]


def read_annotation(annotation_path: str) -> pd.DataFrame:
    return pd.read_csv(annotation_path, sep="\t", compression="gzip")


def test_run_with_major_allele(ref_with_major_allele, fasta_bcr_with_only_major_allele, output_annotation_path):
    cmd = docker_cmd(ref_with_major_allele, fasta_bcr_with_only_major_allele, output_annotation_path)
    run_command(cmd)
    annotation = read_annotation(output_annotation_path)
    logger.info(f'Output annotation path: {output_annotation_path}')
    assert annotation["sequence_id"].to_list() == ["IGLV1-40*01"]
    assert annotation["v_call"].to_list() == ["IGLV1-40*01,IGLV1-50*01"]
    os.remove(output_annotation_path)


def test_run_with_all_alleles(ref_with_all_alleles, fasta_bcr_with_all_alleles, output_annotation_path):
    cmd = docker_cmd(ref_with_all_alleles, fasta_bcr_with_all_alleles, output_annotation_path)
    run_command(cmd)
    annotation = read_annotation(output_annotation_path)
    logger.info(f'Output annotation path: {output_annotation_path}')
    assert annotation["sequence_id"].to_list() == ["IGLV1-40*01", "IGLV1-40*03"]
    assert annotation["v_call"].to_list() == [
        "IGLV1-40*01,IGLV1-50*01",
        "IGLV1-40*01,IGLV1-40*02,IGLV1-40*03",
    ]
    os.remove(output_annotation_path)
