import gzip
import os
import subprocess
import tempfile
import pathlib

from pytest import fixture
import pandas as pd

from logger import set_logger

logger = set_logger(name=__file__)

STEP_DIR = pathlib.Path(__file__).parents[1]


@fixture(scope="module")
def annotation_tcr():
    annotation = tempfile.NamedTemporaryFile().name
    pd.DataFrame(
        data={'sequence': ['GGCGTATATGTGCGGGAGGGCAGTGGTATCAACGCAGAGTACGGCGTTACTGAGTTATGTGCTCTTGGGGGCGACGACCACGTTCCCATCTCCCTAACATGCCTGGTCTTACCTCCATTGGGCCTTGGCACGTGCTGCTCCCTCCTTGGGGACAGCTGTCTCTAGGGCTTCCCATGGCTGCAGCTCAGACCTGGCCTCATGGAGGGGTCTTATCTGAATGCCTGTGGTGGATGTCCTACCCACTCTCCCACCTTACACTCCTGCCCAGGGGCCTTAGCATCCCCTGGCACACCCCAGGCACCCTCCTGCCCGGG'],
              'locus': ['TRD'],
              'stop_codon': ['F'],
              'vj_in_frame': ['F'],
              'v_frameshift': ['F'],
              'productive': ['F'],
              'v_call': ['TRDV1*01'],
              'j_call': ['TRAJ1*01,TRAJ47*01,TRAJ47*02'],
              'junction': ['TGTGCTCTTGGGGGCGACGACCACGTTCCCATCTCCCTAACATGCCTGGTCTTACCTCCATT'],
              'junction_aa': ['CALGGDDHVPISLTCLVLPP'],
              'v_support': [1.573],
              'j_support': [19.65],
              'v_sequence_start': [55],
              'v_sequence_end': [70],
              'j_sequence_start': [109],
              'j_sequence_end': [116]
              }
    ).to_csv(annotation, sep='\t')
    return annotation


@fixture(scope="module")
def annotation_bcr():
    annotation = tempfile.NamedTemporaryFile().name
    pd.DataFrame(
        data={'sequence': ['TCCCATCATAAGGGGGAGGGCAGTGGTATCAACGCAGAGTACTCCCATTCTGAGCTATAAGGTCTTGGGGGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCAGTGGCCTCCAGTCTGACGATGAGGCTGATTATTATTGTGCGTCATGGGATGACAGCCTGAATGGCCGGCTGTTCGGCGGAGGGACCAAGTTGACCGTCCTGGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAAGGGACTTCTACCCGG'],
              'locus': ['IGL'],
              'stop_codon': ['F'],
              'vj_in_frame': ['T'],
              'v_frameshift': ['F'],
              'productive': ['T'],
              'v_call': ['IGLV1-36*01,IGLV1-44*01,IGLV1-44*02'],
              'j_call': ['IGLJ3*02'],
              'junction': ['TGTGCGTCATGGGATGACAGCCTGAATGGCCGGCTGTTC'],
              'junction_aa': ['CASWDDSLNGRLF'],
              'v_support': [2.915E-036],
              'j_support': [0.00000000001269],
              'v_sequence_start': [71],
              'v_sequence_end': [171],
              'j_sequence_start': [177],
              'j_sequence_end': [207]
              }
    ).to_csv(annotation, sep='\t')
    return annotation


@fixture(scope="module")
def olga_models() -> str:
    return os.path.join(STEP_DIR, "olga-models.tar.gz")


@fixture(scope="function")
def output_annotation_path() -> str:
    annotation_path = tempfile.NamedTemporaryFile().name
    with open(annotation_path, "w") as f:
        f.write("")
    return annotation_path


@fixture(scope="function")
def output_json_path() -> str:
    json_path = tempfile.NamedTemporaryFile().name
    with open(json_path, "w") as f:
        f.write("")
    return json_path


@fixture(scope="function")
def output_archive_path() -> str:
    archive_path = tempfile.NamedTemporaryFile().name
    with gzip.open(archive_path, 'wb') as f:
        f.write(''.encode('utf-8'))
    return archive_path


def run_command(command: list[str]):
    try:
        _ = subprocess.run(command, text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e.stderr)
        logger.critical(e.stdout)
        raise Exception(f"Failed to run {command}")
    except Exception as e:
        raise Exception(f"Undefined error: {e}")


def docker_cmd(olga_models, annotation_bcr, annotation_tcr, input_json, output_annotation_path,
               output_json_path, output_archive_path) -> list[str]:
    olga_models_basename = "olga-models.tar.gz"
    input_TCR_annotation_basename = "annotation.TCR.tsv"
    input_BCR_annotation_basename = "annotation.BCR.tsv"
    input_json_basename = "input.json"
    output_json_basename = "stat.json"
    output_annotation_basename = "corrected_annotation.tsv"
    output_archive_basename = "pyigmap.tar.gz"
    return [
        "docker", "run",
        "-v", f"{olga_models}:/root/{olga_models_basename}",
        "-v", f"{annotation_tcr}:/root/{input_TCR_annotation_basename}",
        "-v", f"{annotation_bcr}:/root/{input_BCR_annotation_basename}",
        "-v", f"{input_json}:/root/{input_json_basename}",
        "-v", f"{output_annotation_path}:/root/{output_annotation_basename}",
        "-v", f"{output_json_path}:/root/{output_json_basename}",
        "-v", f"{output_archive_path}:/root/{output_archive_basename}",
        "cdr3nt-error-corrector-tool",
        "--in-tcr-annotation", f"/root/{input_TCR_annotation_basename}",
        "--in-bcr-annotation", f"/root/{input_BCR_annotation_basename}",
        "--pgen-threshold", str(0),
        "--only-productive",
        "--only-functional",
        "--remove-chimeras",
        "--clonotype-collapse-factor", str(0.05),
        "--olga-models", f"/root/{olga_models_basename}",
        "--out-corrected-annotation", f"/root/{output_annotation_basename}",
        "--in-json", f"/root/{input_json_basename}",
        "--out-json", f"/root/{output_json_basename}",
        "--out-archive", f"/root/{os.path.basename(output_archive_path)}"]


def read_annotation(annotation_path: str) -> pd.DataFrame:
    return pd.read_csv(annotation_path, sep="\t")


def test_run_with_calib(olga_models, annotation_bcr, annotation_tcr, calib_json, output_annotation_path,
               output_json_path, output_archive_path):
    cmd = docker_cmd(olga_models, annotation_bcr, annotation_tcr, calib_json, output_annotation_path,
               output_json_path, output_archive_path)
    run_command(cmd)
    annotation = read_annotation(output_annotation_path)
    logger.info(f'Output annotation path: {annotation}')
    assert not len(annotation)


