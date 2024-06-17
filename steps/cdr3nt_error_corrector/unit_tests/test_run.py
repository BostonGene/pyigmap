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

container_engine = "podman" if os.environ.get("USE_PODMAN", False) else "docker"


@fixture(scope="module")
def annotation_tcr():
    annotation = tempfile.NamedTemporaryFile().name
    pd.DataFrame(
        data={'sequence': ['TAGAAGTCTTTTTTATGAGACGGTGACCGTGGAAACGGGAGTTACACAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCTTTGAAATGTGAACAACATCTGGGGCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAGAAGCCACTGGAGCTCATGTTTGTCTACAACTTTAAAGAACAGACTGAAAACAACAGTGTGCCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATTCCGTCACCTACACACCCTGCAGCCAGAAGACTCGGCCCTGTATCTCTGTGCCAGCAGCCAAGTCACATCGGGTAGCGGGGAGCAGTTCTTCGGGCCAGGGACACGGCTCACCGTGCTAGAGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCCTTCCCGACTCCATCACCCCCAAGAAAAAAGAACTCAGCATTCTAGTACTCTGCGTTGATACCACTGCCCTCGGAAAAAGACTTCTA'],
              'locus': ['TRB'],
              'stop_codon': ['F'],
              'vj_in_frame': ['T'],
              'v_frameshift': ['F'],
              'productive': ['T'],
              'v_call': ['TRBV4-2*01'],
              'j_call': ['TRBJ2-1*01'],
              'junction': ['TGTGCCAGCAGCCAAGTCACATCGGGTAGCGGGGAGCAGTTCTTC'],
              'junction_aa': ['CASSQVTSGSGEQFF'],
              'v_support': [5.841E-126],
              'j_support': [1.09E-017],
              'v_sequence_start': [32],
              'v_sequence_end': [317],
              'j_sequence_start': [335],
              'j_sequence_end': [374]
              }
    ).to_csv(annotation, sep='\t')
    return annotation


@fixture(scope="module")
def annotation_bcr():
    annotation = tempfile.NamedTemporaryFile().name
    pd.DataFrame(
        data={'sequence': ['GTGCCAGACAGATTAGAGGGCAGTGGTATCAACGCAGAGTACGTGCCTACTGAGGTACAGATTCTTGGGGGATGCTTTCTGAGAGTCATGGATCTCATGTGCAAGAAAATGAAGCACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTCCCAACTACAGTTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCAGTGTCTCTGTTGGCTTCATCGACATTGAAGGTTATCACTGGGGCTGGATCCGCCAGTCCCCAGGGGCGGCCCTGGAGGGGCTTGGGAGCATCGATTATCGTGACACTTCCTGGCACAACCCGTCCCTCGGGAGGCGAGTCGCCCTGTCCATGGACACGCCCAAGAACAACTTCTCTCTGCAGTTGACCTCCGTGACCGCCGCAGACACGGCTGTGTATTTCTGTGTGAGACATAAACCTATGGTCCAGGGCGGCGTCGACGTCTGGGGCCAAGGAACCATGGTCACCGTCTCTTATCTGTCTGGCAC',
                           'ATTTTACACTGAAAGTCAGCCGAGGGGAGGCTGAGGATGTTGGACTTTATTACTGCGCACAAGATGCACAAGATCGTCCGCTCACTGTTGGCGGAGGGACCAAGGTGGAGATCAGACGTGAGTGCACTTTCCTAATGCTTTTCTTATACAG',
                           'GAGGATGTTGGACTTTATTACTGCGCATAAGATGCACAAGATCGTCCGCTCACTGTTGGCGGAGGGACCAAGGTGGAGATCAGACGATTTTCTCTGCATCGGTCAGGTTAGTGATATTAACAGCGAAAAGAGACTTTTGTTAAGGACTC',
                           'CAGGATGTTGGACTTTATTGCTGCGCATAAGATGCACAAGATCGTCCGCTCACTGTTGGCGGAGGGACCAAGGTGGAGATCAGACGATTTTCTCTGCATCGGTCAGGTTAGTGATATTAACAGCGAAAAGAGACTTTTGTTAAGGACTCAG'],
              'locus': ['IGH', 'IGK', 'IGK', 'IGK'],
              'stop_codon': ['F', 'F', 'F', 'F'],
              'vj_in_frame': ['T', 'T', 'T', 'T'],
              'v_frameshift': ['F', 'F', 'F', 'F'],
              'productive': ['T', 'T', 'T', 'T'],
              'v_call': ['IGHV4-39*01', 'IGKV2D-26*01', 'IGKV2D-26*01', 'IGKV2D-26*01'],
              'j_call': ['IGHJ6*02', 'IGKJ4*01', 'IGKJ4*01', 'IGKJ4*01'],
              'c_call': [None, 'IGHM', 'IGHM', 'IGHD'],
              'd_call': [None, None, None, None],
              'junction': ['TGTGTGAGACATAAACCTATGGTCCAGGGCGGCGTCGACGTCTGG',
                           'TGCGCACAAGATGCACAAGATCGTCCGCTCACTGTT',
                           'TGCGCACAAGATGCACAAGATCGTCCGCTCACTGTT',
                           'TGCGCACAAGATGCACAAGATCGTCCGCTCACTGTT'],
              'junction_aa': ['CVRHKPMVQGGVDVW', 'CAQDAQDRPLTV', 'CAQDAQDRPLTV', 'CAQDAQDRPLTV'],
              'v_support': [9.911E-086, 0, 0, 0],
              'j_support': [0.000000000005645, 0, 0, 0],
              'v_sequence_start': [166, 1, 1, 1],
              'v_sequence_end': [464, 79, 47, 47],
              'j_sequence_start': [490, 80, 48, 48],
              'j_sequence_end': [524, 114, 82, 82],
              'd_sequence_start': [None, None, None, None],
              'd_sequence_end': [None, None, None, None],
              'c_sequence_start': [None, 153, 102, 100],
              'c_sequence_end': [None, 298, 147, 145]
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
    tcr_annotation_basename = "annotation.TCR.tsv"
    bcr_annotation_basename = "annotation.BCR.tsv"
    input_json_basename = "input.json"
    output_json_basename = "stat.json"
    output_annotation_basename = "corrected_annotation.tsv"
    output_archive_basename = "pyigmap.tar.gz"
    return [
        container_engine, "run",
        "-v", f"{olga_models}:/root/{olga_models_basename}",
        "-v", f"{annotation_tcr}:/root/{tcr_annotation_basename}",
        "-v", f"{annotation_bcr}:/root/{bcr_annotation_basename}",
        "-v", f"{input_json}:/root/{input_json_basename}",
        "-v", f"{output_annotation_path}:/root/{output_annotation_basename}",
        "-v", f"{output_json_path}:/root/{output_json_basename}",
        "-v", f"{output_archive_path}:/root/{output_archive_basename}",
        "cdr3nt_error_corrector-tool",
        "--in-annotation", f"/root/{tcr_annotation_basename}", f"/root/{bcr_annotation_basename}",
        "--remove-chimeras",
        "--top-c-call",
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

    logger.info(f"{annotation['locus']}")
    assert set(annotation['locus']) == {'IGH', 'IGK', 'TRB'}

    annotation['c_call'] = annotation['c_call'].fillna('')
    logger.info(f"{annotation['c_call']}")
    assert set(annotation['c_call']) == {'', 'IGHM', 'IGHD'}
