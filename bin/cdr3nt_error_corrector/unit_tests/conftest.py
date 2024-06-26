import json
import tempfile
import pandas as pd

from pytest import fixture


@fixture(scope='module')
def fastp_json() -> str:
    json_content = {
        "summary": {
            "fastp_version": "0.23.4",
            "before_filtering": {
                "total_reads": 20000  # sum of read1 and read2
            }
        }
    }
    fastp_json_path = tempfile.NamedTemporaryFile().name
    with open(fastp_json_path, 'w') as f:
        f.write(json.dumps(json_content))
    return fastp_json_path


@fixture(scope='module')
def calib_json() -> str:
    json_content = {
        "summary": {
            "before_filtering": {
                "total_reads": 10000  # sum of only read1 (or read2)
            }
        }
    }
    calib_json_path = tempfile.NamedTemporaryFile().name
    with open(calib_json_path, 'w') as f:
        f.write(json.dumps(json_content))
    return calib_json_path


@fixture(scope='module')
def empty_annotation() -> pd.DataFrame:
    columns = [
        'sequence', 'locus', 'stop_codon', 'vj_in_frame', 'v_frameshift', 'productive', 'v_call', 'j_call', 'junction',
        'junction_aa', 'v_support', 'j_support', 'v_sequence_start', 'v_sequence_end', 'j_sequence_start',
        'j_sequence_end', 'j_sequence_alignment_aa', 'pgen', 'duplicate_count'
    ]
    return pd.DataFrame(columns=columns)


@fixture(scope='module')
def empty_annotation_file(empty_annotation) -> str:
    annotation_path = tempfile.NamedTemporaryFile(suffix=".tsv").name
    empty_annotation.to_csv(annotation_path, index=False)
    return annotation_path