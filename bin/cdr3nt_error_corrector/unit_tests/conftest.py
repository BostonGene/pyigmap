import json
import tempfile
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