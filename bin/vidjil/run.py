import argparse
import glob
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from multiprocessing import Pool
from typing import Dict, List, Optional

LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger()

CPU_COUNT = os.cpu_count()
TEMPDIR = os.environ.get('TEMPDIR', tempfile.gettempdir())
REF_DIR = os.path.join(TEMPDIR, 'vdj_ref')
SEQKIT_TEMPDIR = os.path.join(TEMPDIR, 'seqkit')
VIDJIL_TEMPDIR = os.path.join(TEMPDIR, 'vidjil')
VIDJIL_REF_GERMLINE_PATH = os.path.join(REF_DIR, 'homo-sapiens.g')
VIDJIL_OUT_FASTA_PATTERN = '*detected.vdj.fa'
VIDJIL_OUT_JSON_PATTERN = '*.vidjil'
SEQKIT_OUT_EXTENSION = '.fastq.gz'
SEQKIT_OUT_FQ_PATTERN = f'*{SEQKIT_OUT_EXTENSION}'


def configure_logger(fmt: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fastq', help='Input FASTQ')
    parser.add_argument('--vdj-ref', help='V(D)J reference', required=True)
    parser.add_argument(
        '--e-value',
        type=float,
        help='E-value threshold to relax V–J window filtering and capture more potential recombination signals'
    )
    parser.add_argument('--debug', help='Enables saving logs', action="store_true")
    parser.add_argument('--out-json', help='JSON output summary with QC metrics')
    parser.add_argument('--out-fasta', help='Output FASTA with detected V(D)J segments', required=True)

    return parser.parse_args()


def print_error_message(msg: Optional[str]) -> None:
    logger.warning(f"stderr: {msg or '<EMPTY>'}")


def run_and_check_with_message(
    cmd: List[str],
    fail_message: str,
    exit_on_error: bool = True,
    return_proc: bool = False,
    **subprocess_args
) -> Optional[subprocess.CompletedProcess[str]]:
    logger.info(f"Running command {' '.join(cmd)}")
    if 'stderr' not in subprocess_args:
        subprocess_args['stderr'] = subprocess.PIPE
    try:
        proc = subprocess.run(cmd, text=True, check=True, **subprocess_args)
        print_error_message(proc.stderr)
        if return_proc:
            return proc
    except subprocess.CalledProcessError as e:
        logger.critical(f"{fail_message} failed with code {e.returncode}")
        print_error_message(e.stderr)
        if exit_on_error:
            logger.critical(f"{exit_on_error=}, now exiting.")
            sys.exit(1)


@dataclass
class VidjilReads:
    clones: Dict[str, int] = field(default_factory=dict)
    germline: Dict[str, int] = field(default_factory=dict)
    segmented: int = 0
    total: int = 0

    def add(self, other: dict):
        for section in ["clones", "germline"]:
            for k, v in other.get(section, {}).items():
                self_section = getattr(self, section)
                self_section[k] = self_section.get(k, 0) + (v[0] if isinstance(v, list) else v)
        for k in ["segmented", "total"]:
            val = other.get(k, [0])
            setattr(self, k, getattr(self, k) + (val[0] if isinstance(val, list) else val))

    def to_dict(self):
        return {
            "reads": {
                "clones": self.clones,
                "germline": self.germline,
                "detected": self.segmented,
                "total": self.total
            },
            "vdj_detected": self.segmented > 0
        }


def merge_and_compress_fasta(fasta_files: List[str], output_fasta: str) -> None:
    non_empty_files = [f for f in fasta_files if os.path.getsize(f) > 0]
    if not non_empty_files:
        logger.warning("All FASTA files are empty — skipping output file creation")
        return

    logger.info(f"Merging {len(fasta_files)} FASTA files into compressed output: {output_fasta}")

    cmd = ['pigz', '-p', str(CPU_COUNT), '-c']

    with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=open(output_fasta, 'wb')) as proc:
        for fasta in fasta_files:
            logger.debug(f"Appending {fasta}")
            with open(fasta, 'rb') as rfd:
                shutil.copyfileobj(rfd, proc.stdin)
        proc.stdin.close()
        proc.wait()

    if proc.returncode != 0:
        logger.critical(f"pigz compression failed with return code {proc.returncode}")
        sys.exit(1)


def merge_vidjil_json(files: list[str]) -> dict:
    merged = VidjilReads()
    logger.info(f"Merging {len(files)} Vidjil JSON files")
    for path in files:
        logger.debug(f"Reading JSON file {path}")
        with open(path, "rt") as f:
            data = json.load(f)
            merged.add(data.get("reads", {}))
    logger.debug("Merging completed successfully")

    return merged.to_dict()


def save_results(fasta_files: list[str], json_files: list[str], out_fasta: str, out_json: str, debug: bool) -> None:
    logger.info(f"Saving results to {out_fasta} and {out_json}")

    result_dict = merge_vidjil_json(json_files)

    if debug:
        logger.info(f"Writing merged JSON results to {out_json}")
        with open(out_json, 'w') as f:
            json.dump(result_dict, f, indent=2)

    if result_dict["vdj_detected"]:
        merge_and_compress_fasta(fasta_files, out_fasta)
    else:
        logger.warning("Detected cdr3 reads == 0, skipping out_fasta creation.")


def unpack_reference(archive: str) -> None:
    logger.info(f"Unpacking reference archive {archive}")
    os.makedirs(REF_DIR, exist_ok=True)
    cmd = ['tar', '-xvf', archive, '-C', REF_DIR]
    run_and_check_with_message(cmd, 'Unpack reference')


def split_fastq(input_fastq: str) -> List[str]:
    logger.info(f"Splitting FASTQ file {input_fastq}")
    os.makedirs(SEQKIT_TEMPDIR, exist_ok=True)
    run_and_check_with_message([
        'seqkit', 'split2',
        '--threads', str(CPU_COUNT),
        '--by-part', str(CPU_COUNT),
        '-O', SEQKIT_TEMPDIR,
        input_fastq
    ], 'seqkit')
    chunks = sorted(glob.glob(os.path.join(SEQKIT_TEMPDIR, SEQKIT_OUT_FQ_PATTERN)))
    logger.info(f"Created {len(chunks)} FASTQ chunks")
    return chunks


def process_fq_chunk(fastq_chunk_path: str, e_value: float):
    logger.info(f"Processing chunk {fastq_chunk_path}")

    cmd = [
        'vidjil-algo',
        '--germline', VIDJIL_REF_GERMLINE_PATH,
        '--out-detected',
        '--dir', VIDJIL_TEMPDIR,
        '--clean-memory',
        '--base', os.path.basename(fastq_chunk_path),
        '-c', 'clones',
        '--no-airr',
        '--verbose',
        '-e', str(e_value),
        fastq_chunk_path
    ]

    run_and_check_with_message(cmd, 'vidjil')

    logger.info(f"Chunk {fastq_chunk_path} processed.")


def main() -> None:
    configure_logger()
    args = parse_args()
    logger.info(f"Started with args: {vars(args)}")

    unpack_reference(args.vdj_ref)
    fq_chunks = split_fastq(args.in_fastq)

    logger.info(f"FASTQ chunks to process: {fq_chunks}")

    with Pool(CPU_COUNT) as pool:
        pool.starmap(process_fq_chunk, [(chunk, args.e_value) for chunk in fq_chunks])

    out_fasta_files = sorted(glob.glob(os.path.join(VIDJIL_TEMPDIR, VIDJIL_OUT_FASTA_PATTERN)))
    logger.info(f"FASTA outputs collected: {out_fasta_files}")

    out_json_files = sorted(glob.glob(os.path.join(VIDJIL_TEMPDIR, VIDJIL_OUT_JSON_PATTERN)))
    logger.info(f"JSON outputs collected: {out_json_files}")

    save_results(out_fasta_files, out_json_files, args.out_fasta, args.out_json, args.debug)
    logger.info("Run is completed successfully.")


if __name__ == '__main__':
    main()
