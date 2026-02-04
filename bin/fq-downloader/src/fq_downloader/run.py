import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from subprocess import CalledProcessError, CompletedProcess
from typing import Any

from ffq.ffq import ffq_experiment, ffq_run

from fq_downloader.metadata import Experiment, FileSource, Run

LOGGER_FORMAT = '%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s'
logger = logging.getLogger()

CPU_COUNT = os.cpu_count()
TEMPDIR = os.environ.get('TEMPDIR', tempfile.gettempdir())
CHUNK_SIZE = 25 * 1024 * 1024  # 25MB
FFQ_RETRIES = 10
RETRY_DELAY_SECONDS = 3

RUN_PREFIXES = ('SRR', 'ERR', 'DRR')
EXP_PREFIXES = ('SRX', 'ERX', 'DRX')


def configure_logger(logger_format: str = LOGGER_FORMAT) -> logging.Logger:
    """Configure root logger once."""
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--identifier', required=True, help='SRX*/ERX*/DRX* (experiment) or SRR*/ERR*/DRR* (run)')
    parser.add_argument('--out-fq1', required=True, help='Output path for merged/compressed R1 (or single-end)')
    parser.add_argument('--out-fq2', help='Output path for merged/compressed R2 (paired-end only)')
    parser.add_argument('--download-retries', type=int, default=100, help='Number of retries for downloads')
    return parser.parse_args()


def exit_with_error(message: str | None) -> None:
    """Log critical error and exit."""
    if not message:
        message = 'Empty error message'
    logger.critical(message)
    sys.exit(1)


def print_error_message(error_message: str | None) -> None:
    """Log captured stderr of a subprocess as warning (if any)."""
    if not error_message:
        error_message = '<EMPTY>'
    logger.warning(f'stderr: {error_message}')


def run_and_check_with_message(
    cmd: list[str], fail_message: str, exit_on_error: bool = True, return_proc: bool = False, **subprocess_args
) -> CompletedProcess[str] | CompletedProcess[Any] | CalledProcessError | None:
    """
    Run a command, capture stderr, and optionally return the CompletedProcess.
    If the command fails and exit_on_error is True, exit the program.
    """
    logger.info(f'Running command {" ".join(cmd)}')
    if 'stderr' not in subprocess_args:
        subprocess_args['stderr'] = subprocess.PIPE
    try:
        proc = subprocess.run(cmd, text=True, check=True, **subprocess_args)
        print_error_message(proc.stderr)
        if return_proc:
            return proc
    except subprocess.CalledProcessError as e:
        logger.critical(f'{fail_message} failed with code {e.returncode}')
        print_error_message(e.stderr)
        if return_proc:
            return e
        if exit_on_error:
            logger.critical(f'{exit_on_error=}, exiting.')
            sys.exit(1)
    return None


def _all_gz(files: list[str]) -> bool:
    """Return True if all paths end with .gz."""
    return all(f.endswith('.gz') for f in files)


def _any_gz(files: list[str]) -> bool:
    """Return True if any path ends with .gz."""
    return any(f.endswith('.gz') for f in files)


def merge_and_compress_fastqs(fastq_files: list[str], output_fastq_gz: str) -> None:
    """
    Robust merging/compression strategy:
      - If all inputs are .fastq.gz: concatenate gzip members verbatim (valid multi-member gzip).
      - If all inputs are plain .fastq: stream through pigz to produce a single gzip.
      - If mixed: decompress gz inputs on the fly, then compress the unified stream via pigz.
    """
    if not fastq_files:
        exit_with_error('merge_and_compress_fastqs: empty input list')

    logger.info(f'Merging mixed FASTQ set and compressing to {output_fastq_gz} with pigz')
    cmd = ['pigz', '-p', str(CPU_COUNT), '-c']
    with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=open(output_fastq_gz, 'wb')) as proc:
        assert proc.stdin is not None
        for fastq in fastq_files:
            with open(fastq, 'rb') as rfd:
                shutil.copyfileobj(rfd, proc.stdin, length=CHUNK_SIZE)
        proc.stdin.close()
        proc.wait()

    if proc.returncode != 0:
        exit_with_error(f'pigz compression failed with return code {proc.returncode}')


def save_outputs(fq1_runs: list[str], fq2_runs: list[str | None], out_fq1: str, out_fq2: str | None) -> None:
    """Merge and write final outputs. R2 is optional (single-end case)."""
    merge_and_compress_fastqs(fq1_runs, out_fq1)
    if fq2_runs and out_fq2:
        fq2_clean = [p for p in fq2_runs if p]
        if fq2_clean:
            merge_and_compress_fastqs(fq2_clean, out_fq2)


def convert_sra_to_fastq(sra_file: str, accession: str, ngc: str | None = None) -> tuple[str, str | None]:
    """
    Convert a downloaded SRA file to FASTQ using fasterq-dump.
    Returns (r1, r2_or_None). For single-end, r2 is None.
    """
    logger.info(f"Converting '{sra_file}' to FASTQ via fasterq-dump")
    cmd = [
        'fasterq-dump',
        sra_file,
        '--threads',
        str(CPU_COUNT),
        '-O',
        TEMPDIR,
        '--temp',
        TEMPDIR,
        '-x',
        '--format',
        'fastq',
        '--size-check',
        'off',
        '--seq-defline',
        f'@{accession}.$si $sn $ri',
        '--qual-defline',
        '+',
    ]
    if ngc:
        cmd += ['--ngc', ngc]
    run_and_check_with_message(cmd, 'fasterq-dump')

    # fasterq-dump outputs either "<prefix>.fastq" (single) or "<prefix>_1.fastq" + "<prefix>_2.fastq"
    prefix = os.path.basename(sra_file)
    fq_single = os.path.join(TEMPDIR, f'{prefix}.fastq')
    if os.path.exists(fq_single):
        return fq_single, None
    return (os.path.join(TEMPDIR, f'{prefix}_1.fastq'), os.path.join(TEMPDIR, f'{prefix}_2.fastq'))


def download_file(url: str, download_retries: int) -> str:
    with tempfile.NamedTemporaryFile(dir=TEMPDIR, delete=False) as f:
        out = f.name
    for retry in range(1, download_retries + 1):
        cmd = run_and_check_with_message(
            ['axel', '-a', '--num-connections', str(CPU_COUNT), url, f'--output={out}', '--insecure', '--quiet'],
            'axel',
            exit_on_error=False,
            return_proc=True,
        )
        if cmd is None or cmd.returncode != 0:
            if cmd is not None:
                logger.info(cmd.stderr)
            logger.info(f'Attempt {retry}/{download_retries} failed.')
            if retry < download_retries:
                time.sleep(RETRY_DELAY_SECONDS)
                continue
            exit_with_error(f'Failed to download after {download_retries} retries.')
        return out


def handle_run_files(run: Run, download_retries: int) -> tuple[list[str], list[str | None]]:
    """
    Resolve input files for a single run strictly from NCBI source:
      - Find SRA entries in run.files['ncbi'].
      - Download the first available SRA.
      - Convert it to FASTQ via fasterq-dump.
    Returns lists for R1 and (optional) R2 paths.
    """
    fq1_paths: list[str] = []
    fq2_paths: list[str | None] = []

    ncbi_files = run.files.get(FileSource.NCBI, [])
    for file in ncbi_files:
        if file.filetype == 'sra':
            sra_file = download_file(file.url, download_retries)
            fq1, fq2 = convert_sra_to_fastq(sra_file, file.accession)
            fq1_paths.append(fq1)
            if fq2:
                fq2_paths.append(fq2)
            # One SRA per run is enough
            break

    return fq1_paths, fq2_paths


def fetch_metadata(identifier: str) -> dict:
    """Fetch experiment-level metadata with retries (SRX/ERX/DRX)."""
    for retry in range(1, FFQ_RETRIES + 1):
        try:
            return ffq_experiment(identifier)
        except Exception as e:
            logger.warning(f'Retry {retry}/{FFQ_RETRIES} for {identifier} failed with: {e}')
            code = getattr(getattr(e, 'response', None), 'status_code', None)
            if code:
                logger.warning(f'HTTP {code} for {identifier}, retrying...')
            if retry < FFQ_RETRIES:
                time.sleep(RETRY_DELAY_SECONDS)
            else:
                exit_with_error(f'Exceeded retry attempts for {identifier}: {e}')
    return {}


def fetch_metadata_run(identifier: str) -> dict:
    """Fetch run-level metadata with retries (SRR/ERR/DRR)."""
    for retry in range(1, FFQ_RETRIES + 1):
        try:
            # ffq_run returns a dict keyed by the provided id
            return ffq_run(identifier)
        except Exception as e:
            logger.warning(f'Retry {retry}/{FFQ_RETRIES} for run {identifier} failed with: {e}')
            code = getattr(getattr(e, 'response', None), 'status_code', None)
            if code:
                logger.warning(f'HTTP {code} for {identifier}, retrying...')
            if retry < FFQ_RETRIES:
                time.sleep(RETRY_DELAY_SECONDS)
            else:
                exit_with_error(f'Exceeded retry attempts for run {identifier}: {e}')
    return {}


def process_experiment(identifier: str, download_retries: int) -> tuple[list[str], list[str | None]]:
    """
    Process an experiment id:
      - Parse metadata to Experiment object.
      - For each run, resolve NCBI SRA → fasterq-dump.
    """
    metadata = fetch_metadata(identifier)
    experiment = Experiment.from_dict(metadata)

    fq1_paths: list[str] = []
    fq2_paths: list[str | None] = []
    for run in experiment.runs.values():
        r1_list, r2_list = handle_run_files(run, download_retries)
        fq1_paths += r1_list
        fq2_paths += r2_list

    if not fq1_paths and not fq2_paths:
        exit_with_error(f'No FASTQ could be produced for {identifier} from NCBI SRA.')

    return fq1_paths, fq2_paths


def process_run(identifier: str, download_retries: int) -> tuple[list[str], list[str | None]]:
    """
    Process a single run id:
      - Fetch run metadata.
      - Resolve NCBI SRA → fasterq-dump.
    """
    metadata = fetch_metadata_run(identifier)
    if not metadata:
        exit_with_error(f'Empty metadata for run {identifier}')

    run = Run.from_dict(metadata)

    fq1_paths, fq2_paths = handle_run_files(run, download_retries)

    if not fq1_paths and not fq2_paths:
        exit_with_error(f'No FASTQ could be produced for run {identifier} from NCBI SRA.')

    return fq1_paths, fq2_paths


def is_run_id(identifier: str) -> bool:
    """Heuristic: SRR/ERR/DRR → run."""
    up = identifier.upper()
    return up.startswith(RUN_PREFIXES)


def is_experiment_id(identifier: str) -> bool:
    """Heuristic: SRX/ERX/DRX → experiment."""
    up = identifier.upper()
    return up.startswith(EXP_PREFIXES)


def main():
    configure_logger()
    args = parse_args()
    logger.info(f'Starting with args: {vars(args)}')

    fq1_paths: list[str]
    fq2_paths: list[str | None]

    if is_run_id(args.identifier):
        fq1_paths, fq2_paths = process_run(args.identifier, args.download_retries)
    elif is_experiment_id(args.identifier):
        fq1_paths, fq2_paths = process_experiment(args.identifier, args.download_retries)
    else:
        exit_with_error('Unknown id prefix.')
        return  # Unreachable, but helps type checker

    save_outputs(fq1_paths, fq2_paths, args.out_fq1, args.out_fq2)
    logger.info('Run completed successfully.')


if __name__ == '__main__':
    main()
