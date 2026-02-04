import logging
import os
import subprocess
import sys
import tempfile

LOGGER_FORMAT = '%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s'
logger = logging.getLogger()

OLGA_REPO_URL = 'https://github.com/statbiophys/OLGA.git'
OUTPUT_ARCHIVE_PATH = os.path.join('/tmp', 'olga-models.tar.gz')


def configure_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def run_and_check_with_message(
    cmd: list[str], fail_message: str, exit_on_error: bool = True, **subprocess_args
) -> subprocess.CompletedProcess[str] | None:
    logger.info(f'Running command: {" ".join(cmd)}')
    if 'stderr' not in subprocess_args:
        subprocess_args['stderr'] = subprocess.PIPE
    try:
        proc = subprocess.run(cmd, text=True, check=True, **subprocess_args)
        if proc.stderr:
            logger.warning(f'stderr: {proc.stderr.strip()}')
        return proc
    except subprocess.CalledProcessError as e:
        logger.critical(f'{fail_message} failed with code {e.returncode}')
        if e.stderr:
            logger.critical(f'stderr: {e.stderr.strip()}')
        if exit_on_error:
            sys.exit(1)
        return None


def clone_olga_repo(repo_dir: str) -> None:
    run_and_check_with_message(['git', 'clone', '--quiet', OLGA_REPO_URL, repo_dir], 'Cloning OLGA repository')


def archive_olga_models(models_dir: str, output_archive: str) -> None:
    run_and_check_with_message(
        ['tar', '-czf', output_archive, '-C', models_dir, '.'], 'Creating tar.gz archive with OLGA models'
    )
    logger.info(f'Archive with OLGA models here: {output_archive}')


def main():
    configure_logger()
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = os.path.join(tmpdir, 'OLGA')
        clone_olga_repo(repo_dir)

        models_dir = os.path.join(repo_dir, 'olga', 'default_models')
        archive_olga_models(models_dir, OUTPUT_ARCHIVE_PATH)

    logger.info('OLGA model archive creation completed successfully.')


if __name__ == '__main__':
    main()
