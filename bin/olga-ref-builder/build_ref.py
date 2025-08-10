import argparse
import os
import shutil
import subprocess
import tarfile
import tempfile
import logging

OLGA_REPO_URL = "https://github.com/statbiophys/OLGA.git"

def configure_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download OLGA model files and create a tar.gz archive."
    )
    parser.add_argument(
        "--out-archive",
        help="Full path to the output archive (e.g., /tmp/olga-models.tar.gz)"
    )
    return parser.parse_args()


def exit_with_error(message: Optional[str]) -> None:
    if message is None:
        message = "Empty message for error!"
    logger.critical(message)
    sys.exit(1)


def print_error_message(error_message: Optional[str]) -> None:
    if not error_message:
        error_message = '< EMPTY >'
    logger.warning(f'stderr: {error_message}')


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


def main():
    args = parse_args()

    repository_dir = tempfile.mkdtemp(prefix="OLGA_")

    cmd = ["git", "clone", "--quiet", OLGA_REPO_URL, repository_dir]
    run_and_check_with_message(cmd, 'git')

    models_dir = os.path.join(repository_dir, "olga", "default_models")

    with tarfile.open(args.out_archive, "w:gz") as tar:
        tar.add(models_dir, arcname=".")

    logger.info(f"Archive with OLGA models here: {args.output_archive}")


if __name__ == "__main__":
    main()