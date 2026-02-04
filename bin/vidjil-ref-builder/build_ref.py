import argparse
import os
import subprocess
import sys
import tempfile
import shutil
import logging
from typing import List, Optional

LOGGER_FORMAT = "%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s"
logger = logging.getLogger()

VIDJIL_REPO_URL = "https://gitlab.inria.fr/vidjil/vidjil.git"
VIDJIL_VERSION = "release-2024.02"
REFERENCE_ARCHIVE_PATH = os.path.join("/tmp", "vidjil.germline.tar.gz")


def configure_logger(logger_format: str = LOGGER_FORMAT) -> None:
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out-archive",
        help="Full path to the output archive (e.g., /tmp/vidjil.germline.tar.gz)",
        required=True
    )
    return parser.parse_args()


def run_and_check_with_message(
    cmd: List[str],
    fail_message: str,
    exit_on_error: bool = True,
    **subprocess_args
) -> Optional[subprocess.CompletedProcess[str]]:
    logger.info(f"Running command: {' '.join(cmd)}")
    if 'stderr' not in subprocess_args:
        subprocess_args['stderr'] = subprocess.PIPE
    try:
        proc = subprocess.run(cmd, text=True, check=True, **subprocess_args)
        if proc.stderr:
            logger.warning(f"stderr: {proc.stderr.strip()}")
        return proc
    except subprocess.CalledProcessError as e:
        logger.critical(f"{fail_message} failed with code {e.returncode}")
        if e.stderr:
            logger.critical(f"stderr: {e.stderr.strip()}")
        if exit_on_error:
            sys.exit(1)
        return None


def clone_vidjil_repo(repo_dir: str) -> None:
    run_and_check_with_message(
        ["git", "clone", "-b", VIDJIL_VERSION, VIDJIL_REPO_URL, repo_dir],
        "Cloning vidjil repository"
    )


def make_vidjil_germline(repo_dir: str) -> None:
    germline_dir = os.path.join(repo_dir, "germline")
    run_and_check_with_message(
        ["make", "germline"],
        "Building vidjil germline",
        cwd=germline_dir
    )


def merge_presets(germline_dir: str) -> None:
    main = os.path.join(germline_dir, "homo-sapiens.g")
    iso = os.path.join(germline_dir, "homo-sapiens-isotypes.g")
    merged = os.path.join(germline_dir, "homo-sapiens.full.g")
    run_and_check_with_message(
        ["jq", "-s", ".[0] * .[1]", main, iso],
        "Merging homo-sapiens presets",
        stdout=open(merged, "w")
    )
    os.replace(merged, main)


def create_vidjil_archive(germline_dir: str, archive_path: str) -> None:
    files = [
        "gallus-gallus.g", "homo-sapiens-cd.g", "homo-sapiens.g",
        "homo-sapiens-isotypes.g", "homo-sapiens-isoforms.g", "mus-musculus.g",
        "rattus-norvegicus.g", "sus-scrofa.g",
        "homo-sapiens", "gallus-gallus", "mus-musculus", "rattus-norvegicus", "sus-scrofa"
    ]
    cmd = ["tar", "-czf", archive_path, "-C", germline_dir] + files
    run_and_check_with_message(cmd, "Creating tar.gz archive")
    logger.info(f"Archive with vidjil reference here: {archive_path}")


def main():
    configure_logger()
    args = parse_args()
    with tempfile.TemporaryDirectory() as tmpdir:
        repo_dir = os.path.join(tmpdir, "vidjil")
        clone_vidjil_repo(repo_dir)
        make_vidjil_germline(repo_dir)

        germline_dir = os.path.join(repo_dir, "germline")
        merge_presets(germline_dir)
        create_vidjil_archive(germline_dir, args.out_archive)

    logger.info("Vidjil reference generation completed successfully.")


if __name__ == "__main__":
    main()
