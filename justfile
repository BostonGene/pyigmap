# justfile
# -----------------------------------------------------------------------------
# Global settings for just
# - Run all recipes in bash with strict flags (-e exit on error, -u unset vars, -o pipefail)
# - Load environment variables from .env automatically

set shell := ["bash", "-eu", "-o", "pipefail", "-c"]
set dotenv-load := true

# ---- System info

SYSTEM_OS := `uname -s`
SYSTEM_ARCH := `uname -m`

# ---- Versions / flags

JAVA_VERSION := env("JAVA_VERSION", "21")
PYTHON_VERSION := env("PYTHON_VERSION", "3.12")
UV_VERSION := env("UV_VERSION", "latest")

# Container engine: docker or podman

ENGINE := env("ENGINE", "docker")

# Podman-specific flags

PODMAN_PARAM := if ENGINE == "podman" { "--format docker" } else { "" }

# Reference builder stage name

BUILD_REF_STAGE := "build-ref"

# ---- Binaries

ruff := "uv run ruff"
pytest := "uv run pytest"
pyrefly := "uv run pyrefly"
nextflow := "uv run nextflow"

# ---- Install uv into system
[group('setup')]
install-uv:
    curl -LsSf https://astral.sh/uv/install.sh | sh

# ---- Create or refresh virtualenv in .venv
[group('setup')]
create-env:
    uv venv --seed --python {{ PYTHON_VERSION }}

# ---- Install all project dependencies
[group('setup')]
sync:
    uv sync

# ---- Update project dependencies
[group('setup')]
update:
    @echo "Updating packages..."
    uv venv --python {{ PYTHON_VERSION }}
    uv add "nf-core>=2.14.1" "nextflow>=24.04.2"

# ---- Update development dependencies
[group('setup')]
update-dev:
    @echo "Updating dev packages..."
    uv venv --python {{ PYTHON_VERSION }}
    uv sync --group dev
    @echo "Installing component packages in editable mode..."
    uv pip install -e bin/pyumi
    uv pip install -e bin/calib_dedup
    uv pip install -e bin/reporter
    uv pip install -e bin/fastp
    uv pip install -e bin/vidjil
    uv pip install -e bin/igblast
    uv pip install -e bin/cdr3nt-error-corrector
    uv pip install -e bin/fq-downloader

# ---- Build reference image for a specific step
[group('build')]
build-ref-image STEP:
    @echo "Building reference image for {{ STEP }}..."
    {{ ENGINE }} build --target {{ BUILD_REF_STAGE }} -t {{ STEP }}-{{ BUILD_REF_STAGE }} {{ PODMAN_PARAM }} bin/{{ STEP }}/

# ---- Build IgBLAST reference with major alleles only (*01)
[group('build')]
build-igblast-ref-major:
    @echo "Building VDJ reference with major alleles (*01) for IgBLAST..."
    {{ ENGINE }} run --rm -v $(pwd)/bin/igblast:/tmp igblast-{{ BUILD_REF_STAGE }} -a -o /tmp/igblast.reference.major_allele.tar.gz

# ---- Build IgBLAST reference with all alleles (*01, *02, etc.)
[group('build')]
build-igblast-ref-all:
    @echo "Building VDJ reference with all alleles (*01, *02, etc.) for IgBLAST..."
    {{ ENGINE }} run --rm -v $(pwd)/bin/igblast:/tmp igblast-{{ BUILD_REF_STAGE }} -a -o /tmp/igblast.reference.all_alleles.tar.gz

# ---- Build Vidjil germline reference
[group('build')]
build-vidjil-ref:
    @echo "Building VDJ reference for Vidjil..."
    {{ ENGINE }} run --rm -v $(pwd)/bin/vidjil:/tmp vidjil-{{ BUILD_REF_STAGE }} -a -o /tmp/vidjil.germline.tar.gz

# ---- Build OLGA models
[group('build')]
build-olga-models:
    @echo "Building OLGA models for cdr3nt-error-corrector..."
    {{ ENGINE }} run --rm -v $(pwd)/bin/cdr3nt_error_corrector:/tmp cdr3nt_error_corrector-{{ BUILD_REF_STAGE }} -a -o /tmp/olga-models.tar.gz

# ---- Build all references (IgBLAST, Vidjil, OLGA)
[group('build')]
build-ref:
    @echo "Building all references..."
    just build-ref-image igblast
    just build-igblast-ref-major
    just build-igblast-ref-all
    just build-ref-image vidjil
    just build-vidjil-ref
    just build-ref-image cdr3nt_error_corrector
    just build-olga-models

# ---- Build step image (tool or image stage)
[group('build')]
build-step-image STEP STAGE:
    @{{ ENGINE }} version
    {{ ENGINE }} build --target {{ STAGE }} -t {{ STEP }}-{{ STAGE }} {{ PODMAN_PARAM }} bin/{{ STEP }}

# ---- Build all development images
[group('build')]
dev:
    @echo "Building development environment..."
    just clean
    just build-step-image fq-downloader image
    #!/usr/bin/env bash
    set -euo pipefail
    for step in pyumi calib_dedup reporter fastp vidjil igblast cdr3nt_error_corrector; do
    just build-step-image $step image
    just build-step-image $step tool
    done
    just update-dev

# ---- Clean all generated files, images, and caches
[group('quality')]
clean:
    @echo "Removing temporary files and images..."
    -{{ ENGINE }} rmi -f fq-downloader-image \
        pyumi-tool pyumi-image \
        calib_dedup-tool calib_dedup-image \
        reporter-tool reporter-image \
        fastp-tool fastp-image \
        vidjil-tool vidjil-image \
        igblast-tool igblast-image \
        cdr3nt_error_corrector-tool cdr3nt_error_corrector-image
    @rm -rf .venv \
        .nextflow.log* work .nextflow nextflow uv.lock \
        bin/igblast/igblast.reference.all_alleles.tar.gz \
        bin/igblast/igblast.reference.major_allele.tar.gz \
        bin/igblast/vidjil.germline.tar.gz \
        bin/cdr3nt_error_corrector/olga-models.tar.gz
    @rm -rf .pytest_cache .ruff_cache .pyrefly_cache .coverage coverage.xml htmlcov
    @find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
    @find . -type f -name '*.pyc' -delete 2>/dev/null || true

# ---- Run ruff linter
[group('quality')]
lint:
    @echo "Running ruff check..."
    {{ ruff }} check bin/calib_dedup/ bin/fastp/ bin/vidjil/ bin/igblast/ bin/cdr3nt-error-corrector/ bin/pyumi/ bin/reporter/ bin/fq-downloader/

# ---- Auto-format code with ruff
[group('quality')]
format:
    @echo "Running ruff format..."
    {{ ruff }} format bin/calib_dedup/ bin/fastp/ bin/vidjil/ bin/igblast/ bin/cdr3nt-error-corrector/ bin/pyumi/ bin/reporter/ bin/fq-downloader/

# ---- Run pyrefly type checker
[group('quality')]
type:
    @echo "Running pyrefly type checker..."
    {{ pyrefly }} check bin/calib_dedup/
    {{ pyrefly }} check bin/fastp/
    {{ pyrefly }} check bin/vidjil/
    {{ pyrefly }} check bin/igblast/
    {{ pyrefly }} check bin/cdr3nt-error-corrector/
    {{ pyrefly }} check bin/pyumi/
    {{ pyrefly }} check bin/reporter/
    {{ pyrefly }} check bin/fq-downloader/

# ---- Run unit tests (Python function-level tests)
[group('quality')]
tests-unit:
    @echo "Running Python unit tests..."
    {{ pytest }} bin/pyumi/unit_tests -vv
    {{ pytest }} bin/cdr3nt-error-corrector/unit_tests -vv

# ---- Run component tests (individual tool/container tests)
[group('quality')]
tests-component:
    @echo "Running component tests..."
    {{ pytest }} tests -vv --kwdof --tag component-tests

# ---- Run integration tests (full workflow tests)
[group('quality')]
tests-integration:
    @echo "Running integration tests..."
    {{ pytest }} tests -vv --kwdof --tag integration-tests

# ---- Run all tests (unit + component + integration)
[group('quality')]
tests:
    @echo "Running all tests..."
    just tests-unit
    just tests-component
    just tests-integration

# ---- Run full quality suite: format + lint + type check + tests
[group('quality')]
check:
    just format
    just lint
    just type
    just tests
    @echo "✅ All checks passed"

# ---- Install and setup everything (default target)
[group('setup')]
install:
    just update
    just build-ref
    just build-step-image fq-downloader image
    #!/usr/bin/env bash
    set -euo pipefail
    for step in pyumi calib_dedup reporter fastp vidjil igblast cdr3nt_error_corrector; do
    just build-step-image $step image
    done

# ---- Install Docker (Ubuntu/Debian)
[group('ci-cd')]
install-docker:
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    sudo add-apt-repository "deb [arch={{ SYSTEM_ARCH }}] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
    sudo apt-get update
    sudo apt-get -y -o Dpkg::Options::="--force-confnew" install docker-ce

# ---- Install Podman (Ubuntu/Debian)
[group('ci-cd')]
install-podman:
    sudo apt-get update
    sudo apt-get -y install podman

# ---- Install Java
[group('ci-cd')]
install-java:
    curl -fL https://download.oracle.com/java/{{ JAVA_VERSION }}/latest/jdk-{{ JAVA_VERSION }}_linux-x64_bin.tar.gz -o /tmp/java.tar.gz
    sudo mkdir -p /usr/local/bin/java-{{ JAVA_VERSION }}
    sudo tar --strip-components=1 -xzvf /tmp/java.tar.gz -C /usr/local/bin/java-{{ JAVA_VERSION }}
    rm /tmp/java.tar.gz
    @echo "Add to your shell config:"
    @echo "export JAVA_HOME=/usr/local/bin/java-{{ JAVA_VERSION }}"
    @echo 'export PATH=$$JAVA_HOME/bin:$$PATH'

# ---- Install gh CLI
[group('ci-cd')]
install-gh:
    curl -sS https://webi.sh/gh | sh

# ---- Start Docker
[group('ci-cd')]
docker-up:
    @if command -v brew >/dev/null; then \
        open -a Docker; \
    else \
        sudo systemctl start docker; \
    fi

# ---- Stop Docker
[group('ci-cd')]
docker-down:
    @if command -v brew >/dev/null; then \
        osascript -e 'quit app "Docker"'; \
    else \
        sudo systemctl stop docker; \
    fi

# ---- Install git-cliff for changelog generation
[group('release')]
install-cliff:
    @echo "Installing git-cliff..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    ~/.cargo/bin/cargo install --version 2.6.1 git-cliff

# ---- Generate changelog with git-cliff
[group('release')]
changelog TAG:
    @echo "Updating CHANGELOG.md for {{ TAG }}..."
    git-cliff --tag {{ TAG }} --bump -o CHANGELOG.md
    git add CHANGELOG.md && git commit -m "chore(release): prepare for {{ TAG }}"
    git show

# ---- Create and push release tag
[group('release')]
release TAG:
    just install-gh
    just changelog {{ TAG }}
    git tag {{ TAG }} && git push origin {{ TAG }}
    git push
    git cliff --strip all --latest | gh release create {{ TAG }} --notes-file -

# ---- Show help for all recipes
help:
    @just --list --unsorted
