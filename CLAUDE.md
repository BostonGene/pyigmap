# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

pyIgMap is a Nextflow DSL2 pipeline for extracting and summarizing antigen receptor gene rearrangements (BCR/TCR) from sequencing data. It processes both bulk RNASeq and targeted AIRR-Seq protocols, outputting results in AIRR-compliant format.

The pipeline combines multiple bioinformatics tools (Fastp, Vidjil, IgBLAST, OLGA) wrapped in Python scripts, executed through Nextflow workflows, and packaged in Docker/Podman containers.

## Build and Setup Commands

```bash
# Initial setup (installs Nextflow, builds references and images)
just install

# For Podman instead of Docker
ENGINE=podman just install

# Development setup (builds dev images, creates venv, installs dependencies)
# This also runs `uv sync --group dev` to install dev dependencies
just dev

# Sync dependencies after pulling changes
just sync           # Production dependencies
just update-dev     # Development dependencies (includes uv sync --group dev)

# Build only V(D)J reference archives
just build-ref

# Clean all generated files and images
just clean
```

## Testing Commands

The project has three levels of tests:

```bash
# Run all tests (unit + component + integration)
just tests

# Run unit tests only (Python function-level tests)
just tests-unit

# Run component tests only (individual tool/container tests)
just tests-component

# Run integration tests only (full pipeline end-to-end tests)
just tests-integration

# Code quality checks
just check       # full suite: format + lint + type check + tests
just lint        # ruff linter
just format      # ruff formatter
just check-types # pyrefly type checking
```

### Test Levels Explained

1. **Unit Tests** (`just tests-unit`): Python function/class tests in `bin/*/unit_tests/`
   - Test individual Python functions in isolation
   - Fast, no Docker required
   - Example: `bin/pyumi/unit_tests/test_pattern.py`

2. **Component Tests** (`just tests-component`): Individual tool/container tests in `tests/*/`
   - Test single Docker containers/tools with real data
   - Tagged with `component-tests` in pytest-workflow YAML files
   - Example: `tests/fastp/test_fastp.yml` (runs `docker run fastp-tool`)

3. **Integration Tests** (`just tests-integration`): Full pipeline tests in `tests/test_main.yml`
   - Test complete Nextflow workflows end-to-end
   - Tagged with `integration-tests` in pytest-workflow YAML files
   - Example: `tests/test_main.yml` (runs `nextflow run main.nf`)

## Running the Pipeline

```bash
# Basic RNASeq-bulk run
uv run nextflow main.nf -profile docker \
    --library rnaseq \
    --fq1 R1.fastq.gz \
    --fq2 R2.fastq.gz \
    --outdir ./results

# AIRR-Seq targeted with UMI
uv run nextflow main.nf -profile docker \
    --library amplicon \
    --fq1 R1.fastq.gz \
    --fq2 R2.fastq.gz \
    --fq2_pattern "^ADAPTER(UMI:N{19})ADAPTER" \
    --outdir ./results

# From public database (SRA, GEO, etc.)
uv run nextflow main.nf -profile docker \
    --library rnaseq \
    --sample_id SRR3743469 \
    --outdir ./results

# Resume a failed run
uv run nextflow main.nf -resume -profile docker [other params...]
```

## Architecture

### Pipeline Flow

There are two main workflows based on library type:

1. **PYIGMAP_AMPLICON** (targeted AIRR-Seq):
   - Optional UMI extraction (PyUMI) → clustering/deduplication (CalibDedup) → read merging (Fastp) → V(D)J mapping (IgBLAST) → error correction (CDR3ErrorCorrector)
   - Without UMI: directly to read merging → V(D)J mapping → error correction

2. **PYIGMAP_RNASEQ** (bulk RNASeq):
   - Read merging (Fastp) → junction detection (Vidjil) → V(D)J mapping (IgBLAST on FASTA) → error correction with pgen filtering (CDR3ErrorCorrector + OLGA)

### Directory Structure

- `main.nf` - Main Nextflow entrypoint, routes to appropriate subworkflow based on `--library` parameter
- `subworkflows/` - Reusable workflow definitions (pyigmap_amplicon.nf, pyigmap_rnaseq.nf, download_fastq.nf)
- `modules/local/` - Individual Nextflow process definitions, one per tool
- `bin/` - Python tools and their Dockerfiles:
  - `pyumi/` - UMI extraction using regex patterns
  - `calib_dedup/` - UMI-based clustering and consensus generation
  - `fastp/` - Read merging and quality control
  - `vidjil/` - Seed-based V(D)J junction detection
  - `igblast/` - V(D)J mapping against IMGT reference
  - `cdr3nt_error_corrector/` - Clonotype filtering, aggregation, pgen calculation
  - `reporter/` - UMI metrics visualization
  - `fq-downloader/` - Public database fetching
- `conf/` - Nextflow configuration (base.config for resource limits, modules.config for per-process settings)
- `tests/` - Pytest-workflow integration tests

### Key Python Components

Each tool in `bin/*/` follows this pattern:
- `src/{component}/` - Python source code following standard src layout
  - `run.py` - CLI entrypoint with main() function
  - `utils.py` - Core logic functions (if needed)
  - `logger.py` - Logging setup
  - `__init__.py` - Package initialization
- `Dockerfile` - Multi-stage build using UV for dependency management
- `pyproject.toml` - Python dependencies, metadata, and build config (hatchling for most, maturin for fastp)
- `uv.lock` - Locked dependency versions (generated by uv)
- `unit_tests/` - Pytest unit tests (if applicable)

**Package Structure:**
- All components use absolute imports (e.g., `from pyumi.logger import set_logger`)
- Entry points defined in pyproject.toml `[project.scripts]` section
- Source code organized in `src/{component}/` directories

**Dependency Management:**
- All Python tools use `pyproject.toml` + `uv.lock` for reproducible builds
- Requirements Python version: `>=3.12`
- To add dependencies: edit `pyproject.toml`, then run `uv lock` in the tool directory
- Tools without Python dependencies (binary wrappers) still have pyproject.toml for packaging

**Container Architecture:**
- Base images: `python:3.12-slim-bullseye` (Debian 11) for most tools
- Special cases:
  - `fq-downloader`: Uses `python:3.12-slim-bookworm` (Debian 12) for Ubuntu 24.04-like features
  - `fastp`: Uses `ubuntu:24.04` base with Rust toolchain for native extension builds via maturin
- Dependency installation: `uv sync --frozen --no-dev` using pyproject.toml + uv.lock
- UV environment variables:
  - `UV_SYSTEM_PYTHON=1` - Use system Python, don't download/manage Python versions
  - `UV_COMPILE_BYTECODE=1` - Pre-compile .pyc files for faster imports
  - `UV_PYTHON_DOWNLOADS=never` - Never auto-download Python (use base image Python)
- Entry point: `python3.12 -m {component}.run` for all components

### Reference Data

The pipeline requires three reference archives built at setup:
- `igblast.reference.major_allele.tar.gz` or `igblast.reference.all_alleles.tar.gz` - IMGT V(D)J genes for IgBLAST
- `vidjil.germline.tar.gz` - V(D)J germline sequences for Vidjil
- `olga-models.tar.gz` - OLGA models for pgen calculation (RNASeq only)

These are built via `just build-ref` which creates temporary Docker containers to generate and package the references.

### Configuration

- Pipeline parameters are set via CLI (`--param value`) or `-params-file`, NOT via `-c` custom configs
- Use `-profile docker` or `-profile podman` to select container engine
- `--library amplicon` or `--library rnaseq` determines which workflow executes
- `--all_alleles` flag switches between major-allele-only (*01) or all-alleles (*01, *02, etc.) reference
- Resource limits: `params.max_cpus`, `params.max_memory`, `params.max_time` (see nextflow.config)

## Common Development Tasks

### Adding a New Python Tool Step

1. Create `bin/newtool/` directory structure:
   ```
   bin/newtool/
   ├── src/
   │   └── newtool/
   │       ├── __init__.py
   │       ├── run.py (with main() function)
   │       ├── utils.py (optional)
   │       └── logger.py (optional)
   ├── unit_tests/ (optional)
   ├── pyproject.toml
   ├── Dockerfile
   └── README.md (optional)
   ```

2. Create `pyproject.toml` with:
   - `[build-system]` section (hatchling for most tools, maturin for Rust extensions)
   - `[project]` section with `requires-python = ">=3.12"` and dependencies
   - `[project.scripts]` section defining CLI entry point: `newtool = "newtool.run:main"`
   - `[tool.hatch.build.targets.wheel]` with `packages = ["src/newtool"]` (for hatchling)
   - For maturin projects: `[tool.maturin]` with `python-source = "src"` and `module-name = "newtool"`

3. Create `Dockerfile` using UV and Python 3.12-slim:
   ```dockerfile
   FROM python:3.12-slim-bullseye AS image

   # Install uv
   COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

   ENV UV_SYSTEM_PYTHON=1
   ENV UV_COMPILE_BYTECODE=1
   ENV UV_PYTHON_DOWNLOADS=never

   WORKDIR /app
   COPY pyproject.toml uv.lock ./
   RUN uv sync --frozen --no-dev

   COPY src/ ./src/

   FROM image AS tool
   ENTRYPOINT ["python3.12", "-m", "newtool.run"]
   ```

4. Run `uv lock` in the tool directory to generate `uv.lock`
5. Add Nextflow process in `modules/local/newtool.nf` with appropriate container and resource labels
6. Integrate process into relevant subworkflow in `subworkflows/`
7. Add unit tests in `bin/newtool/unit_tests/` (use absolute imports: `from newtool.module import ...`)
8. Add workflow test in `tests/newtool/test_newtool.yml`
9. Build images: `just build-step-image newtool image && just build-step-image newtool tool`

### Running a Single Test

```bash
# Single component test (individual tool/container)
uv run pytest tests/fastp/test_fastp.yml -vv

# Single Python unit test file
uv run pytest bin/pyumi/unit_tests/test_pattern.py -vv

# Specific Python unit test function
uv run pytest bin/pyumi/unit_tests/test_pattern.py::test_function_name -vv

# Single integration test (full workflow)
uv run pytest tests/test_main.yml::test_main_rnaseq_sample_id -vv
```

### Modifying Container Images

After changing Python code or Dockerfile in `bin/STEP/`:
```bash
# If you modified pyproject.toml dependencies
cd bin/STEP && uv lock

# Rebuild specific tool image
just build-step-image igblast image

# Rebuild tool image
just build-step-image igblast tool

# Or rebuild everything
just dev
```

### Working with References

```bash
# Rebuild IgBLAST references
just build-ref-image igblast
just build-igblast-ref-major  # or build-igblast-ref-all

# Rebuild Vidjil reference
just build-ref-image vidjil
just build-vidjil-ref

# Rebuild OLGA models
just build-ref-image cdr3nt_error_corrector
just build-olga-models
```

## Important Patterns

### Error Handling in Python Tools

All tools use consistent patterns:
- Argparse for CLI with `check_argument_consistency()` validation
- Logger setup via `logger.py` with module-level logging
- Structured error messages via `exit_with_error()` or similar
- Return proper exit codes

### Nextflow Process Conventions

- Use `label` directives for resource classes: `process_single`, `process_low`, `process_medium`, `process_high`
- Container names follow pattern: `TOOLNAME-image` (e.g., `igblast-image`)
- Save important intermediate files when `params.save_all` is true
- Emit named outputs for channel passing

### Testing Patterns

- Unit tests: pytest with fixtures, parametrized tests, clear arrange-act-assert structure
- Integration tests: pytest-workflow with YAML test definitions containing `command`, `files`, `exit_code` checks
- Tag tests: `--tag unit-tests` or `--tag integration-tests` for selective running

## Container Engine Notes

- Pipeline supports both Docker and Podman via `ENGINE` variable
- Images use multi-stage builds: `build-ref` stage for reference generation, `image` stage for runtime tools
- Podman requires `--format docker` flag added automatically via `PODMAN_PARAM` variable
- All containers are rootless-compatible

## Output Format

Final output is `pyigmap.tar.gz` containing:
- `corrected_annotation.tsv` - AIRR-formatted clonotype table with V/D/J calls, CDR3 sequences, duplicate counts
- `stat.json` - Summary metrics (total reads, productive clonotypes, etc.)
- Intermediate files if `--save_all` is enabled
