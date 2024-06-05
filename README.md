# pyigmap

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11103554.svg)](https://doi.org/10.5281/zenodo.11103554)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with podman](https://img.shields.io/badge/run%20with-podman-892CA0?labelColo=000000&logo=podman)](https://podman.io/)

## Pipeline summary

`pyIgMap` is a [Nextflow](https://github.com/nextflow-io/nextflow)-driven and Python-based workflow for extracting and summarizing antigen receptor gene rearrangements from sequencing data.

<p align="center">
    <img title="Pyigmap Workflow" src="docs/images/pyigmap_workflow.svg">
</p>

## Quick start

1. Clone a repository:

```bash
git clone https://github.com/BostonGene/pyigmap.git
cd pyigmap
```

2. This workflow requires [Docker](https://docs.docker.com/engine/install/) or [Podman](https://podman.io/), Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html). You can install it manually or execute:

```bash
make install-docker # requirements: Ubuntu x64
make install-podman # requirements: Ubuntu 20.10 and newer
make install-java # requirements: Linux x64
```

3. Run this command to install all dependencies:

```bash
make # will use Docker as container engine
```

or

```bash
make ENGINE=podman # will use Podman as container engine
```

4. Start running your own analysis!

```bash
# for amplicon data
./pyigmap --mode amplicon --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz -profile docker

# for RNASeq data
./pyigmap --mode rnaseq --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz -profile docker

# for public data by sample id
./pyigmap --mode rnaseq --sample SRR3743469 --reads 200000 -profile docker

# for public data from ZENODO
./pyigmap --mode rnaseq --zenodo --fq1 SRR3743469_R1.fastq.gz --fq2 SRR3743469_R2.fastq.gz --reads 200000 -profile docker
```

## Contributing

Contributions are more than welcome. See the [`CONTRIBUTING`](CONTRIBUTING.md) for details.
