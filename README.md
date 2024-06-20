# pyigmap

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11103554.svg)](https://doi.org/10.5281/zenodo.11103554)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with podman](https://img.shields.io/badge/run%20with-podman-892CA0?labelColo=000000&logo=podman)](https://podman.io/)

## Introduction

`pyIgMap` is a [Nextflow](https://github.com/nextflow-io/nextflow)-driven and Python-based pipeline for extracting and summarizing antigen receptor gene rearrangements from sequencing data.

> [!NOTE]
> The pipeline is built upon open source components commonly used for AIRR-seq data processing. The output is provided in [AIRR format](https://docs.airr-community.org/en/stable/) enabling downstream analysis with AIRR-compliant software such as [Immcantation](https://immcantation.readthedocs.io/en/stable/index.html).

<p align="center">
    <img title="Pyigmap Workflow" src="docs/images/pyigmap_workflow.svg">
</p>

## Quick start

1. This pipeline requires [Docker](https://docs.docker.com/engine/install/) or [Podman](https://podman.io/), [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html), [Bash](https://www.gnu.org/software/bash/) and [Make](https://www.gnu.org/software/make/) tool.

> [!NOTE]
> To install the [Make](https://www.gnu.org/software/make/) tool execute: ```sudo apt install make```

2. Clone the repository and go inside:

```bash
git clone https://github.com/BostonGene/pyigmap.git
cd pyigmap
```
 
> [!TIP]
> If you have Ubuntu 20.10 (amd64) or higher, you can install requirements above using:  
> ```make install-docker``` or ```make install-podman```  
> ```make install-java```

3. Install Nextflow, build V(D)J references and Docker container images with a single command:

```bash
make
```
> [!TIP]
> To use a [Podman](https://podman.io/) as a container engine, run: ```make ENGINE=podman```

4. Start running your own analysis!

```bash
./pyigmap -profile <docker/podman> \
    --library <amplicon/rnaseq> \
    --fq1 "R1.fastq.gz" \
    --fq2 "R2.fastq.gz" \
    --outdir "./results"
```

## Pipeline summary

...

## Pipeline parameters

...

## Contributing

Contributions are more than welcome. See the [`CONTRIBUTING`](CONTRIBUTING.md) for details.
