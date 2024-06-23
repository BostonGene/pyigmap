# pyigmap

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11103554.svg)](https://doi.org/10.5281/zenodo.11103554)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
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

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Pipeline summary

`pyigmap` allows the processing of raw BCR/TCR sequencing data from **bulk** and **targeted** sequencing protocols.
For more details on the supported protocols, please refer to the [usage](#Usage) documentation.

### 1. FASTQ pre-processing

* RNASeq-bulk:
  * Merging overlapping reads, joining non-overlapping reads with a selected insert size, and raw read quality control (`Fastp`).

* AIRR-Seq (target):
  * Extracting the UMI from the reads (`PyUMI`).
  * Alignment-free clustering of UMI tagged reads with subsequent consensus generation (`Calib`).
  * Merging overlapping reads, saving not-overlapping reads, and raw read quality control (`Fastp`).

### 2. V(D)J mapping

* RNASeq-bulk:
  * Fast detecting V(D)J junctions from FASTQ reads using a seed-based heuristic without initial alignment to database germline sequences (`Vidjil`).
  * Mapping previously identified junctions against [IMGT reference](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/) and producing AIRR-formatted table (`IgBLAST`).

* AIRR-Seq (target):
  * Mapping FASTQ reads against [IMGT reference](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/) and producing AIRR-formatted table (`IgBLAST`).

### 3. Cleansing and aggregating clonotypes

* RNASeq-bulk and AIRR-Seq (target):
  * Filter out chimeric clonotypes, that have different locus in V-/J-segments (except for TRA and TRD).
  * Store the best aligned V, D, J and C genes call.
  * Discard clonotypes with undefined nucleotide or amino acid in CDR3 sequence.
  * Aggregate clonotypes based on [Levenstein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of 1 and read count ratio and subsequent `duplicate_count` column calculation.

* Only RNASeq-bulk:
  * Compute generation probabilities (`pgen` column) of CDR3 amino acid sequences and remove clonotypes with `pgen` values less than or equal to the selected threshold (`OLGA`).
  * Store clonotypes with the most weighted and frequent C-gene call and V-gene alignment call.

## Usage

A typical command to run the pipeline from **RNASeq-bulk** sequencing data is:

```bash
./pyigmap -profile <docker/podman> \
    --library rnaseq \
    --fq1 "R1.fastq.gz" \
    --fq2 "R2.fastq.gz" \
    --outdir "./results"
```

For common **AIRR-Seq targeted** sequencing protocols we provide pre-set parameters, including a parameter for specifying a UMI barcode pattern.  

Here is an example command to process the data from the **AIRR-Seq targeted** protocol, where there is a 19-base pair UMI located between two adapters in the reverse FASTQ file:

```bash
./pyigmap -profile <docker/podman> \
    --library amplicon \
    --fq1 "R1.fastq.gz" \
    --fq2 "R2.fastq.gz" \
    --fq2_pattern "^TGGTATCAACGCAGAGTAC(UMI:N{19})TCTTGGGGG" \
    --outdir "./results"
```

You can also use public data from these databases by using a sample ID: [GEO](https://www.ncbi.nlm.nih.gov/geo/), [SRA](https://www.ncbi.nlm.nih.gov/sra), [EMBL-EBI](https://www.ebi.ac.uk/), [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html), [NIH Biosample](https://www.ncbi.nlm.nih.gov/biosample) and [ENCODE](https://www.encodeproject.org/):

```bash
./pyigmap -profile <docker/podman> \
    --library rnaseq \
    --sample_id SRR3743469 \
    --outdir "./results"
```

Alternatively, you can provide an HTTP/HTTPS/FTP link to your FASTQ files.

```bash
./pyigmap -profile <docker/podman> \
    --library amplicon \
    --fq1 https://zenodo.org/records/11103555/files/fmba_TRAB_R1.fastq.gz \
    --fq2 https://zenodo.org/records/11103555/files/fmba_TRAB_R2.fastq.gz \
    --outdir "./results"
```

## Contributing

Contributions are more than welcome. See the [`CONTRIBUTING`](CONTRIBUTING.md) for details.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
