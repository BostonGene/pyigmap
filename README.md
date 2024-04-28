# pyigmap

`PyIgMap` is a [Nextflow](https://github.com/nextflow-io/nextflow)-driven and Python-based workflow for mapping and annotating (in [AIRR](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields) format) TCR/BCR repertoire sequencing data. 

## Quick start

1. This pipeline requires Docker, Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html).
2. Download repository:
```bash
git clone https://github.com/BostonGene/pyigmap.git
```
3. Run this command to install nextflow and build container steps:
```bash
cd pyigmap
make # sudo apt install make
```
4. Start running your own analysis!

For amplicon data:
```bash
./nextflow run main.nf --mode amplicon --outdir ./results --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz
```

For RNASeq data:
```bash
./nextflow run main.nf --mode rnaseq --outdir ./results --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz
```

For public data:
```bash
./nextflow run main.nf --mode rnaseq --outdir ./results --sample SRR3743469 --reads 10000
```