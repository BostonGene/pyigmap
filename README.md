# pyigmap

`pyIgMap` is a [Nextflow](https://github.com/nextflow-io/nextflow)-driven and Python-based workflow for extracting and summarizing antigen receptor gene rearrangements from sequencing data.

## Quick start

1. This pipeline requires [Docker](https://docs.docker.com/engine/install/), Bash 3.2 (or later) and [Java 11 (or later, up to 22)](http://www.oracle.com/technetwork/java/javase/downloads/index.html).

```bash
sudo make install-docker # requires linux ubuntu x64
sudo make install-java # requires linux x64
```

2. Clone repository:

```bash
git clone https://github.com/BostonGene/pyigmap.git
```

3. Run this command to install nextflow and build container steps:

```bash
cd pyigmap
make # sudo apt install make
chmod +x pyigmap
```

4. Start running your own analysis!

```bash
# for amplicon data
./pyigmap --mode amplicon --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz

# for RNASeq data
./pyigmap --mode rnaseq --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz

# for public data by sample id
./pyigmap --mode rnaseq --sample SRR3743469 --reads 200000

# for public data from ZENODO
./pyigmap --mode rnaseq --zenodo --fq1 SRR3743469_R1.fastq.gz --fq2 SRR3743469_R2.fastq.gz --reads 200000
```

## Benchmark datasets

[Datasets](https://zenodo.org/records/11103555) for benchmarking a `pyIgMap` pipeline.

## Contributing

Contributions are more than welcome. See the [CONTRIBUTING.md](CONTRIBUTING.md) file for details.
