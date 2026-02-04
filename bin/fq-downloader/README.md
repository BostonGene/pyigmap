# fq-downloader component

> Owner: nikita.syzrantsev@bostongene.com

This tool retrieves experiment metadata using [ffq](https://github.com/pachterlab/ffq), downloads the corresponding sequencing data, and merges the FASTQ files from multiple runs into combined forward and reverse FASTQ files.

## Input

* `--identifier`: experiment identifier (accession number).

## Output

* `--out-fq1`: basename of output forward FASTQ
* `--out-fq2`: basename of output reverse FASTQ

The resulting files will contain all forward and reverse reads from the experiment's runs, respectively, merged into two consolidated FASTQ files.
