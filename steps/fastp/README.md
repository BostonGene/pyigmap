# Fastp step

This step is a wrapping of [fastp](https://github.com/OpenGene/fastp) tool designed to provide fast all-in-one preprocessing for FastQ files.

## Overview

By default:
* Trims Illumina adapters
* Trims polyG for 2ch SBS
* Filters out reads if >40% bases have quality < 20
* Filters out reads shorter than 15bp
* Filters out reads with more than 5 N's
* Generates JSON & HTML report

## Parameters:
* `--merge`: merges forward and reverse fastq.
* `--disable`: disables modes. Example:
    ```bash
    --disable "length_filtering" # disables filtering reads shorter than 15bp
    --disable "quality_filtering" # disables quality filtering (if >40% bases have quality <20)
    --disable "adapter_trimming" # disables trimming Illumina adapters
    ```
* `--mock-merge`: enables mock merging of not overlapped forward and reverse reads with a selected insert size (distance)
    * `--insert-size`: insert distance between reads in mock merge 

Example:

Not overlapped forward and reverse reads:
```
CCCAAA ->
          <- GGGTTT
```
Mock merging with `--insert-size 1`:
```
CCCAAANAAACCC
```

## Input

  * `--in-fq1`: path to the forward fastq (`path/to/R1.fastq.gz`)
  * `--in-fq2` (**optional**): path to the reverse fastq (`path/to/R2.fastq.gz`)

## Output

  * `--out-fq1` (**optional**): output forward fastq
  * `--out-fq2` (**optional**): output reverse fastq
  * `--out-fq12` (**optional**): output merged fastq
  * `--out-json`: output json with metrics
  * `--out-html`: output html with metrics

## How to run

```bash
docker build --target tool -t fastp .

FOLDER_WITH_DATA=path/to/your/folder # should contain: R1.fastq.gz (or cR1.fastq.gz) and R2.fastq.gz (or cR2.fastq.gz)

docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   fastp \
   --in-fq1 /root/R1.fastq.gz \
   --in-fq2 /root/R2.fastq.gz \
   --out-fq1 /root/mR1.fastq.gz \
   --out-fq2 /root/mR2.fastq.gz \
   --out-fq12 /root/mR12.fastq.gz \
   --disable "length_filtering" "adapter_trimming" "quality_filtering" \
   --merge \
   --out-html /root/fastp.html \
   --out-json /root/fastp.json
```
