# Fastp component

> Owner: nikita.syzrantsev@bostongene.com

This component is a wrapping of [fastp](https://github.com/OpenGene/fastp) tool designed to provide fast all-in-one preprocessing for FASTQ files with additional logic to merge non overlapped reads.

## Input

  * `--in-fq1`: Path to the forward FASTQ.
  * `--in-fq2` (**optional**): Path to the reverse FASTQ.
  * `--merge-reads`: Enable merging of forward and reverse FASTQ.
  * `--disable-filters`: List of modes that we need to disable.  
   Example:
      ```bash
      --disable-filters "length_filtering" # disables filtering reads shorter than 15bp
      --disable-filters "quality_filtering" # disables quality filtering (if >40% bases have quality <20)
      --disable-filters "adapter_trimming" # disables trimming Illumina adapters
      ```
  * `--mock-merge-reads`: Enable mock merging of not overlapped forward and reverse reads with a selected insert size (distance)
  * `--inner-distance-size`: Inner distance between reads (used in mock merging). Default: `1`.
  * `--reads-chunk-size`: The maximum number of reads that a chunk can contain to perform mock merging. Default: `5000000`.

Example of merging overlapped reads:

```
|------------------------| Full sequence
|----------ATG>            R1
          <ATG-----------| R2
|----------ATG-----------| Consensus (overlap)
```

Example of mock merging not-overlapped reads with `--innner-distance-size 3`:

```
|------------------------| Full sequence
|--------->                R1
              <----------| R2
|----------NNN-----------| Concatenated (no overlap)
```

## Output

  * `--out-fq1` (**optional**): Path to the output forward FASTQ. Default: `/outputs/R1.fastq.gz`.
  * `--out-fq2` (**optional**): Path to the output reverse FASTQ. Default: `/outputs/R2.fastq.gz`.
  * `--out-fq12` (**optional**): Path to the output merged FASTQ. Default: `/outputs/R12.fastq.gz`.
  * `--json`: Path to the output json with metrics. Default: `/outputs/fastp.json`.
  * `--html`: Path to the output html with metrics. Default: /outputs/fastp.html.