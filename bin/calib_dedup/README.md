# calib-dedup step

This step is a wrapping of [Calib](https://github.com/vpc-ccg/calib) alignment-free UMI-based deduplication tool.

## Parameters

### UMI-clustering params

See [here](https://github.com/vpc-ccg/calib?tab=readme-ov-file#clustering-parameters) for more information.

* `--fq1-umi-length` — length of UMI barcode in forward FASTQ
* `--fq2-umi-length` — length of UMI barcode in reverse FASTQ
* `--kmer-size`
* `--minimizer-count`
* `--minimizer-threshold`
* `--error-tolerance` — a [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) between two barcodes

### UMI error correction params

See [here](https://github.com/vpc-ccg/calib?tab=readme-ov-file#error-correction-parameters) for more information

* `--min-reads-per-cluster` — the minimum number of reads required in a cluster to output the cluster consensus; default is 1.
* `--max-reads-per-cluster` — the maximum number of reads required in a cluster to output the cluster consensus; default is 1000.

## Input

* `--in-fq1`: path to the forward FASTQ (`path/to/R1.fastq.gz`)
* `--in-fq2`: path to the reverse FASTQ (`path/to/R2.fastq.gz`)

## Output

* `--out-fq1`: path to the output deduplicated forward FASTQ (`path/to/cR1.fastq.gz`)
* `--out-fq2`: path to the output deduplicated reverse FASTQ (`path/to/cR2.fastq.gz`)

## How to run

```bash
docker build --target tool -t calib-dedup .

docker run \
   -v ./unit_tests/test_data:/root/ \
   calib-dedup \
   --in-fq1 /root/age_ig_s1_R1_umi.fastq.gz \
   --in-fq2 /root/age_ig_s1_R2_rc_umi.fastq.gz \
   --out-fq1 /root/cR1.fastq.gz \
   --out-fq2 /root/cR2.fastq.gz \
   --fq1-umi-length 13
```