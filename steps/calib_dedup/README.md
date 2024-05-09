# calib-dedup step

This step is a wrapping of [calib](https://github.com/vpc-ccg/calib) alignment-free UMI-based deduplication tool.

## Parameters

### UMI-extraction params

Works something like [mixcr barcode pattern extraction](https://mixcr.com/mixcr/reference/ref-tag-pattern/).

* `--fq1-barcode-pattern` — barcode pattern of the fq1
* `--fq2-barcode-pattern` — barcode pattern of the fq2
* `--pattern-max-error-budget` — maximum error size (budget) between pattern and read substring
* `--find-umi-in-reverse-complement` — enable finding umi in reverse complement reads

### UMI-clustering params

See [here](https://github.com/vpc-ccg/calib?tab=readme-ov-file#clustering-parameters) for more information.

* `--kmer-size`
* `--minimizer-count`
* `--minimizer-threshold`
* `--error-tolerance`

### UMI error correction params

See [here](https://github.com/vpc-ccg/calib?tab=readme-ov-file#error-correction-parameters) for more information

* `--min-reads-per-cluster` — the minimum number of reads required in a cluster to output the cluster consensus; default is 1.
* `--max-reads-per-cluster` — the maximum number of reads required in a cluster to output the cluster consensus; default is 1000.

## Input

* `--in-fq1`: path to the forward fastq (`path/to/R1.fastq.gz`)
* `--in-fq2`: path to the reverse fastq (`path/to/R2.fastq.gz`)

## Output

* `--out-fq1`: path to the output deduplicated forward fastq (`path/to/cR1.fastq.gz`)
* `--out-fq2`: path to the output deduplicated reverse fastq (`path/to/cR2.fastq.gz`)
* `--out-json`: path to the output json with metrics (`path/to/calib.json`)

## How to run

```bash
docker build --target tool -t calib-dedup .

FOLDER_WITH_DATA=path/to/your/folder # should contain: R1.fastq.gz and R2.fastq.gz

docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   calib-dedup \
   --in-fq1 /root/R1.fastq.gz \
   --in-fq2 /root/R2.fastq.gz \
   --out-fq1 /root/cR1.fastq.gz \
   --out-fq2 /root/cR2.fastq.gz \
   --out-json /root/calib.json \
   --kmer-size 8 \
   --minimizer-count 7 \
   --minimizer-threshold 7 \
   --error-tolerance 2 \
   --fq1-barcode-pattern "^UMI:N{12}" # UMI - the first 12 nucleotides
```