# pyumi step

This step is a Python-based tool to preprocess UMI and CELL barcodes in FASTQ reads.

## Parameters

* `--fq1-pattern`: barcode pattern of the first FASTQ
* `--fq2-pattern`: barcode pattern of the second FASTQ
* `--max-error`: maximum error size (budget) between pattern and read substring
* `--find-in-reverse-complement`: enable finding umi in reverse complement reads

## Input

* `--in-fq1`: path to the forward FASTQ (`path/to/R1.fastq.gz`)
* `--in-fq2`: path to the reverse FASTQ (`path/to/R2.fastq.gz`)

## Output

* `--out-fq1`: path to the output deduplicated forward FASTQ (`path/to/cR1.fastq.gz`)
* `--out-fq2`: path to the output deduplicated reverse FASTQ (`path/to/cR2.fastq.gz`)
* `--out-json`: path to the output json with umi metrics and total reads (`path/to/calib.json`)

## How to run

```bash
docker build -t pyumi .

docker run \
   -v ./unit_tests/test_data:/root/ \
   pyumi \
   --in-fq1 /root/age_ig_s1_R1_umi.fastq.gz \
   --in-fq2 /root/age_ig_s1_R2_rc_umi.fastq.gz \
   --out-fq1 /root/bR1.fastq.gz \
   --out-fq2 /root/bR2.fastq.gz \
   --out-json /root/calib.json \
   --fq1-pattern "^UMI:N{13}" # UMI - the first 13 nucleotides
```