# Vidjil step

Detects V(D)J segments in fastq files via **[Vidjil](https://www.vidjil.org/doc/)** tool and saves them into fasta file.

## Parameters

* `--mode`: there are three modes:
  * "detect" - returns .fa.gz file with detected V(D)J segments, 
  * "annotate" - returns .tsv.gz file with V(D)J annotation (in AIRR format, but with less column number that in IgBLAST)
  * "all" - return .fa.gz and .tsv.gz files
* `--receptor`: receptor, that V(D)J we want to detect or annotate. Should be "TCR", "BCR" or "all" ("TCR" and "BCR").
* `--organism`: organism, that V(D)J we want to detect or annotate. Should be "human", "rat" or "mouse".

## Input

* `--in-fq1`: path to the forward fastq (`path/to/mR1.fq.gz`)
* `--in-fq2`: path to the reverse fastq (`path/to/mR2.fq.gz`)
* `--in-fq12`: path to the merged fastq (`path/to/mR12.fq.gz`)
* `--in-ref`: reference vidjil germline (`path/to/vidjil.germline.tar.gz`)

## Output

* `--out-annotation`: path to the annotation (`path/to/annotation.tsv.gz`)
* `--out-fasta`: path to the fasta with found V(D)J segments (`path/to/vidjil.fa.gz`)
* `--logs`: path to the vidjil logs (`path/to/vidjil.log`)

## Build archive with V(D)J reference in vidjil format:

Run script:
```bash
bash build_ref.sh
``` 

## How to run

```bash
docker build -t vidjil .

# should contain: mR1.fastq.gz, mR2.fastq.gz, mR12.fastq.gz and vidjil.germline.tar.gz
FOLDER_WITH_DATA=path/to/your/folder

docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   vidjil \
   --in-fq1 /root/mR1.fastq.gz \
   --in-fq2 /root/mR2.fastq.gz \
   --in-fq12 /root/mR12.fastq.gz \
   --mode 'all' \
   --receptor 'TCR' \
   --organism 'human' \
   --ref /root/vidjil.germline.tar.gz \
   --out-fasta /root/vidjil.TCR.fasta.gz \
   --out-annotation /root/vidjil.TCR.tsv.gz \
   --logs /root/vidjil.TCR.log

docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   vidjil \
   --in-fq1 /root/mR1.fastq.gz \
   --in-fq2 /root/mR2.fastq.gz \
   --mode 'all' \
   --receptor 'BCR' \
   --organism 'human' \
   --ref /root/vidjil.germline.tar.gz \
   --out-fasta /root/vidjil.BCR.fasta.gz \
   --out-annotation /root/vidjil.BCR.tsv.gz \
   --logs /root/vidjil.BCR.log
```