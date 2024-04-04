# IgBLAST step

This step is a wrapping of [IgBlast](https://ncbi.github.io/igblast/) V(D)J mapping tool.

## Parameters

* `--receptor`: receptor name: "TCR" or "BCR"
* `--organism`: organism name: "human" or "mouse"

## Input

* `--in_fq1` (**optional**): path to the first fastq (`path/to/mR1.fq.gz`)
* `--in-fq2` (**optional**): path to the second fastq (`path/to/mR2.fq.gz`)
* `--in-fq12` (**optional**): path to the merged fastq (`path/to/mR12.fq.gz`)
* `--in-fasta` (**optional**): path to the preprocessed fasta with detected V(D)J segments (`path/to/vidjil.fq.gz`)
* `--in-ref`: path to the archive with reference V(D)J segments (`path/to/iblast.reference.tar.gz`)

## Output

* `out_annotation`: path to the output IgBLAST annotation (`path/to/igblast_annotation.tsv.gz`)

## Build archive with V(D)J reference in IgBLAST format:

Run script:
```bash
bash build_ref.sh

cp /tmp/igblast.reference.tar.gz .
``` 

## How run

```bash
docker build -t igblast .

# should contain: vidjil.TCR.fasta.gz, vidjil.BCR.fasta.gz and igblast.reference.tar.gz
FOLDER_WITH_DATA=path/to/your/folder

docker run \
    -v ${FOLDER_WITH_DATA}:/root/ \
    igblast \
    --in-fasta /root/vidjil.TCR.fasta.gz \
    --in-ref /root/igblast.reference.tar.gz \
    --receptor 'TCR' \
    --organism human \
    --out-annotation /root/raw_annotation.TCR.tsv.gz

docker run \
    -v ${FOLDER_WITH_DATA}:/root/ \
    igblast \
    --in-fasta /root/vidjil.BCR.fasta.gz \
    --in-ref /root/igblast.reference.tar.gz \
    --receptor 'BCR' \
    --organism human \
    --out-annotation /root/raw_annotation.BCR.tsv.gz
```