# IgBLAST component

This component is a wrapping of [IgBLAST](https://ncbi.github.io/igblast/) V(D)J mapping tool.

## Input

* `--in-fastq` (**optional**): Path to the FASTQ file.
* `--in-fasta` (**optional**): Path to the FASTA file with detected V(D)J segments.
* `--ref`: Path to the archive with reference V(D)J segments.
* `--receptor`: Receptor type: "BCR", "TCR" or "all". Default: `"all"`.
* `--organism`: organism name: "human" or "mouse"
* `--reads-chunk-size`: Count of sequences processed in one run of IgBLAST tool. Default: `50_000`.

## Output

* `--out-annotation`: Path to the output IgBLAST annotation. Default: `igblast_annotation.tsv.gz`

## How to build V(D)J reference in IgBLAST format

Build docker image with tool to generate IgBLAST reference:
```bash
docker build --target build-ref -t build-ref .
```

### Only major alleles

If you need to keep only *01 (major) alleles, execute:
```bash
docker run --rm \
    -v ./:/tmp \
    build-ref \
    --out-archive igblast.reference.major_allele.tar.gz
```

Reference will contain sequences with only major allele (*01):
```
>IGHD2-2*01
aggatattgtagtagtaccagctgctatgcc
```

### All alleles

Or you can keep all alleles (*01, *02, etc.) by specifying `--all-alleles` flag:
```bash
docker run --rm \
    -v ./:/tmp \
    build-ref \
    --all-alleles --out-archive igblast.reference.all_alleles.tar.gz
```

Reference will contain sequences with all alleles (*01, *02, *03, etc.):
```
>IGHD2-2*01
aggatattgtagtagtaccagctgctatgcc
>IGHD2-2*02
aggatattgtagtagtaccagctgctatacc
>IGHD2-2*03
tggatattgtagtagtaccagctgctatgcc
```