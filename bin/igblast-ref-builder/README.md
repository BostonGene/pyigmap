## How to build V(D)J reference in IgBLAST format

Build docker image with tool to generate IgBLAST reference:
```bash
docker build -t igblast-ref-builder .
```

### Only major alleles

If you need to keep only *01 (major) alleles, execute:
```bash
docker run --rm \
    -v ./:/tmp \
    igblast-ref-builder \
    --out-archive /tmp/igblast.reference.major_allele.tar.gz
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
    igblast-ref-builder \
    --all-alleles --out-archive /tmp/igblast.reference.all_alleles.tar.gz
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