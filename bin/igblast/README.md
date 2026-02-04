# IgBLAST component

This component is a wrapping of [IgBLAST](https://ncbi.github.io/igblast/) V(D)J mapping tool.

## Input

* `--in-fastq` (**optional**): Path to the FASTQ file.
* `--in-fasta` (**optional**): Path to the FASTA file with detected V(D)J segments.
* `--ref`: Path to the archive with reference V(D)J segments.
* `--receptor`: Receptor type: "BCR", "TCR" or "all". Default: `"all"`.
* `--reads-chunk-size`: Count of sequences processed in one run of IgBLAST tool. Default: `50_000`.

## Output

* `--out-annotation`: Path to the output IgBLAST annotation. Default: `igblast_annotation.tsv.gz`