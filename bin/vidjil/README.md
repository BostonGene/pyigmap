# Vidjil component

Detects reads with CDR3 in FASTQ files via **[Vidjil](https://www.vidjil.org/doc/)** tool and saves them into FASTA file.

## Input

* `--in-fastq`: Path to the FASTQ.
* `--vdj-ref`: Path to the vidjil reference.
* `--debug`: Enable debug mode to save vidjil logs. Default: `false`.
* `--organism`: organism, that V(D)J we want to detect or annotate. Should be "human", "rat" or "mouse".

## Output

* `--out-fasta`: Path to the fasta with found V(D)J segments. Default: `vidjil.fasta.gz`.
* `--logs`: Path to the vidjil logs. Default: `vidjil.log`.
* `--reads-chunk-size`: The maximum number of reads in the chunk to process using Vidjil tool. Default: `5000000`.

## How to build vidjil reference

Run script:
```bash
bash build_ref.sh
```