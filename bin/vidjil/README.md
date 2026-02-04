# Vidjil component

Detects reads with CDR3 in FASTQ files via **[Vidjil](https://www.vidjil.org/doc/)** tool and saves them into FASTA file.

## Input

* `--in-fastq`: Path to the FASTQ.
* `--vdj-ref`: Path to the vidjil reference.
* `--debug`: Enable debug mode to save vidjil logs. Default: `false`.

## Output

* `--out-fasta`: Path to the fasta with found V(D)J segments. Default: `vidjil.fasta.gz`.
* `--out-json`: Path to the JSON with QC metrics. Default: `vidjil.json`.