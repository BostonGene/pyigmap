# pyigmap


## Quick start

1. This pipleine requires Docker, Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html). Be sure, that you 
2. Run this command to install nextflow and build container steps:
```bash
sudo apt install make
make
```
3. Start running your own analysis!

```
./nextflow run main.nf --mode amplicon --outdir results --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz
```