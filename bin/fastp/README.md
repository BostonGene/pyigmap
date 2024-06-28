# Fastp step

This step is a wrapping of [fastp](https://github.com/OpenGene/fastp) tool designed to provide fast all-in-one preprocessing for FASTQ files.

## Parameters:
* `--merge-reads`: merges forward and reverse FASTQ.
* `--disable-filters`: disables modes. Example:
    ```bash
    --disable-filters "length_filtering" # disables filtering reads shorter than 15bp
    --disable-filters "quality_filtering" # disables quality filtering (if >40% bases have quality <20)
    --disable-filters "adapter_trimming" # disables trimming Illumina adapters
    ```
* `--mock-merge-reads`: enables mock merging of not overlapped forward and reverse reads with a selected insert size (distance)
    * `--inner-distance-size`: inner distance between reads (used in mock merging) 

Example of mering overlapped reads:

```
|------------------------| Full sequence
|----------ATG>            R1
          <ATG-----------| R2
|----------ATG-----------| Concatenated (overlap)
```

Example of mock mering not-overlapped reads with `--insert-distance-size 3`:

```
|------------------------| Full sequence
|--------->                R1
              <----------| R2
|----------NNN-----------| Concatenated (no overlap)
```

## Input

  * `--in-fq1`: path to the forward fastq (`path/to/R1.fastq.gz`)
  * `--in-fq2` (**optional**): path to the reverse fastq (`path/to/R2.fastq.gz`)

## Output

  * `--out-fq1` (**optional**): output forward fastq
  * `--out-fq2` (**optional**): output reverse fastq
  * `--out-fq12` (**optional**): output merged fastq
  * `--json`: output json with metrics
  * `--html`: output html with metrics

## Testing

To run unit tests, execute:

```bash
pip install pytest
docker build -t fastp .
pytest unit_tests -vv
```