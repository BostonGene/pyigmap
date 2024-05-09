# cdr3nt-error-corrector step

This step removes spurious rearrangements via [OLGA](https://github.com/statbiophys/OLGA) tool, corrects and filters out bad clones.

## Parameters
* `--pgen-threshold` (**optional**): probability generation (pgen) threshold; all clones with `pgen <= pgen_threshold` will be removed. If you need disable this filtration, remove **pgen_threshold** from your `values.yml`.
* `--calculate-pgen` (**optional**): calculate generation probability of clonotypes. Automatically on if `--pgen-threshold` is set and not 0
* `--only-functional` (**optional**): filter out non-functional clonotypes
* `--only-productive` (**optional**): filter out non-productive clonotypes
* `--filter-pgen-singletons` (**optional**): Filter out singleton clones with duplicate_count=1 and pgen<=pgen_threshold
* `--clonotype-collapse-factor` (**optional**): value, that involved in clonotypes collapsing
* `--remove-chimeras` (**optional**): filter out chimeras, that have different locus in V-/J-segments

## Input

* `--in-tcr-annotation`: path to the raw TCR annotation, which we need to correct (`path/to/raw_annotation.TCR.tsv.gz`)
* `--in-bcr-annotation`: path to the raw BCR annotation, which we need to correct (`path/to/raw_annotation.BCR.tsv.gz`)
* `--olga-models`: path to the archive with [OLGA](https://github.com/statbiophys/OLGA/tree/master/olga/default_models) models (`path/to/olga-models.tar.gz`)
* `--in-json`: path to the json files with total reads count (`path/to/stat.json`)

## Output

* `--out-archive`: path to the output archive with corrected annotation and json with total reads count (`path/to/pyigmap.tar.gz`)
  * `pyigmap.tar.gz`:
    * `corrected_annotation.tsv`
    * `stat.json`:
      ```json5
        {
           "total_reads": 100000, // total reads count in the input fastq
           "IGH_aligned_reads": 100, // reads that aligned to IGH locus
           ..., // reads that aligned to ... locus
           "TRB_aligned_reads": 20, // reads that aligned to TRB locus
           "no_v_call": 10, // reads with v_call = None
           "no_j_call": 100 // reads with j_call = None
        }
      ```


## Build archive with OLGA models:

Run script:
```bash
bash build_ref.sh
``` 

## How to run

1. Build an image
```bash
docker build --target tool -t cdr3nt-error-corrector .
```

2. Be sure, that you have all input files

```bash
$ FOLDER_WITH_DATA=path/to/your/folder
$ ls $FOLDER_WITH_DATA
raw_annotation.TCR.tsv.gz
raw_annotation.BCR.tsv.gz
```

3. Run container
```bash
docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   -v ./olga-models.tar.gz:/root/olga-models.tar.gz \
   cdr3nt-error-corrector \
   --in-tcr-annotation /root/raw_annotation.TCR.tsv.gz \
   --in-bcr-annotation /root/raw_annotation.BCR.tsv.gz \
   --pgen-threshold 0 \
   --only-functional \
   --remove-chimeras \
   --clonotype-collapse-factor 0.05 \
   --olga-models /root/olga-models.tar.gz \
   --out-corrected-annotation /root/corrected_annotation.tsv \
   --in-json /root/test_fastp.json \
   --out-json /root/stat.json \
   --out-archive /root/pyigmap.tar.gz # archive with final results will be saved into ./unit_tests/test_data/
```

4. Outputs will be here
```bash
$ ls $FOLDER_WITH_DATA
raw_annotation.TCR.tsv.gz
raw_annotation.BCR.tsv.gz
corrected_annotation.tsv
stat.json
pyigmap.tar.gz # <-- OUTPUT
```