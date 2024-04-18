# CDR3nt-error-corrector step

This step removes spurious rearrangements via [OLGA](https://github.com/statbiophys/OLGA) tool, corrects and filters out bad clones.

## Parameters
* `--pgen-threshold` (**optional**): probability generation (pgen) threshold; all clones with `pgen <= pgen_threshold` will be removed. If you need disable this filtration, remove **pgen_threshold** from your `values.yml`.
* `--calculate-pgen` (**optional**): calculate generation probability of clonotypes. Automatically on if `--pgen-threshold` is set and not 0
* `--only-functional` (**optional**): filter out non-functional clonotypes
* `--only-productive` (**optional**): filter out non-productive clonotypes
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
    * `stat.json` (structure: `{total_reads: 100000}`, where `100000` -- total reads count parsed from `--in-json`)


## Build archive with OLGA models:

Run script:
```bash
bash build_ref.sh
``` 

## How to run

```bash
docker build -t cdr3nt-error-corrector .

# should contain: raw_annotation.TCR.tsv.gz and raw_annotation.BCR.tsv.gz, olga-models.tar.gz and calib.json (or fastp.json)
FOLDER_WITH_DATA=path/to/your/folder

docker run \
   -v ${FOLDER_WITH_DATA}:/root/ \
   cdr3nt-error-corrector \
   --in-tcr-annotation /root/raw_annotation.TCR.tsv.gz \
   --in-bcr-annotation /root/raw_annotation.BCR.tsv.gz \
   --pgen-threshold 0 \
   --only-functional \
   --remove-chimeras \
   --clonotype-collapse-factor 0.05 \
   --olga-models /root/olga-models.tar.gz \
   --out-corrected-annotation /root/corrected_annotation.tsv \
   --in-json /root/calib.json \
   --out-json /root/stat.json \
   --out-archive /root/pyigmap.tar.gz # archive with final results
```

## Run tests

```bash
python3 -m venv venv
. env/bin/activate
pip3 install -r requirements.txt

pytest unit_tests/ -v
```