# cdr3nt-error-corrector step

This step removes spurious rearrangements via [OLGA](https://github.com/statbiophys/OLGA) tool and/or by means of comparing igblast alignment scores. The script here operates on clonotype abundance tables by filtering out bad records and assigning erroneous records ('artificial diversity' produced by PCR and sequencing errors) to their parent records based on a certain parent-to-child ratio.

## Definition of erroneous and spurious records

In case only Vidjil tool was used to extract clonotypes, they are defined based on V-J mapping and junction sequence (here it is an alias for CDR3 sequence with bounding Cys and Phe/Trp codons included). Junction sequences not containing stop codons ``*`` and frameshifts ``_`` are defined as ``functional`` as otherwise they give rise to mRNA molecules that could not be ``productive`` and translated to proteins.
Additionally, junction sequences are required to start with a Cys codon and end with a Phe/Trp codon (actually ``[FW]G.G`` [motif](https://www.pnas.org/doi/10.1073/pnas.121101598)) to form a proper structure called [omega-loop](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753249/). These junction sequences are called ``canonical``.

In case IgBlast is used following Vidjil, more detailed V/J mappings are provided, including ``stop_codon``, ``vj_in_frame``, ``v_frameshift`` and ``productive``. The ``stop_codon``, ``vj_in_frame`` fields are the same as ``*`` and ``_`` in Vidjil junction, i.e. ``functional = !stop_codon || vj_in_frame``. The ``v_frameshift`` and ``productive`` are defined in terms of entire sequence, including checking that V/J genes are not pseudogenes and give rise to a productive TCR or antibody protein (e.g. human TRBV21-1 is a pseudogene).

By default, all the rules for filtering out spurious sequences are defined in the following way:

* among all chains where a given well-defined ``junction`` is present, only the ``chain`` having lowest ``v_support`` (checked first) or ``j_support`` is selected to prevent miss-mappings (N.B. lowest is due to the fact that support is a proxy for match P-value in igblast)
* chimeras defined as mapping of V and J to different loci are dropped, except for TRA and TRD chains (due to presence of TRAV-X/DV-Y genes)
* V(D)J junction is ``functional``:
  * well-defined (non-empty entry for ``junction``)
  * does not contain ``*`` or ``_``
* junction is canonical, starting with ``C`` and ending with either ``F`` or ``W``
* either supported by several reads (``duplicate_count > 1``), or having ``Pgen > 0``
* the entire sequence is productive: V and J genes are not pseudogenes (currently not checked), ``junction`` is ``functional`` for Vidjil and ``productive`` for IgBlast
* additional filtering involves ``pgen_threshold`` calculated using OLGA software that estimates the theoretical likelihood of a given V(D)J rearrangements, can be used for filtering spurious rearrangements that are likely to be replicate artefacts


## Parameters

* `--filter-pgen-all <pgen_threshold>` (**optional**): calculates generation probability of junctions, all clonotypes with `pgen <= pgen_threshold` will be removed
* `--filter-pgen-singletons <pgen_threshold>` (**optional**): calculate generation probability of junctions, all clonotypes with `duplicate_count == 1 && pgen <= pgen_threshold` will be removed
* `--skip-pgen-calculation` (**optional**): skip calculation of generation probability of junctions
* `--clonotype-collapse-factor` (**optional**): value, that specifies the parent-to-child ratio in order to define erroneous (based on Levenstein distance of 1 and read count ratio) records that should be clustered to their parent clonotypes. Defaults to ``0.05``
* `--remove-chimeras` (**optional**): filter out chimeras, that have different locus in V-/J-segments (except for TRA and TRD)
* `--only-functional` (**optional**): filter out non-functional clonotypes
* `--only-canonical` (**optional**): filter out non-canonical clonotypes
* `--only-productive` (**optional**): filter out non-productive clonotypes (if IgBlast fields are present, for Vidjil same as ``--only-functional``)

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
           "no_j_call": 100, // reads with j_call = None
           "no_junction": 1000 // reads with junction = None
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
   --filter-pgen-all 0 \
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
