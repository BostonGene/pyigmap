process CDR3ErrorCorrector {
    publishDir "${params.outdir}", mode: 'copy'
    container 'cdr3nt-error-corrector'

    input:
        path bcr_annotation
        path tcr_annotation
        path olga_models
        path json
    output:
        path "pyigmap.tar.gz", emit: archive
    script:
        """
        python3.9 /usr/local/run.py \
            --in-tcr-annotation $bcr_annotation \
            --in-bcr-annotation $tcr_annotation \
            --pgen-threshold 0 \
            --only-functional \
            --remove-chimeras \
            --clonotype-collapse-factor 0.05 \
            --olga-models $olga_models \
            --out-corrected-annotation corrected_annotation.tsv \
            --in-json $json \
            --out-json stat.json \
            --out-archive pyigmap.tar.gz
        """
}