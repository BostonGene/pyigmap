process CDR3ErrorCorrector {
    input:
        path raw_annotation
        path olga_models
        path json
    output:
        path "pyigmap.tar.gz", emit: archive
    script:
        """
        python3.9 /usr/local/run.py \
            --in-annotation $raw_annotation \
            --filter-pgen-all 0 \
            --only-functional \
            --only-canonical \
            --remove-chimeras \
            --clonotype-collapse-factor 0.05 \
            --olga-models $olga_models \
            --out-corrected-annotation corrected_annotation.tsv \
            --in-json $json \
            --out-json stat.json \
            --out-archive pyigmap.tar.gz
        """
}