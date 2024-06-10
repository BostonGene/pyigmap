process CDR3ErrorCorrector {
    input:
        path raw_annotation
        path olga_models
        path json
    output:
        path params.out_corrected_annotation, emit: archive
    script:
        """
        python3.9 /usr/local/run.py \
            --in-annotation $raw_annotation \
            --filter-pgen-all ${params.pgen_threshold} \
            ${params.enabled_filters} \
            --clonotype-collapse-factor ${params.clonotype_collapse_factor} \
            --olga-models $olga_models \
            --out-corrected-annotation ${params.out_corrected_annotation} \
            --in-json $json \
            --out-json ${params.out_stat_json} \
            --out-archive ${params.out_archive}
        """
}