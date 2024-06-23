process CDR3ErrorCorrector {
    // labels are defined in conf/base.config
    label "process_low"

    input:
        path raw_annotation
        path olga_models
        path json
    output:
        path params.out_archive, emit: archive
    script:
        def librarySpecifiedOptions = params.library == "rnaseq" ? params.rnaseq_corrector_options : params.amplicon_corrector_options
        """
        python3.9 /usr/local/run.py \
            --in-annotation $raw_annotation \
            ${params.default_corrector_options} \
            $librarySpecifiedOptions \
            --olga-models $olga_models \
            --out-corrected-annotation ${params.out_corrected_annotation} \
            --in-json $json \
            --out-json ${params.out_stat_json} \
            --out-archive ${params.out_archive}
        """
}