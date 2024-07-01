process FastpMerge {
    // labels are defined in conf/base.config
    label "process_low"

    input:
        path fq1
        path fq2
    output:
        path params.out_fastp_fq1, emit: fq1
        path params.out_fastp_fq2, emit: fq2
        path params.out_fastp_fq12, emit: fq12
        path params.out_fastp_json, emit: json
        path params.out_fastp_html, emit: html
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq1 ${params.out_fastp_fq1} \
            --out-fq2 ${params.out_fastp_fq2} \
            --out-fq12 ${params.out_fastp_fq12} \
            --disable-filters ${params.disable} \
            --merge-reads \
            --html ${params.out_fastp_html} \
            --json ${params.out_fastp_json}
        """
}

process FastpMockMerge {
    // labels are defined in conf/base.config
    label "process_low"

    input:
        path fq1
        path fq2
    output:
        path params.out_fastp_fq12, emit: fq12
        path params.out_fastp_json, emit: json
        path params.out_fastp_html, emit: html
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq12 ${params.out_fastp_fq12} \
            --disable-filters ${params.disable} \
            --mock-merge-reads \
            --inner-distance-size ${params.insert_size} \
            --reads-chunk-size 5000000 \
            --html ${params.out_fastp_html} \
            --json ${params.out_fastp_json}
        """
}