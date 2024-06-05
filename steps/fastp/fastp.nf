process Fastp {
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
            --disable "length_filtering" "adapter_trimming" "quality_filtering" \
            --merge \
            --out-html ${params.out_fastp_html} \
            --out-json ${params.out_fastp_json}
        """
}