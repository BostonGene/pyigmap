process PyUMI {
    input:
        path fq1
        path fq2
    output:
        path params.out_pyumi_fq1, emit: fq1
        path params.out_pyumi_fq2, emit: fq2
        path params.out_pyumi_json, emit: json
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --fq1-pattern "${params.fq1_pattern ? params.fq1_pattern : ''}" \
            --fq2-pattern "${params.fq2_pattern ? params.fq2_pattern : ''}" \
            --out-fq1 ${params.out_pyumi_fq1} \
            --out-fq2 ${params.out_pyumi_fq2} \
            --out-json ${params.out_pyumi_json}
        """
}