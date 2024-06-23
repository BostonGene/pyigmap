process CalibDedup {
    // labels are defined in conf/base.config
    label "process_low"

    input:
        path fq1
        path fq2
        path pyumi_json
    output:
        path params.out_calib_dedup_fq1, emit: fq1
        path params.out_calib_dedup_fq2, emit: fq2
    script:
        """
        fq1_umi_len=\$(cat $pyumi_json | jq -r '.fq1_umi_length')
        fq2_umi_len=\$(cat $pyumi_json | jq -r '.fq2_umi_length')

        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --fq1-umi-length \$fq1_umi_len \
            --fq2-umi-length \$fq2_umi_len \
            --out-fq1 ${params.out_calib_dedup_fq1} \
            --out-fq2 ${params.out_calib_dedup_fq2} \
            --kmer-size ${params.kmer_size} \
            --minimizer-count ${params.minimizer_count} \
            --minimizer-threshold ${params.minimizer_threshold} \
            --error-tolerance ${params.error_tolerance}
        """
}