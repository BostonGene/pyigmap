process CalibDedup {
    input:
        path fq1
        path fq2
    output:
        path params.out_calib_dedup_fq1, emit: fq1
        path params.out_calib_dedup_fq2, emit: fq2
        path params.out_calib_dedup_json, emit: json
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq1 ${params.out_calib_dedup_fq1} \
            --out-fq2 ${params.out_calib_dedup_fq2} \
            --out-json ${params.out_calib_dedup_json} \
            --kmer-size ${params.kmer_size} \
            --minimizer-count ${params.minimizer_count} \
            --minimizer-threshold ${params.minimizer_threshold} \
            --error-tolerance ${params.error_tolerance} \
            --fq1-barcode-pattern ${params.fq1_barcode_pattern}
        """
}