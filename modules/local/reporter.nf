process Reporter {
    // labels are defined in conf/base.config
    // label "process_low"

    input:
        path fq1_pyumi
        path fq2_pyumi
        path fq1_calib
        path fq2_calib
    output:
        path params.out_report_file, emit: html
    script:
        """
        python3.10 /usr/local/src/run.py \
            --in-fq1-pyumi $fq1_pyumi \
            --in-fq2-pyumi $fq2_pyumi \
            --in-fq1-calib $fq1_calib \
            --in-fq2-calib $fq2_calib \
            --report-file ${params.out_report_file}
        """
}