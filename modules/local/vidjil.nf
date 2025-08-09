process Vidjil {
    // labels are defined in conf/base.config
    // label "process_medium"

    input:
        path fq12
        path ref
    output:
        path params.out_vidjil_fasta, emit: fasta
        path params.out_vidjil_logs, emit: logs
    script:
        """
        python3 /usr/local/src/run.py \
            --in-fastq $fq12 \
            --vdj-ref $ref \
            --out-fasta ${params.out_vidjil_fasta} \
            --debug \
            --logs ${params.out_vidjil_logs}
        """
}