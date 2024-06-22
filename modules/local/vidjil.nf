process Vidjil {
    // labels are defined in conf/base.config
    label "process_medium"

    input:
        path fq12
        path ref
    output:
        path params.out_vidjil_fasta, emit: fasta
        path params.out_vidjil_logs, emit: logs
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq12 $fq12 \
            --ref $ref \
            --out-fasta ${params.out_vidjil_fasta} \
            --logs ${params.out_vidjil_logs}
        """
}