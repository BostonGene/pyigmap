process IgBlastFASTA {
    // labels are defined in conf/base.config
    label "process_medium"

    input:
        path fasta
        path ref
    output:
        path params.out_igblast_annotation, emit: annotation
    script:
        """
        python3.9 /usr/local/src/run.py \
            --in-fasta \${PWD}/$fasta \
            --ref \${PWD}/$ref \
            --receptor ${params.igblast_receptor} \
            --organism ${params.igblast_organism} \
            --out-annotation \${PWD}/${params.out_igblast_annotation}
        """
}

process IgBlastFASTQ {
    // labels are defined in conf/base.config
    label "process_medium"

    input:
        path fq1
        path fq2
        path fq12
        path ref
    output:
        path params.out_igblast_annotation, emit: annotation
    script:
        """
        python3.9 /usr/local/src/run.py \
            --in-fastq \${PWD}/$fq1 \${PWD}/$fq2 \${PWD}/$fq12 \
            --ref \${PWD}/$ref \
            --receptor ${params.igblast_receptor} \
            --organism ${params.igblast_organism} \
            --out-annotation \${PWD}/${params.out_igblast_annotation}
        """
}

process IgBlastMockFASTQ {
    // labels are defined in conf/base.config
    label "process_medium"

    input:
        path fq12
        path ref
    output:
        path params.out_igblast_annotation, emit: annotation
    script:
        """
        python3.9 /usr/local/src/run.py \
            --in-fastq \${PWD}/$fq12 \
            --in-ref \${PWD}/$ref \
            --receptor ${params.igblast_receptor} \
            --organism ${params.igblast_organism} \
            --out-annotation \${PWD}/${params.out_igblast_annotation}
        """
}