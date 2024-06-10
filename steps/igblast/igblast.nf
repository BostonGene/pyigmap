process IgBlastFASTA {
    input:
        path fasta
        path ref
    output:
        path params.out_igblast_annotation, emit: annotation
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fasta \${PWD}/$fasta \
            --in-ref \${PWD}/$ref \
            --receptor ${params.igblast_receptor} \
            --organism ${params.igblast_organism} \
            --out-annotation \${PWD}/${params.out_igblast_annotation}
        """
}

process IgBlastFASTQ {
    input:
        path fq1
        path fq2
        path fq12
        path ref
    output:
        path params.out_igblast_annotation, emit: annotation
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 \${PWD}/$fq1 \
            --in-fq2 \${PWD}/$fq2 \
            --in-fq12 \${PWD}/$fq12 \
            --in-ref \${PWD}/$ref \
            --receptor ${params.igblast_receptor} \
            --organism ${params.igblast_organism} \
            --out-annotation \${PWD}/${params.out_igblast_annotation}
        """
}