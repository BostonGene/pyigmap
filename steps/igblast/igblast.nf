process IgBlast {
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