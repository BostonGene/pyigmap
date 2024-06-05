process IgBlast {
    input:
        path fasta
        path ref
    output:
        path "raw_annotation.tsv.gz", emit: annotation
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fasta \${PWD}/$fasta \
            --in-ref \${PWD}/$ref \
            --receptor "all" \
            --organism human \
            --out-annotation \${PWD}/raw_annotation.tsv.gz
        """
}