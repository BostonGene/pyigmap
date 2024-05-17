process IgBlast {
//     publishDir "${params.outdir}/igblast", mode: 'copy'

    input:
        path fasta
        path ref
        val receptor
    output:
        path "raw_annotation.${receptor}.tsv.gz", emit: annotation
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fasta \${PWD}/$fasta \
            --in-ref \${PWD}/$ref \
            --receptor ${receptor} \
            --organism human \
            --out-annotation \${PWD}/raw_annotation.${receptor}.tsv.gz
        """
}