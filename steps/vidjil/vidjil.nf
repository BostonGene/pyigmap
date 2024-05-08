process Vidjil {
//     publishDir "${params.outdir}/vidjil", mode: 'copy'
    container 'vidjil'

    input:
        path fq1
        path fq2
        path fq12
        path ref
        val receptor
    output:
        path "vidjil.${receptor}.fasta.gz", emit: fasta
        path "vidjil.${receptor}.tsv.gz", emit: tsv
        path "vidjil.${receptor}.log", emit: logs
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --in-fq12 $fq12 \
            --mode 'all' \
            --receptor $receptor \
            --organism 'human' \
            --ref $ref \
            --out-fasta vidjil.${receptor}.fasta.gz \
            --out-annotation vidjil.${receptor}.tsv.gz \
            --logs vidjil.${receptor}.log
        """
}