process Vidjil {
    input:
        path fq1
        path fq2
        path fq12
        path ref
    output:
        path "vidjil.fasta.gz", emit: fasta
        path "vidjil.log", emit: logs
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --in-fq12 $fq12 \
            --ref $ref \
            --out-fasta vidjil.fasta.gz \
            --logs vidjil.log
        """
}