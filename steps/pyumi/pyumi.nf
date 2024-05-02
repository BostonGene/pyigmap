process PyUMI {
//     publishDir "${params.outdir}/pyumi", mode: 'copy'
    container 'pyumi'

    input:
        path fq1
        path fq2
    output:
        path "bR1.fastq.gz", emit: fq1
        path "bR2.fastq.gz", emit: fq2
        path "pyumi.json", emit: json
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq1 cR1.fastq.gz \
            --out-fq2 cR2.fastq.gz \
            --out-json pyumi.json \
            --fq1-pattern "^UMI:N{12}"
        """
}