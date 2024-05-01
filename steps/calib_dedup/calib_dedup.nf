process CalibDedup {
//     publishDir "${params.outdir}/calib", mode: 'copy', overwrite: false
    container 'calib-dedup'

    input:
        path fq1
        path fq2
    output:
        path "cR1.fastq.gz", emit: fq1
        path "cR2.fastq.gz", emit: fq2
        path "calib.json", emit: json
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq1 cR1.fastq.gz \
            --out-fq2 cR2.fastq.gz \
            --out-json calib.json \
            --kmer-size 7 \
            --minimizer-count 6 \
            --minimizer-threshold 6 \
            --error-tolerance 2 \
            --fq1-barcode-pattern "^UMI:N{12}"
        """
}