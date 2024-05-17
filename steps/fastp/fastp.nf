process Fastp {
//     publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
        path fq1
        path fq2
    output:
        path "mR1.fastq.gz", emit: fq1
        path "mR2.fastq.gz", emit: fq2
        path "mR12.fastq.gz", emit: fq12
        path "fastp.json", emit: json
        path "fastp.html", emit: html
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq1 mR1.fastq.gz \
            --out-fq2 mR2.fastq.gz \
            --out-fq12 mR12.fastq.gz \
            --disable "length_filtering" "adapter_trimming" "quality_filtering" \
            --merge \
            --out-html fastp.html \
            --out-json fastp.json
        """
}