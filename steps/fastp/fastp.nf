process Fastp {
//     publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
        path fq1
        path fq2
    output:
        path "mR12.fastq.gz", emit: fq12
        path "fastp.json", emit: json
        path "fastp.html", emit: html
    script:
        """
        python3.9 /usr/local/run.py \
            --in-fq1 $fq1 \
            --in-fq2 $fq2 \
            --out-fq12 mR12.fastq.gz \
            --disable "length_filtering" "adapter_trimming" "quality_filtering" \
            --mock-merge \
            --insert-size ${params.insert_size} \
            --out-html fastp.html \
            --out-json fastp.json
        """
}