process GetLinks {
    publishDir "${params.outdir}/calib", mode: 'copy', overwrite: false
    container 'downloader'

    input:
        val sample_id

    output:
        env link_1, emit: link_1
        env link_2, emit: link_2

    script:
        """
        link_1=\$(ffq --ftp $sample_id | jq -r '.[0] | .url')
        link_2=\$(ffq --ftp $sample_id | jq -r '.[1] | .url')
        """
}

process Download {
    publishDir "${params.outdir}/calib", mode: 'copy', overwrite: false
    container 'downloader'

    input:
        val sample_id
        val link
        val reads_to_save
        val read

    output:
        path "${sample_id}_${read}.fastq.gz", emit: fq

    script:
        """
        if [[ "${reads_to_save}" == "all" ]]; then
            wget -q "${link}" -O "${sample_id}_${read}.fastq.gz"
        else
            lines_to_save=\$(( ${reads_to_save} * 4 ))
            wget -qO- "${link}" | zcat | head -n "\${lines_to_save}" | gzip > "${sample_id}_${read}.fastq.gz"
        fi
        """
}

process Downsample {
    publishDir "${params.outdir}/calib", mode: 'copy', overwrite: false
    container 'downloader'

    input:
        val fastq
        val reads_to_save
        val read

    output:
        path "R${read}.fastq.gz", emit: fq

    script:
        """
        lines_to_save=\$( expr $reads_to_save * 4 )
        zcat $fastq | head -\$lines_to_save | gzip > R${read}.fastq.gz
        """
}