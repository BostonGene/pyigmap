process GetLinks {
    // labels are defined in conf/base.config
    // label "Get_links_from_SRA_ID"

    input:
        val sample_id

    output:
        env link_1, emit: link_1
        env link_2, emit: link_2

    script:
        """
        urls=\$(ffq --ftp $sample_id | jq -r '.[].url')

        link_1=\$(echo "\$urls" | sed -n '1p')
        if [ -z "\$link_1" ]; then
            echo "Error: failed to retrieve FASTQ link for $sample_id" >&2
            exit 1
        fi

        link_2=\$(echo "\$urls" | sed -n '2p' || true)
        """
}

process Download {
    // labels are defined in conf/base.config
    // label "Download_FASTQ_by_link"

    input:
        val sample_id
        val link
        val read

    output:
        path "${sample_id}_R${read}.fastq.gz", emit: fq

    script:
        """
        if [[ "${params.first_reads}" == "all" ]]; then
            wget "${link}" -O "${sample_id}_R${read}.fastq.gz"
        else
            let lines_to_save=${params.first_reads}*4
            wget -O - "${link}" | zcat | head -n "\${lines_to_save}" | gzip > "${sample_id}_R${read}.fastq.gz" || true
        fi
        """
}

process Downsample {
    // labels are defined in conf/base.config
    label "process_single"

    input:
        path fastq
        val read

    output:
        path "R${read}.fastq.gz", emit: fq

    script:
        """
        #!/bin/bash
        let lines_to_save=${params.first_reads}*4
        zcat $fastq | head -\${lines_to_save} | gzip > R${read}.fastq.gz
        """
}