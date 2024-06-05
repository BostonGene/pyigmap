process GetLinks {

    input:
        val sample_id

    output:
        env link_1, emit: link_1
        env link_2, emit: link_2

    script:
        """
        link_1=\$(ffq --ftp $sample_id | jq -r '.[0] | .url')
        if [ -z "\$link_1" ]; then
            echo "Error: failed to retrieve fastq1 link for $sample_id" >&2
            exit 1
        fi

        link_2=\$(ffq --ftp $sample_id | jq -r '.[1] | .url')
        if [ -z "\$link_2" ]; then
            echo "Error: failed to retrieve fastq2 link for $sample_id" >&2
            exit 1
        fi
        """
}

process Download {
    input:
        val sample_id
        val link
        val reads_to_save
        val read

    output:
        path "${sample_id}_R${read}.fastq.gz", emit: fq

    script:
        """
        if [[ "${reads_to_save}" == "all" ]]; then
            wget "${link}" -O "${sample_id}_R${read}.fastq.gz"
        else
            let lines_to_save=$reads_to_save*4
            wget -O - "${link}" | zcat | head -n "\${lines_to_save}" | gzip > "${sample_id}_R${read}.fastq.gz"
        fi
        """
}

process Downsample {
    input:
        path fastq
        val reads_to_save
        val read

    output:
        path "R${read}.fastq.gz", emit: fq

    script:
        """
        #!/bin/bash
        let lines_to_save=$reads_to_save*4
        zcat $fastq | head -\$lines_to_save | gzip > R${read}.fastq.gz
        """
}