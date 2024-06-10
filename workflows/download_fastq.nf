include { GetLinks; Download as DownloadRead1; Download as DownloadRead2 } from '../steps/downloader/downloader.nf'

workflow DOWNLOAD_FASTQ_BY_SAMPLE_ID {
    take:
        sample_id
        reads_to_process

    main:
        GetLinks(sample_id)
        DownloadRead1(sample_id, GetLinks.out.link_1, reads_to_process, Channel.from('1'))
        DownloadRead2(sample_id, GetLinks.out.link_2, reads_to_process, Channel.from('2'))

    emit:
        fq1 = DownloadRead1.out
        fq2 = DownloadRead2.out
}

workflow DOWNLOAD_FASTQ_BY_LINK {
    take:
        link_1
        link_2
        sample_id
        reads_to_process

    main:
        DownloadRead1(sample_id, link_1, reads_to_process, Channel.from('1'))
        DownloadRead2(sample_id, link_2, reads_to_process, Channel.from('2'))

    emit:
        fq1 = DownloadRead1.out
        fq2 = DownloadRead2.out
}