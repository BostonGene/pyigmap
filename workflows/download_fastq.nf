include { GetLinks; Download as DownloadRead1; Download as DownloadRead2 } from '../steps/downloader/downloader.nf'

workflow DOWNLOAD_FASTQ {
    take:
        sample_id
        reads_to_save

    main:
        GetLinks(sample_id)
        DownloadRead1(sample_id, GetLinks.out.link_1, reads_to_save, Channel.from('1'))
        DownloadRead2(sample_id, GetLinks.out.link_2, reads_to_save, Channel.from('2'))

    emit:
        fq1 = DownloadRead1.out
        fq2 = DownloadRead2.out
}