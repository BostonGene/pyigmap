include { GetLinks; Download as DownloadRead1; Download as DownloadRead2 } from '../steps/downloader/downloader.nf'

workflow DOWNLOAD_FASTQ_BY_SAMPLE_ID {
    main:
        GetLinks(params.sample_id)
        DownloadRead1(params.sample_id, GetLinks.out.link_1, Channel.from('1'))
        DownloadRead2(params.sample_id, GetLinks.out.link_2, Channel.from('2'))

    emit:
        fq1 = DownloadRead1.out
        fq2 = DownloadRead2.out
}

workflow DOWNLOAD_FASTQ_BY_LINK {
    main:
        fq1_link = params.zenodo_link + params.fq1
        fq2_link = params.zenodo_link + params.fq2
        sample_id = params.fq1.replace("_R1.fastq.gz", "").replace("_1.fastq.gz", "")
        DownloadRead1(sample_id, fq1_link, Channel.from('1'))
        DownloadRead2(sample_id, fq2_link, Channel.from('2'))

    emit:
        fq1 = DownloadRead1.out
        fq2 = DownloadRead2.out
}