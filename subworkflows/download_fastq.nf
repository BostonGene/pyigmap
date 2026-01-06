include { GetLinks; Download as DownloadRead1; Download as DownloadRead2 } from '../modules/local/downloader.nf'
include { Downsample as DownsampleRead1; Downsample as DownsampleRead2 } from '../modules/local/downloader.nf'

workflow DOWNLOAD_FASTQ_BY_LINK {
    take:
        link_1
        link_2

    main:
        DownloadRead1(params.sample_id, link_1, Channel.of('1'))

        if (params.paired) {
            DownloadRead2(params.sample_id, link_2, Channel.of('2'))
        }

    emit:
        fq1 = DownloadRead1.out
        fq2 = params.paired ? DownloadRead2.out : Channel.empty()
}

workflow DOWNLOAD_FASTQ_BY_SAMPLE_ID {
    take:
        sample_id

    main:
        GetLinks(sample_id)

        if (params.single) {
            // только R1
            DOWNLOAD_FASTQ_BY_LINK(GetLinks.out.link_1, "")
        } else {
            // R1 + R2
            DOWNLOAD_FASTQ_BY_LINK(GetLinks.out.link_1, GetLinks.out.link_2)
        }

    emit:
        fq1 = DOWNLOAD_FASTQ_BY_LINK.out.fq1
        fq2 = DOWNLOAD_FASTQ_BY_LINK.out.fq2
}

workflow DOWNSAMPLE_FASTQ {
    take:
        fq1
        fq2
    main:
        DownsampleRead1(fq1, Channel.from("1"))
        DownsampleRead2(fq2, Channel.from("2"))

    emit:
        fq1 = DownsampleRead1.out
        fq2 = DownsampleRead2.out
}