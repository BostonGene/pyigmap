#!/usr/bin/env nextflow

include { PYIGMAP_AMPLICON } from './workflows/pyigmap_amplicon.nf'
include { PYIGMAP_RNASEQ } from './workflows/pyigmap_rnaseq.nf'
include { DOWNLOAD_FASTQ_BY_SAMPLE_ID; DOWNLOAD_FASTQ_BY_LINK } from './workflows/download_fastq.nf'
include { Downsample as DownsampleRead1; Downsample as DownsampleRead2 } from './steps/downloader/downloader.nf'

params.zenodo = false
params.vidjil_ref = './steps/vidjil/vidjil.germline.tar.gz'
params.igblast_ref = './steps/igblast/igblast.reference.tar.gz'
params.olga_models = './steps/cdr3nt_error_corrector/olga-models.tar.gz'
params.reads = 'all'

log.info ""
log.info "                     P Y I G M A P                     "
log.info "======================================================="
log.info "Mode                 : ${params.mode}"
log.info "Sample               : ${params.sample}"
log.info "Enable zenodo        : ${params.zenodo}"
log.info "Reads to process     : ${params.reads}"
log.info "Fastq1               : ${params.fq1}"
log.info "Fastq2               : ${params.fq2}"
log.info "Output directory     : ${params.outdir}"
log.info "Vidjil reference     : ${params.vidjil_ref}"
log.info "IgBLAST reference    : ${params.igblast_ref}"
log.info "OLGA models          : ${params.olga_models}"
log.info ""

workflow {
    reads_to_save = Channel.from(params.reads)
    if (params.sample) {
        if (params.fq1 || params.fq2) {
            error "Error: --sample cannot be used with --fq1 and --fq2 at the same time, exiting..."
        }
        ArrayList sample = params.sample.split(",")
        sample_ch = Channel.from(sample)
        DOWNLOAD_FASTQ_BY_SAMPLE_ID(sample_ch, reads_to_save)
        fq1 = DOWNLOAD_FASTQ_BY_SAMPLE_ID.out.fq1
        fq2 = DOWNLOAD_FASTQ_BY_SAMPLE_ID.out.fq2
    } else {
        if (!params.fq1 || !params.fq2) {
            error "Error: single-end is not supported, exiting..."
        }
        if (params.zenodo) {
            fq1_link = params.zenodo_link + params.fq1
            fq2_link = params.zenodo_link + params.fq2
            sample = params.fq1.replace("_R1.fastq.gz", "").replace("_1.fastq.gz", "")
            DOWNLOAD_FASTQ_BY_LINK(fq1_link, fq2_link, sample, reads_to_save)
            fq1 = DOWNLOAD_FASTQ_BY_LINK.out.fq1
            fq2 = DOWNLOAD_FASTQ_BY_LINK.out.fq2
        } else {
            fq1 = Channel.fromPath(params.fq1)
            fq2 = Channel.fromPath(params.fq2)
        }
        if (params.reads.toString().isInteger()) {
            DownsampleRead1(fq1, reads_to_save, Channel.from("1"))
            DownsampleRead2(fq2, reads_to_save, Channel.from("2"))
            fq1 = DownsampleRead1.out.fq
            fq2 = DownsampleRead2.out.fq
        }
    }
    igblast_ref = file(params.igblast_ref)
    vidjil_ref = file(params.vidjil_ref)
    olga_models = file(params.olga_models)
    if (params.mode == "amplicon") {
        PYIGMAP_AMPLICON(fq1, fq2, vidjil_ref, igblast_ref, olga_models)
    } else if (params.mode == "rnaseq") {
        PYIGMAP_RNASEQ(fq1, fq2, vidjil_ref, igblast_ref, olga_models)
    } else {
        error "Error: unexpected --mode '${params.mode}'. Supported only two modes: --mode 'amplicon' and --mode 'rnaseq'."
    }
}
