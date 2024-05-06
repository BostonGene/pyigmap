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
if (params.help) { exit 0, help_message() }

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

def help_message() {
    log.info """
    pyIgMap is a workflow for extracting and summarizing antigen receptor gene rearrangements from sequencing data.

        Usage example:
    1. For amplicon data
    ./pyigmap --mode amplicon --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz

    2. For RNASeq data
    ./pyigmap --mode rnaseq --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz

    3. For public data by sample id
    ./pyigmap --mode rnaseq --sample SRR3743469 --reads 200000

    4. For public data from ZENODO
    ./pyigmap --mode rnaseq --zenodo --fq1 SRR3743469_R1.fastq.gz --fq2 SRR3743469_R2.fastq.gz --reads 200000

        Optional input:
    --fq1                       path to the read1 fastq (default: none)
    --fq2                       path to the read2 fastq (default: none)
    --zenodo                    enables downloading from ZENODO (default: false)

    or

    --sample                    sample id name (default: none)

        Output:
    --outdir                    path to the output directory (default: ${params.outdir})

        Output file:

    pyigmap.tar.gz archive that contains:
    * corrected_annotation.tsv  the corrected annotation in AIRR-format
    * stat.json                 the file with different statistics


        Workflow Options:
    --mode                      the mode of pipeline "rnaseq" (without umi-preprocessing) or "amplicon" (with umi-preprocessing)

        Nextflow options:
    -resume                     resume the workflow where it stopped
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)
    """
}
workflow {
    reads_to_save = Channel.from(params.reads)

    // Running pipeline by sample id
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

workflow.onComplete {
  log.info ( workflow.success ? "\nDone! Results are stored here --> ${params.outdir} \n": "Oops .. something went wrong" )
}