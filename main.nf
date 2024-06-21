#!/usr/bin/env nextflow

include { PYIGMAP_AMPLICON } from './workflows/pyigmap_amplicon.nf'
include { PYIGMAP_RNASEQ } from './workflows/pyigmap_rnaseq.nf'
include { DOWNLOAD_FASTQ_BY_SAMPLE_ID; DOWNLOAD_FASTQ_BY_LINK; DOWNSAMPLE_FASTQ } from './workflows/download_fastq.nf'


if (params.help) { exit 0, help_message() }

log.info ""
log.info "                     P Y I G M A P                     "
log.info "======================================================="
log.info "Library type         : ${params.library}"
log.info "Sample               : ${params.sample_id}"
log.info "Reads to process     : ${params.reads_to_process}"
log.info "All alleles          : ${params.all_alleles}"
log.info "Forward FASTQ        : ${params.fq1}"
log.info "Reverse FASTQ        : ${params.fq2}"
log.info "Output directory     : ${params.outdir}"
log.info "Vidjil reference     : ${params.vidjil_ref}"
log.info "IgBLAST reference    : ${params.igblast_ref}"
log.info "OLGA models          : ${params.olga_models}"
log.info ""

def help_message() {
    log.info """
    pyIgMap is a workflow for extracting and summarizing antigen receptor gene rearrangements from sequencing data.

        Usage example:
    1. To process AIRR-Seq target data with UMI in the first twelve letters of the forward read
    ./pyigmap --library amplicon --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz --fq1_pattern "^UMI:N{12}"

    2. To process RNASeq-bulk data
    ./pyigmap --library rnaseq --fq1 /path/to/R1.fastq.gz --fq2 /path/to/R2.fastq.gz

    3. To process public RNASeq data by sample id
    ./pyigmap --library rnaseq --sample_id SRR3743469 --reads_to_process 200000

    4. To process RNASeq data downloaded by HTTP/HTTPS/FTP link
    ./pyigmap --library rnaseq --fq1 https://zenodo.org/records/11103555/files/SRR3743469_R1.fastq.gz --fq2 https://zenodo.org/records/11103555/files/SRR3743469_R1.fastq.gz --reads_to_process 200000

        Optional input:
    --fq1                       path to the forward FASTQ (default: none)
    --fq2                       path to the reverse FASTQ (default: none)

    or

    --sample_id                    sample id name (default: none)

        Output:
    --outdir                    path to the output directory (default: ${params.outdir})

        Output file:

    pyigmap.tar.gz archive that contains:
    * corrected_annotation.tsv  corrected AIRR-formatted annotation
    * stat.json                 file with common statistics


        Workflow Options:
    --library                   the library type of input data: "rnaseq" (RNASeq-bulk) or "amplicon" (AIRR-Seq target)
    --all_alleles               will use all alleles provided in the antigen receptor segment database (*01, *02, etc. according to IMGT);
                                only major (*01) allele for each gene will be used otherwise (default: false)

        Nextflow options:
    -resume                     resume the workflow where it stopped
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)
    """
}

workflow {
    if (params.sample_id && !params.fq1 && !params.fq2) {

        DOWNLOAD_FASTQ_BY_SAMPLE_ID(params.sample_id)

        fq1 = DOWNLOAD_FASTQ_BY_SAMPLE_ID.out.fq1
        fq2 = DOWNLOAD_FASTQ_BY_SAMPLE_ID.out.fq2
    } else {

        if (!params.fq1 || !params.fq2) {
            error "Error: single-end is not supported, exiting..."
        }

        def urlPatterns = ["ftp://", "http://", "https://"]

        if (urlPatterns.any { params.fq1.startsWith(it) } || urlPatterns.any { params.fq2.startsWith(it) }) {
            DOWNLOAD_FASTQ_BY_LINK(params.fq1, params.fq2)
            fq1 = DOWNLOAD_FASTQ_BY_LINK.out.fq1
            fq2 = DOWNLOAD_FASTQ_BY_LINK.out.fq2
        } else if (params.reads_to_process.toString().isInteger()) {
            DOWNSAMPLE_FASTQ(params.fq1, params.fq2)
            fq1 = DOWNSAMPLE_FASTQ.out.fq1
            fq2 = DOWNSAMPLE_FASTQ.out.fq2
        }
    }

    if (params.library == "amplicon") {
        PYIGMAP_AMPLICON(fq1, fq2)
    } else if (params.library == "rnaseq") {
        PYIGMAP_RNASEQ(fq1, fq2)
    } else {
        error "Error: unexpected --library '${params.library}'. Supported only two libraries: 'amplicon' and 'rnaseq'."
    }
}

workflow.onComplete {
  log.info ( workflow.success ? "\nDone! Results are stored here --> ${params.outdir}/${params.out_archive} \n": "Oops .. something went wrong" )
}