#!/usr/bin/env nextflow

include { PYIGMAP_AMPLICON } from './workflows/pyigmap_amplicon.nf'


params.vidjil_ref = './steps/vidjil/vidjil.germline.tar.gz'
params.igblast_ref = './steps/igblast/igblast.reference.tar.gz'
params.olga_models = './steps/cdr3nt-error-corrector/olga-models.tar.gz'

log.info ""
log.info "                     P Y I G M A P                     "
log.info "======================================================="
log.info "Mode                 : ${params.mode}"
log.info "Fastq1               : ${params.fq1}"
log.info "Fastq2               : ${params.fq2}"
log.info "Output directory     : ${params.outdir}"
log.info "Vidjil reference     : ${params.vidjil_ref}"
log.info "IgBLAST reference    : ${params.igblast_ref}"
log.info "OLGA models          : ${params.olga_models}"
log.info ""


workflow {
    fq1 = file(params.fq1)
    fq2 = file(params.fq2)
    igblast_ref = file(params.igblast_ref)
    vidjil_ref = file(params.vidjil_ref)
    olga_models = file(params.olga_models)
    if (params.mode == "amplicon") {
        PYIGMAP_AMPLICON(fq1, fq2, vidjil_ref, igblast_ref, olga_models)
    } else if (params.mode == "rnaseq") {
    } else {
        error "Unexpected --mode '${params.mode}'. Please make sure to choose either --mode 'amplicon' or --mode 'rnaseq'."
    }
}
