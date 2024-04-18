include { Fastp } from '../steps/fastp/fastp.nf'
include { VDJ_MAPPING } from './vdj-mapping.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt-error-corrector.nf'

workflow PYIGMAP_RNASEQ {
    take:
        fq1
        fq2
        vidjil_ref
        igblast_ref
        olga_models

    main:
        Fastp(fq1, fq2)
        VDJ_MAPPING(Fastp.out.fq1, Fastp.out.fq2, Fastp.out.fq12, vidjil_ref, igblast_ref, olga_models)
        CDR3ErrorCorrector(VDJ_MAPPING.out.tcr_annotation, VDJ_MAPPING.out.bcr_annotation, olga_models, Fastp.out.json)
}