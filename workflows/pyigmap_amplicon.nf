include { CalibDedup } from '../steps/calib_dedup/calib_dedup.nf'
include { Fastp } from '../steps/fastp/fastp.nf'
include { VDJ_MAPPING } from './vdj_mapping.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt_error_corrector.nf'

workflow PYIGMAP_AMPLICON {
    take:
        fq1
        fq2
        vidjil_ref
        igblast_ref
        olga_models

    main:
        CalibDedup(fq1, fq2)
        Fastp(CalibDedup.out.fq1, CalibDedup.out.fq2)
        VDJ_MAPPING(Fastp.out.fq1, Fastp.out.fq2, Fastp.out.fq12, vidjil_ref, igblast_ref, olga_models)
        CDR3ErrorCorrector(VDJ_MAPPING.out.raw_annotation, olga_models, CalibDedup.out.json)
}