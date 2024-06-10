include { Fastp } from '../steps/fastp/fastp.nf'
include { VDJ_MAPPING } from './vdj_mapping.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt_error_corrector.nf'

workflow PYIGMAP_RNASEQ {
    take:
        fq1
        fq2

    main:
        Fastp(fq1, fq2)
        VDJ_MAPPING(Fastp.out.fq12)

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(VDJ_MAPPING.out.raw_annotation, olga_models, Fastp.out.json)
}