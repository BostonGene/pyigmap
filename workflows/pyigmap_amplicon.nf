include { CalibDedup } from '../steps/calib_dedup/calib_dedup.nf'
include { FastpMerge } from '../steps/fastp/fastp.nf'
include { IgBlastFASTQ } from '../steps/igblast/igblast.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt_error_corrector.nf'

workflow PYIGMAP_AMPLICON {
    take:
        fq1
        fq2

    main:
        CalibDedup(fq1, fq2)

        FastpMerge(CalibDedup.out.fq1, CalibDedup.out.fq2)

        igblast_ref = file(params.igblast_ref)
        IgBlastFASTQ(FastpMerge.out.fq1, FastpMerge.out.fq2, FastpMerge.out.fq12, igblast_ref)

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(IgBlastFASTQ.out.annotation, olga_models, CalibDedup.out.json)
}