include { CalibDedup } from '../steps/calib-dedup/calib-dedup.nf'
include { Fastp } from '../steps/fastp/fastp.nf'
include { Vidjil as VidjilBCR } from '../steps/vidjil/vidjil.nf'
include { Vidjil as VidjilTCR } from '../steps/vidjil/vidjil.nf'
include { IgBlast as IgBlastBCR } from '../steps/igblast/igblast.nf'
include { IgBlast as IgBlastTCR } from '../steps/igblast/igblast.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt-error-corrector/cdr3nt-error-corrector.nf'

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

        VidjilTCR(Fastp.out.fq1, Fastp.out.fq2, Fastp.out.fq12, vidjil_ref, 'TCR')
        VidjilBCR(Fastp.out.fq1, Fastp.out.fq2, Fastp.out.fq12, vidjil_ref, 'BCR')

        IgBlastTCR(VidjilTCR.out.fasta, igblast_ref, 'TCR')
        IgBlastBCR(VidjilBCR.out.fasta, igblast_ref, 'BCR')

        CDR3ErrorCorrector(IgBlastBCR.out.annotation, IgBlastTCR.out.annotation, olga_models, CalibDedup.out.json)
}