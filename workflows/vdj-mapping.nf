include { Vidjil as VidjilBCR } from '../steps/vidjil/vidjil.nf'
include { Vidjil as VidjilTCR } from '../steps/vidjil/vidjil.nf'
include { IgBlast as IgBlastBCR } from '../steps/igblast/igblast.nf'
include { IgBlast as IgBlastTCR } from '../steps/igblast/igblast.nf'

workflow VDJ_MAPPING {
    take:
        fq1
        fq2
        fq12
        vidjil_ref
        igblast_ref
        olga_models

    main:
        VidjilTCR(fq1, fq2, fq12, vidjil_ref, 'TCR')
        VidjilBCR(fq1, fq2, fq12, vidjil_ref, 'BCR')

        IgBlastTCR(VidjilTCR.out.fasta, igblast_ref, 'TCR')
        IgBlastBCR(VidjilBCR.out.fasta, igblast_ref, 'BCR')

    emit:
        tcr_annotation = IgBlastTCR.out
        bcr_annotation = IgBlastBCR.out
}