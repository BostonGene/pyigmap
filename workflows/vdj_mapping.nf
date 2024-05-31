include { Vidjil } from '../steps/vidjil/vidjil.nf'
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
        Vidjil(fq1, fq2, fq12, vidjil_ref)

        IgBlastTCR(Vidjil.out.fasta, igblast_ref, 'TCR')
        IgBlastBCR(Vidjil.out.fasta, igblast_ref, 'BCR')

    emit:
        tcr_annotation = IgBlastTCR.out
        bcr_annotation = IgBlastBCR.out
}