include { Vidjil } from '../steps/vidjil/vidjil.nf'
include { IgBlast } from '../steps/igblast/igblast.nf'

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
        IgBlast(Vidjil.out.fasta, igblast_ref)

    emit:
        raw_annotation = IgBlast.out
}