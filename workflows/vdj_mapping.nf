include { Vidjil } from '../steps/vidjil/vidjil.nf'
include { IgBlast } from '../steps/igblast/igblast.nf'

workflow VDJ_MAPPING {
    take:
        fq12

    main:
        vidjil_ref = file(params.vidjil_ref)
        Vidjil(fq12, vidjil_ref)

        igblast_ref = file(params.igblast_ref)
        IgBlast(Vidjil.out.fasta, igblast_ref)

    emit:
        raw_annotation = IgBlast.out
}