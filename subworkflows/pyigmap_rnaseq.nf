include { FastpMockMerge } from '../modules/local/fastp.nf'
include { Vidjil } from '../modules/local/vidjil.nf'
include { IgBlastFASTA } from '../modules/local/igblast.nf'
include { CDR3ErrorCorrector } from '../modules/local/cdr3nt_error_corrector.nf'

workflow PYIGMAP_RNASEQ {
    take:
        fq1
        fq2

    main:
        FastpMockMerge(fq1, fq2)

        vidjil_ref = file(params.vidjil_ref)
        Vidjil(FastpMockMerge.out.fq12, vidjil_ref)

        igblast_ref = file(params.igblast_ref)
        IgBlastFASTA(Vidjil.out.fasta, igblast_ref)

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(IgBlastFASTA.out.annotation, olga_models, FastpMockMerge.out.json)
}