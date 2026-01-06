include { FastpMockMerge; FastpSingle } from '../modules/local/fastp.nf'
include { Vidjil } from '../modules/local/vidjil.nf'
include { IgBlastFASTA } from '../modules/local/igblast.nf'
include { CDR3ErrorCorrector } from '../modules/local/cdr3nt_error_corrector.nf'

workflow PYIGMAP_RNASEQ {
    take:
        fq1
        fq2

    main:

        vidjil_ref   = file(params.vidjil_ref)
        igblast_ref  = file(params.igblast_ref)
        olga_models  = file(params.olga_models)

        if (params.paired) {

            FastpMockMerge(fq1, fq2)
            Vidjil(FastpMockMerge.out.fq12, vidjil_ref)
            IgBlastFASTA(Vidjil.out.fasta, igblast_ref)
            CDR3ErrorCorrector(
                IgBlastFASTA.out.annotation,
                olga_models,
                FastpMockMerge.out.json
            )

        } else {

            FastpSingle(fq1)
            Vidjil(FastpSingle.out.fq1, vidjil_ref)
            IgBlastFASTA(Vidjil.out.fasta, igblast_ref)
            CDR3ErrorCorrector(
                IgBlastFASTA.out.annotation,
                olga_models,
                FastpSingle.out.json
            )
        }
}