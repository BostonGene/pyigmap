include { Fastp } from '../steps/fastp/fastp.nf'
include { Vidjil } from '../steps/vidjil/vidjil.nf'
include { IgBlastFASTA } from '../steps/igblast/igblast.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt_error_corrector.nf'

workflow PYIGMAP_RNASEQ {
    take:
        fq1
        fq2

    main:
        Fastp(fq1, fq2)

        vidjil_ref = file(params.vidjil_ref)
        Vidjil(Fastp.out.fq1, Fastp.out.fq2, Fastp.out.fq12, vidjil_ref)

        igblast_ref = file(params.igblast_ref)
        IgBlastFASTA(Vidjil.out.fasta, igblast_ref)

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(IgBlastFASTA.out.annotation, olga_models, Fastp.out.json)
}