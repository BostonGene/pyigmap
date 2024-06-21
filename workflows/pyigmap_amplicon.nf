include { PyUMI } from '../steps/pyumi/pyumi.nf'
include { CalibDedup } from '../steps/calib_dedup/calib_dedup.nf'
include { FastpMerge; FastpMockMerge } from '../steps/fastp/fastp.nf'
include { IgBlastFASTQ; IgBlastMockFASTQ } from '../steps/igblast/igblast.nf'
include { CDR3ErrorCorrector } from '../steps/cdr3nt_error_corrector/cdr3nt_error_corrector.nf'

workflow PYIGMAP_AMPLICON {
    take:
        fq1
        fq2

    main:
        if (params.fq1_pattern == null && params.fq2_pattern == null) {
            PYIGMAP_AMPLICON_WITHOUT_UMI(fq1, fq2)
        } else {
            PYIGMAP_AMPLICON_WITH_UMI(fq1, fq2)
        }
}

workflow PYIGMAP_AMPLICON_WITH_UMI {
    take:
        fq1
        fq2

    main:
        PyUMI(fq1, fq2)
        CalibDedup(PyUMI.out.fq1, PyUMI.out.fq2, PyUMI.out.json)

        igblast_ref = file(params.igblast_ref)

        if (params.mock_merge_amplicon) {
            FastpMockMerge(CalibDedup.out.fq1, CalibDedup.out.fq2)
            IgBlastMockFASTQ(FastpMockMerge.out.fq12, igblast_ref)
            raw_annotation = IgBlastMockFASTQ.out.annotation
        } else {
            FastpMerge(CalibDedup.out.fq1, CalibDedup.out.fq2)
            IgBlastFASTQ(FastpMerge.out.fq1, FastpMerge.out.fq2, FastpMerge.out.fq12, igblast_ref)
            raw_annotation = IgBlastFASTQ.out.annotation
        }

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(raw_annotation, olga_models, PyUMI.out.json)
}

workflow PYIGMAP_AMPLICON_WITHOUT_UMI {
    take:
        fq1
        fq2

    main:
        igblast_ref = file(params.igblast_ref)

        if (params.mock_merge_amplicon) {
            FastpMockMerge(fq1, fq2)
            IgBlastMockFASTQ(FastpMockMerge.out.fq12, igblast_ref)
            raw_annotation = IgBlastMockFASTQ.out.annotation
        } else {
            FastpMerge(fq1, fq2)
            IgBlastFASTQ(FastpMerge.out.fq1, FastpMerge.out.fq2, FastpMerge.out.fq12, igblast_ref)
            raw_annotation = IgBlastFASTQ.out.annotation
        }

        olga_models = file(params.olga_models)
        CDR3ErrorCorrector(raw_annotation, olga_models, FastpMerge.out.json)
}