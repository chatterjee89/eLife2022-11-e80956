include { DESEQ2_DIAGNOSTICS } from '../modules/bulk_de/deseq2_diagnostics.nf'
include { EDGER_DE }           from '../modules/bulk_de/edger_de.nf'

workflow BULK_DE {
    take:
    counts        // path: featureCounts-style matrix CSV
    sample_info   // path: sample metadata CSV

    main:
    mqc_ch = Channel.empty()

    EDGER_DE(counts, sample_info)
    mqc_ch = mqc_ch.mix(EDGER_DE.out.mqc)

    if (!params.skip_deseq2_diagnostics) {
        DESEQ2_DIAGNOSTICS(counts, sample_info)
        mqc_ch = mqc_ch.mix(DESEQ2_DIAGNOSTICS.out.mqc)
    }

    emit:
    degenes = EDGER_DE.out.degenes
    mqc     = mqc_ch
}
