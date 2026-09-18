include { SCENIC_CORE }       from '../modules/scenic/scenic_core.nf'
include { SCENIC_DOWNSTREAM } from '../modules/scenic/scenic_downstream.nf'

workflow SCENIC {
    take:
    loom         // expression matrix as .loom
    seurat_rds   // Seurat object supplying the RSS/CSI grouping column; only needed if run_scenic_downstream

    main:
    mqc_ch = Channel.empty()

    SCENIC_CORE(loom)
    mqc_ch = mqc_ch.mix(SCENIC_CORE.out.mqc)

    if (params.run_scenic_downstream) {
        SCENIC_DOWNSTREAM(SCENIC_CORE.out.scenic_options, seurat_rds)
        mqc_ch = mqc_ch.mix(SCENIC_DOWNSTREAM.out.mqc)
    }

    emit:
    scenic_options = SCENIC_CORE.out.scenic_options
    mqc            = mqc_ch
}
