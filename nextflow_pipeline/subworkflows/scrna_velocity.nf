include { SEURAT_PROCESS }         from '../modules/scrna_velocity/seurat_process.nf'
include { SEURAT_VELOCITY_EXPORT } from '../modules/scrna_velocity/seurat_velocity_export.nf'
include { SCVELO_ANALYSIS }        from '../modules/scrna_velocity/scvelo_analysis.nf'

// 2-condition (control vs. treatment) scRNA-seq QC/clustering/integration,
// with an optional RNA-velocity leg. The velocity leg is decoupled from the
// integration above (matching the original eLife scripts, where the
// velocity input is a separately prepared, already-subset Seurat object) --
// it only runs if velocity_rds is supplied.
workflow SCRNA_VELOCITY {
    take:
    ctrl_10x_dir
    treatment_loom
    cell_cycle_genes
    velocity_rds   // optional (pass Channel.empty() to skip the velocity leg)

    main:
    mqc_ch = Channel.empty()

    SEURAT_PROCESS(ctrl_10x_dir, treatment_loom, cell_cycle_genes)
    mqc_ch = mqc_ch.mix(SEURAT_PROCESS.out.mqc)

    velocity_h5ad = Channel.empty()

    if (params.scrna_velocity_rds) {
        SEURAT_VELOCITY_EXPORT(velocity_rds)
        SCVELO_ANALYSIS(SEURAT_VELOCITY_EXPORT.out.h5ad)
        mqc_ch = mqc_ch.mix(SCVELO_ANALYSIS.out.mqc)
        velocity_h5ad = SCVELO_ANALYSIS.out.h5ad
    }

    emit:
    integrated_rds = SEURAT_PROCESS.out.integrated_rds
    velocity_h5ad
    mqc            = mqc_ch
}
