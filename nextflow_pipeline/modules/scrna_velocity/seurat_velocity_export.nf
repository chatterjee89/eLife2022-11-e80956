process SEURAT_VELOCITY_EXPORT {
    tag "scrna_velocity"
    label 'process_medium'
    conda "${projectDir}/envs/seurat.yml"
    publishDir "${params.outdir}/scrna_velocity/velocity_export", mode: 'copy'

    input:
    path velocity_rds

    output:
    path "velocity.loom", emit: loom
    path "velocity.h5ad", emit: h5ad

    script:
    """
    seurat_velocity_export.R \\
        --velocity_rds ${velocity_rds} \\
        --assays '${params.scrna_velocity_export_assays}' \\
        --loom_out velocity.loom \\
        --h5ad_out velocity.h5ad
    """

    stub:
    """
    touch velocity.loom
    touch velocity.h5ad
    """
}
