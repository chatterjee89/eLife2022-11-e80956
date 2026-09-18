process SCENIC_DOWNSTREAM {
    tag "scenic"
    label 'process_medium'
    conda "${projectDir}/envs/scenic.yml"
    publishDir "${params.outdir}/scenic/downstream", mode: 'copy'

    input:
    path scenic_options_rds
    path seurat_rds

    output:
    path "scenic_downstream/svg/*.svg", emit: plots
    path "scenic_downstream/*.csv",     emit: tables
    path "scenic_downstream/*_mqc.txt", emit: mqc

    script:
    """
    scenic_downstream_analysis.R \\
        --scenic_options_rds ${scenic_options_rds} \\
        --seurat_rds ${seurat_rds} \\
        --cell_type_column ${params.scenic_cell_type_column} \\
        --rss_threshold ${params.scenic_rss_threshold} \\
        --csi_nclust ${params.scenic_csi_nclust} \\
        --outdir scenic_downstream
    """

    stub:
    """
    mkdir -p scenic_downstream/svg
    touch scenic_downstream/svg/stub.svg
    touch scenic_downstream/regulon_specificity_score.csv
    echo "# plot_type: 'html'" > scenic_downstream/scenic_downstream_summary_mqc.txt
    """
}
