process SEURAT_PROCESS {
    tag "scrna_velocity"
    label 'process_high'
    conda "${projectDir}/envs/seurat.yml"
    publishDir "${params.outdir}/scrna_velocity/seurat_process", mode: 'copy'

    input:
    path ctrl_10x_dir
    path treatment_loom
    path cell_cycle_genes

    output:
    path "seurat_process/integrated_final.rds",  emit: integrated_rds
    path "seurat_process/svg/*.svg",             emit: plots
    path "seurat_process/*.csv",                 emit: cluster_tables
    path "seurat_process/*_mqc.txt",             emit: mqc

    script:
    """
    seurat_process.R \\
        --ctrl_10x_dir ${ctrl_10x_dir} \\
        --ctrl_sample_name ${params.scrna_ctrl_sample_name} \\
        --treatment_loom ${treatment_loom} \\
        --treatment_sample_name ${params.scrna_treatment_sample_name} \\
        --cell_cycle_genes ${cell_cycle_genes} \\
        --n_g2m_genes ${params.scrna_n_g2m_genes} \\
        --n_s_genes ${params.scrna_n_s_genes} \\
        --mito_pattern '${params.scrna_mito_pattern}' \\
        --ribo_pattern '${params.scrna_ribo_pattern}' \\
        --ctrl_max_count ${params.scrna_ctrl_max_count} \\
        --ctrl_max_feature ${params.scrna_ctrl_max_feature} \\
        --ctrl_min_feature ${params.scrna_ctrl_min_feature} \\
        --ctrl_max_mito ${params.scrna_ctrl_max_mito} \\
        --treat_max_count ${params.scrna_treat_max_count} \\
        --treat_max_feature ${params.scrna_treat_max_feature} \\
        --treat_min_feature ${params.scrna_treat_min_feature} \\
        --treat_max_mito ${params.scrna_treat_max_mito} \\
        --npcs ${params.scrna_npcs} \\
        --dims ${params.scrna_dims} \\
        --nfeatures ${params.scrna_nfeatures} \\
        --ctrl_resolution ${params.scrna_ctrl_resolution} \\
        --ctrl_umap_neighbors ${params.scrna_ctrl_umap_neighbors} \\
        --ctrl_umap_min_dist ${params.scrna_ctrl_umap_min_dist} \\
        --ctrl_drop_clusters '${params.scrna_ctrl_drop_clusters}' \\
        --treat_resolution ${params.scrna_treat_resolution} \\
        --treat_umap_neighbors ${params.scrna_treat_umap_neighbors} \\
        --treat_umap_min_dist ${params.scrna_treat_umap_min_dist} \\
        --treat_drop_clusters '${params.scrna_treat_drop_clusters}' \\
        --anchor_features ${params.scrna_anchor_features} \\
        --anchor_dims ${params.scrna_anchor_dims} \\
        --integration1_resolution ${params.scrna_integration1_resolution} \\
        --integration1_umap_neighbors ${params.scrna_integration1_umap_neighbors} \\
        --integration1_umap_min_dist ${params.scrna_integration1_umap_min_dist} \\
        --integration1_drop_clusters '${params.scrna_integration1_drop_clusters}' \\
        --integration2_resolution ${params.scrna_integration2_resolution} \\
        --integration2_umap_neighbors ${params.scrna_integration2_umap_neighbors} \\
        --integration2_umap_min_dist ${params.scrna_integration2_umap_min_dist} \\
        --outdir seurat_process
    """

    stub:
    """
    mkdir -p seurat_process/svg
    touch seurat_process/integrated_final.rds
    touch seurat_process/svg/stub.svg
    touch seurat_process/integration2_cluster_by_sample.csv
    echo "# plot_type: 'html'" > seurat_process/seurat_process_summary_mqc.txt
    """
}
