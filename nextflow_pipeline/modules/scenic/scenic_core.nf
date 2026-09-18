process SCENIC_CORE {
    tag "scenic"
    label 'process_high'
    conda "${projectDir}/envs/scenic.yml"
    publishDir "${params.outdir}/scenic/core", mode: 'copy'

    input:
    path loom

    output:
    path "scenic/int/scenicOptions.Rds", emit: scenic_options
    path "scenic/int/*.Rds",             emit: intermediates
    path "scenic/*_mqc.txt",             emit: mqc

    script:
    """
    run_scenic_core.R \\
        --loom ${loom} \\
        --org ${params.scenic_org} \\
        --db_dir ${params.scenic_cistarget_dir} \\
        --db_10kb ${params.scenic_db_10kb} \\
        --db_500bp ${params.scenic_db_500bp} \\
        --dataset_title '${params.scenic_dataset_title}' \\
        --min_counts_per_cell ${params.scenic_min_counts_per_cell} \\
        --min_pct_cells ${params.scenic_min_pct_cells} \\
        --n_cores ${task.cpus} \\
        --seed ${params.scenic_seed} \\
        --outdir scenic
    """

    stub:
    """
    mkdir -p scenic/int
    touch scenic/int/scenicOptions.Rds
    touch scenic/int/exprMat.Rds
    echo "# plot_type: 'html'" > scenic/scenic_core_summary_mqc.txt
    """
}
