process SCVELO_ANALYSIS {
    tag "scrna_velocity"
    label 'process_high'
    conda "${projectDir}/envs/scvelo.yml"
    publishDir "${params.outdir}/scrna_velocity/scvelo", mode: 'copy'

    input:
    path h5ad

    output:
    path "scvelo_analysis/*.svg",                  emit: plots
    path "scvelo_analysis/velocity_processed.h5ad", emit: h5ad
    path "scvelo_analysis/*_mqc.txt",              emit: mqc

    script:
    """
    scvelo_analysis.py \\
        --h5ad ${h5ad} \\
        --groupby ${params.scvelo_groupby} \\
        --min_r2 ${params.scvelo_min_r2} \\
        --min_corr_rank_genes ${params.scvelo_min_corr_rank_genes} \\
        --top_n_genes ${params.scvelo_top_n_genes} \\
        --outdir scvelo_analysis
    """

    stub:
    """
    mkdir -p scvelo_analysis
    touch scvelo_analysis/stub.svg
    touch scvelo_analysis/velocity_processed.h5ad
    echo "# plot_type: 'html'" > scvelo_analysis/scvelo_summary_mqc.txt
    """
}
