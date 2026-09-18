process DESEQ2_DIAGNOSTICS {
    tag "bulk_de"
    label 'process_medium'
    conda "${projectDir}/envs/deseq2_edger.yml"
    publishDir "${params.outdir}/bulk_de/deseq2_diagnostics", mode: 'copy'

    input:
    path counts
    path sample_info

    output:
    path "deseq2_diagnostics/svg/*.svg",              emit: plots
    path "deseq2_diagnostics/normalized_count_stats.csv", emit: stats
    path "deseq2_diagnostics/*_mqc.txt",               emit: mqc

    script:
    """
    deseq2_diagnostics.R \\
        --counts ${counts} \\
        --sample_info ${sample_info} \\
        --sample_id_col ${params.bulk_sample_id_col} \\
        --group_columns ${params.deseq2_group_columns} \\
        --alpha ${params.deseq2_alpha} \\
        --outdir deseq2_diagnostics
    """

    stub:
    """
    mkdir -p deseq2_diagnostics/svg
    touch deseq2_diagnostics/svg/stub.svg
    touch deseq2_diagnostics/normalized_count_stats.csv
    echo "# plot_type: 'html'" > deseq2_diagnostics/deseq2_diagnostics_summary_mqc.txt
    """
}
