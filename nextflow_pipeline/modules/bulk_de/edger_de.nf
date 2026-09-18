process EDGER_DE {
    tag "bulk_de"
    label 'process_medium'
    conda "${projectDir}/envs/deseq2_edger.yml"
    publishDir "${params.outdir}/bulk_de/edger_de", mode: 'copy'

    input:
    path counts
    path sample_info

    output:
    path "edger_de/DGEgenes_*.csv", emit: degenes
    path "edger_de/*_mqc.txt",      emit: mqc

    script:
    """
    edger_de.R \\
        --counts ${counts} \\
        --sample_info ${sample_info} \\
        --sample_id_col ${params.bulk_sample_id_col} \\
        --contrast_col ${params.bulk_contrast_column} \\
        --contrast_levels ${params.bulk_contrast_levels} \\
        --min_cpm ${params.bulk_min_cpm} \\
        --min_samples ${params.bulk_min_samples} \\
        --lfc_threshold ${params.bulk_lfc_threshold} \\
        --fdr ${params.bulk_fdr} \\
        --outdir edger_de
    """

    stub:
    """
    mkdir -p edger_de
    touch edger_de/DGEgenes_stub.csv
    echo "# plot_type: 'html'" > edger_de/edger_de_summary_mqc.txt
    """
}
