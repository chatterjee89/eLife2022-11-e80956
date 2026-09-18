process MULTIQC {
    tag "report"
    label 'process_low'
    conda "${projectDir}/envs/multiqc.yml"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path mqc_files, stageAs: "mqc_inputs/*"

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc mqc_inputs -n multiqc_report.html
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    """
}
