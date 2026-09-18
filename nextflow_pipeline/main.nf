#!/usr/bin/env nextflow
/*
 * nf-multiomics-de-grn
 * A generalized Nextflow DSL2 pipeline combining:
 *   - bulk RNA-seq differential expression (DESeq2 diagnostics + edgeR DE),
 *     ported from eLife2022-11-e80956/bulkSeq_code.R
 *   - scRNA-seq QC/integration + RNA velocity (Seurat + scVelo),
 *     ported from eLife2022-11-e80956/seurat_code.R and scvelo_code.py
 *   - R-SCENIC gene regulatory network inference (+ optional RSS/CSI
 *     downstream analysis), ported from
 *     BPTF_mammary_GRN_analyses/02_SCENIC/{01_run_SCENIC.R,02_SCENIC_downstream_analysis.R}
 *
 * All three started as separate, dataset-specific analyses; here they're
 * exposed as independent subworkflows selected via --mode, sharing one
 * config/execution/reporting layer. Pre-processed inputs only -- this
 * pipeline does not run alignment (STAR) or velocyto itself.
 *
 * AI usage note: this pipeline's structure and Nextflow wiring were written
 * by Claude (Anthropic), generalizing R/Python analysis scripts the author
 * originally wrote per-dataset; reviewed by the author. See AI_DISCLOSURE.md.
 */

nextflow.enable.dsl = 2

include { BULK_DE }        from './subworkflows/bulk_de.nf'
include { SCRNA_VELOCITY }  from './subworkflows/scrna_velocity.nf'
include { SCENIC }          from './subworkflows/scenic.nf'
include { MULTIQC }         from './modules/report/multiqc.nf'

def VALID_MODES = ['bulk_de', 'scrna_velocity', 'scenic', 'full']

workflow {

    if (!params.mode || !(params.mode in VALID_MODES)) {
        error "Please set --mode to one of: ${VALID_MODES.join(', ')} (got: ${params.mode})"
    }

    mqc_ch = Channel.empty()
    run_bulk   = params.mode in ['bulk_de', 'full']
    run_scrna  = params.mode in ['scrna_velocity', 'full']
    run_scenic = params.mode in ['scenic', 'full']

    // -------------------------------------------------------------------
    // BULK RNA-SEQ DE
    // -------------------------------------------------------------------
    if (run_bulk) {
        if (!params.bulk_counts || !params.bulk_sample_info)
            error "--mode ${params.mode} requires --bulk_counts and --bulk_sample_info"

        BULK_DE(
            Channel.fromPath(params.bulk_counts, checkIfExists: true),
            Channel.fromPath(params.bulk_sample_info, checkIfExists: true)
        )
        mqc_ch = mqc_ch.mix(BULK_DE.out.mqc)
    }

    // -------------------------------------------------------------------
    // SCRNA-SEQ + VELOCITY (2-condition: control vs. treatment)
    // -------------------------------------------------------------------
    scrna_integrated_rds = Channel.empty()
    if (run_scrna) {
        if (!params.scrna_ctrl_10x_dir || !params.scrna_treatment_loom || !params.scrna_cell_cycle_genes)
            error "--mode ${params.mode} requires --scrna_ctrl_10x_dir, --scrna_treatment_loom, and --scrna_cell_cycle_genes"

        velocity_rds_ch = params.scrna_velocity_rds
            ? Channel.fromPath(params.scrna_velocity_rds, checkIfExists: true)
            : Channel.empty()

        SCRNA_VELOCITY(
            Channel.fromPath(params.scrna_ctrl_10x_dir, checkIfExists: true),
            Channel.fromPath(params.scrna_treatment_loom, checkIfExists: true),
            Channel.fromPath(params.scrna_cell_cycle_genes, checkIfExists: true),
            velocity_rds_ch
        )
        mqc_ch = mqc_ch.mix(SCRNA_VELOCITY.out.mqc)
        scrna_integrated_rds = SCRNA_VELOCITY.out.integrated_rds
    }

    // -------------------------------------------------------------------
    // R-SCENIC GRN INFERENCE (+ optional RSS/CSI downstream)
    // -------------------------------------------------------------------
    if (run_scenic) {
        if (!params.scenic_loom)
            error "--mode ${params.mode} requires --scenic_loom"
        if (!params.scenic_cistarget_dir || !params.scenic_db_10kb || !params.scenic_db_500bp)
            error "--mode ${params.mode} requires --scenic_cistarget_dir, --scenic_db_10kb, and --scenic_db_500bp"

        // Downstream RSS/CSI needs a Seurat object for its cell-type/cluster
        // column. Prefer an explicit --scenic_seurat_rds; in `full` mode,
        // fall back to the scrna_velocity subworkflow's integrated object
        // (same run, plausibly the same cells) if one isn't given.
        seurat_rds_for_scenic = Channel.empty()
        if (params.run_scenic_downstream) {
            if (params.scenic_seurat_rds) {
                seurat_rds_for_scenic = Channel.fromPath(params.scenic_seurat_rds, checkIfExists: true)
            } else if (params.mode == 'full') {
                seurat_rds_for_scenic = scrna_integrated_rds
            } else {
                error "--run_scenic_downstream requires --scenic_seurat_rds when --mode is not 'full'"
            }
        }

        SCENIC(
            Channel.fromPath(params.scenic_loom, checkIfExists: true),
            seurat_rds_for_scenic
        )
        mqc_ch = mqc_ch.mix(SCENIC.out.mqc)
    }

    // -------------------------------------------------------------------
    // REPORT
    // -------------------------------------------------------------------
    if (params.run_multiqc) {
        MULTIQC(mqc_ch.collect())
    }
}
