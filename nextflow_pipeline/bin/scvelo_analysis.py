#!/usr/bin/env python
# RNA velocity analysis via scVelo, generalized from the eLife2022-11-e80956
# scvelo_code.py (steady-state pass -> confidence/rank checks -> full
# dynamical model -> latent time -> likelihood-ranked gene heatmap).
# AI usage note: ported/generalized with Claude (Anthropic) from code the
# author originally wrote for one fixed input file; reviewed by the author.

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import scvelo as scv

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True, help="AnnData with spliced/unspliced layers (e.g. from seurat_velocity_export.R)")
parser.add_argument("--groupby", default="seurat_clusters", help="obs column used to color/group every plot [default: %(default)s]")
parser.add_argument("--min_r2", type=float, default=1, help="min_r2 for the steady-state velocity pass [default: %(default)s]")
parser.add_argument("--min_corr_rank_genes", type=float, default=0.3, help="min_corr for rank_velocity_genes [default: %(default)s]")
parser.add_argument("--top_n_genes", type=int, default=15, help="Number of top-likelihood genes shown as scatter panels [default: %(default)s]")
parser.add_argument("--outdir", default="scvelo_analysis")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
scv.settings.figdir = args.outdir
scv.settings.autosave = True

# Pre-processed AnnData object (spliced/unspliced counts + cluster labels
# already attached upstream).
adata = scv.read(args.h5ad)

# Steady-state model: fast, closed-form velocity estimate used as a first pass.
scv.tl.velocity(adata, mode="steady_state", min_r2=args.min_r2)
n_velocity_genes = int(adata.var["velocity_genes"].sum())
print(f"velocity_genes={n_velocity_genes} n_vars={adata.n_vars} n_obs={adata.n_obs}")

# Spliced/unspliced read proportions per cluster -- sanity check that there's
# enough unspliced signal to support velocity estimation before trusting it.
scv.pl.proportions(adata, groupby=args.groupby, save="proportions.svg")

# Steady-state velocity field, projected onto the UMAP embedding.
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color=args.groupby, save="velocity_stream_steady_state.svg")
scv.pl.velocity_embedding(adata, basis="umap", color=args.groupby, arrow_length=3, arrow_size=2, dpi=120, save="velocity_embedding_steady_state.svg")
scv.pl.velocity_graph(adata, threshold=5, color=args.groupby, save="velocity_graph_steady_state.svg")

# Per-cell velocity length/confidence -- flags cells/regions where the
# steady-state model's directionality is unreliable.
scv.tl.velocity_confidence(adata)
keys = "velocity_length", "velocity_confidence"
scv.pl.scatter(adata, c=keys, cmap="coolwarm", perc=[5, 95], save="velocity_confidence.svg")

# Genes whose velocity best explains each cluster's transition, under the
# steady-state model.
scv.tl.rank_velocity_genes(adata, groupby=args.groupby, min_corr=args.min_corr_rank_genes)

# Switch to the full dynamical model: fits per-gene splicing kinetics
# (transcription/splicing/degradation rates) instead of assuming a single
# steady-state ratio -- slower, but more accurate, and required for latent time.
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color=args.groupby, save="velocity_stream_dynamical.svg")

# Latent time: each cell's estimated position along the inferred dynamical
# process, independent of real sampling time -- orders cells by progression
# through the transition rather than by collection timepoint.
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", cmap="coolwarm", size=80, save="latent_time.svg")

# Genes best explained by the dynamical model (highest fit likelihood), shown
# ordered by latent time to reveal expression cascades across the transition.
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color=args.groupby,
                n_convolve=100, cmap="coolwarm", save="latent_time_heatmap.svg")
scv.pl.scatter(adata, basis=top_genes[: args.top_n_genes], ncols=5, color=args.groupby,
                frameon=False, save="top_likelihood_genes.svg")

# Genes whose dynamical-model fit best distinguishes each cluster.
scv.tl.rank_dynamical_genes(adata, groupby=args.groupby)

adata.write(os.path.join(args.outdir, "velocity_processed.h5ad"))

with open(os.path.join(args.outdir, "scvelo_summary_mqc.txt"), "w") as f:
    f.write("# plot_type: 'html'\n# section_name: 'scVelo RNA velocity'\n")
    f.write(f"<p>Cells: <b>{adata.n_obs}</b>, genes: <b>{adata.n_vars}</b></p>\n")
    f.write(f"<p>Steady-state velocity genes: <b>{n_velocity_genes}</b></p>\n")

print("scVelo analysis done.")
