# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
# A generalized, parameterized Nextflow version of this script (`--mode
# scrna_velocity`) is available in nextflow_pipeline/ -- see its README, or
# the "Reproducing this analysis via Nextflow" section of the repo-root
# README.
import scvelo as scv

# Pre-processed AnnData object (spliced/unspliced counts + Seurat cluster
# labels already attached from the upstream Seurat/loom pipeline).
adata = scv.read("lgl.h5ad")

# Steady-state model: fast, closed-form velocity estimate used as a first pass.
scv.tl.velocity(adata, mode='steady_state', min_r2=1)
print(adata.var['velocity_genes'].sum(), adata.n_vars, adata.n_obs)

# Spliced/unspliced read proportions per cluster -- sanity check that there's
# enough unspliced signal to support velocity estimation before trusting it.
scv.pl.proportions(adata, groupby='seurat_clusters')

# Steady-state velocity field, projected onto the UMAP embedding.
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
scv.pl.velocity_graph(adata, threshold=5, color='seurat_clusters')

# Per-cell velocity length/confidence -- flags cells/regions where the
# steady-state model's directionality is unreliable.
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

# Genes whose velocity best explains each cluster's transition, under the
# steady-state model.
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)

# Switch to the full dynamical model: fits per-gene splicing kinetics
# (transcription/splicing/degradation rates) instead of assuming a single
# steady-state ratio -- slower, but more accurate, and required for latent time.
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters')

# Latent time: each cell's estimated position along the inferred dynamical
# process, independent of real sampling time -- orders cells by progression
# through the transition rather than by collection timepoint.
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='coolwarm', size=80)

# Genes best explained by the dynamical model (highest fit likelihood),
# shown ordered by latent time to reveal expression cascades across the
# transition.
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, cmap='coolwarm')
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, color='seurat_clusters', frameon=False)

# Genes whose dynamical-model fit best distinguishes each cluster.
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
