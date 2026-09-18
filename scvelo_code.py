# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
import scvelo as scv
adata = scv.read("lgl.h5ad")
scv.tl.velocity(adata, mode='steady_state', min_r2=1)
print(adata.var['velocity_genes'].sum(), adata.n_vars, adata.n_obs)

scv.pl.proportions(adata, groupby='seurat_clusters')

scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
scv.pl.velocity_graph(adata, threshold=5, color='seurat_clusters')
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters')
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='coolwarm', size=80)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, cmap='coolwarm')
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, color='seurat_clusters', frameon=False)
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
