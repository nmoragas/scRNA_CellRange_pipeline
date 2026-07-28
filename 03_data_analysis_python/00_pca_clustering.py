#!/usr/bin/env python3
# ==============================================================================
# 00_pca_clustering.py
# PCA, correccio de lot (Harmony) opcional, calcul del graf de veins,
# clustering (Leiden) i projeccio UMAP per visualitzacio.
#
# Input : results/python/qc_filtering/adata_normalized_hvg_scaled.h5ad
#         results/python/qc_filtering/adata_normalized.h5ad (per recuperar tots els gens)
# Output: results/python/data_analysis/adata_clustered.h5ad
#         results/python/data_analysis/plots/umap_clusters.png, umap_sample.png
#
# Us: python 3_data_analysis/00_pca_clustering.py
# Requereix: scanpy, harmonypy (nomes si INTEGRATE_BATCH=True), leidenalg
# ==============================================================================

import os
import scanpy as sc

PROJECT_DIR = os.getcwd()
QC_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "python", "data_analysis")
PLOTS_DIR = os.path.join(OUT_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# ---- Activar si hi ha efecte de lot fort entre mostres/condicions ----
INTEGRATE_BATCH = False   # posar True per aplicar Harmony (veure 03_integration_batch_correction.py)

adata_hvg = sc.read_h5ad(os.path.join(QC_DIR, "adata_normalized_hvg_scaled.h5ad"))

# ---- PCA ----
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

rep_to_use = "X_pca"
if INTEGRATE_BATCH:
    import scanpy.external as sce
    sce.pp.harmony_integrate(adata_hvg, key="sample")
    rep_to_use = "X_pca_harmony"

# ---- Graf de veins + clustering ----
sc.pp.neighbors(adata_hvg, n_neighbors=15, use_rep=rep_to_use)
sc.tl.leiden(adata_hvg, resolution=1.0, key_added="leiden")
sc.tl.umap(adata_hvg)

print("Nombre de clusters detectats:", adata_hvg.obs["leiden"].nunique())
print(adata_hvg.obs["leiden"].value_counts())

# ---- Plots ----
sc.pl.umap(adata_hvg, color="leiden", show=False, save="_clusters.png")
sc.pl.umap(adata_hvg, color="sample", show=False, save="_sample.png")
os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

# ---- Recuperar tots els gens (no nomes HVG) i afegir-hi els resultats de clustering ----
adata_full = sc.read_h5ad(os.path.join(QC_DIR, "adata_normalized.h5ad"))
adata_full.obs["leiden"] = adata_hvg.obs["leiden"].values
adata_full.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata_full.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
if INTEGRATE_BATCH:
    adata_full.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]

out_h5ad = os.path.join(OUT_DIR, "adata_clustered.h5ad")
adata_full.write(out_h5ad)

print(f"\nFet. AnnData amb clustering guardat a: {out_h5ad}")
