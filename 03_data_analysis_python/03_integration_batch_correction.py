#!/usr/bin/env python3
# ==============================================================================
# 03_integration_batch_correction.py
# Corregeix l'efecte de lot entre mostres/experiments amb Harmony, i genera
# una comparativa visual abans/despres per decidir si cal aplicar-ho.
#
# Quan es necessari: si al UMAP de 00_pca_clustering.py les cellules se separen
# principalment per "sample" en comptes de per tipus cel·lular, hi ha efecte
# de lot i cal integrar abans de interpretar els clusters biologicament.
#
# Input : results/python/qc_filtering/adata_normalized_hvg_scaled.h5ad
# Output: results/python/data_analysis/adata_integrated.h5ad
#         results/python/data_analysis/plots/umap_before_after_integration.png
#
# Us: python 3_data_analysis/03_integration_batch_correction.py
# Requereix: scanpy, harmonypy
# ==============================================================================

import os
import scanpy as sc
import scanpy.external as sce

PROJECT_DIR = os.getcwd()
QC_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")
OUT_DIR = os.path.join(PROJECT_DIR, "results", "python", "data_analysis")
PLOTS_DIR = os.path.join(OUT_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

adata = sc.read_h5ad(os.path.join(QC_DIR, "adata_normalized_hvg_scaled.h5ad"))

# ---- PCA sense integrar (referencia "abans") ----
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)
adata.obsm["X_umap_before"] = adata.obsm["X_umap"].copy()

# ---- Integracio amb Harmony ----
sce.pp.harmony_integrate(adata, key="sample")
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.obsm["X_umap_after"] = adata.obsm["X_umap"].copy()

# ---- Clustering sobre l'espai integrat ----
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_integrated")

# ---- Comparativa visual abans/despres ----
adata.obsm["X_umap"] = adata.obsm["X_umap_before"]
sc.pl.umap(adata, color="sample", show=False, save="_before_integration.png")

adata.obsm["X_umap"] = adata.obsm["X_umap_after"]
sc.pl.umap(adata, color="sample", show=False, save="_after_integration.png")
sc.pl.umap(adata, color="leiden_integrated", show=False, save="_after_integration_clusters.png")

os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

out_h5ad = os.path.join(OUT_DIR, "adata_integrated.h5ad")
adata.write(out_h5ad)

print(f"Fet. Comparativa abans/despres a: {PLOTS_DIR}")
print(f"AnnData integrat guardat a: {out_h5ad}")
print("\nSi el UMAP 'despres' mostra mostres ben barrejades (i no per sample),")
print("la integracio ha funcionat. Continuar l'analisi (01, 02) amb aquest objecte")
print("en comptes de 'adata_clustered.h5ad' si l'efecte de lot era rellevant.")
