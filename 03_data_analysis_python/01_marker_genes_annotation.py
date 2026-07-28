#!/usr/bin/env python3
# ==============================================================================
# 01_marker_genes_annotation.py
# Identifica gens marcadors per cada cluster (rank_genes_groups) i genera
# les visualitzacions habituals per a l'anotacio manual de tipus cel·lulars.
# Inclou tambe un exemple d'anotacio automatica per score de marcadors coneguts.
#
# Input : results/python/data_analysis/adata_clustered.h5ad
# Output: results/python/data_analysis/marker_genes_per_cluster.csv
#         results/python/data_analysis/adata_annotated.h5ad
#         results/python/data_analysis/plots/dotplot_markers.png, umap_celltype.png
#
# Us: python 3_data_analysis/01_marker_genes_annotation.py
# ==============================================================================

import os
import pandas as pd
import scanpy as sc

PROJECT_DIR = os.getcwd()
IN_DIR = os.path.join(PROJECT_DIR, "results", "python", "data_analysis")
PLOTS_DIR = os.path.join(IN_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

adata = sc.read_h5ad(os.path.join(IN_DIR, "adata_clustered.h5ad"))

# ---- Gens marcadors per cluster (test de Wilcoxon, estandard en scRNA-seq) ----
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# ---- Exportar taula de marcadors (top 25 per cluster) ----
markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_df = markers_df.groupby("group").head(25)
markers_df.to_csv(os.path.join(IN_DIR, "marker_genes_per_cluster.csv"), index=False)

sc.pl.rank_genes_groups_dotplot(
    adata, n_genes=5, show=False, save="_top_markers.png"
)
os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

# ==============================================================================
# ANOTACIO: exemple amb marcadors coneguts (AJUSTAR AL TEIXIT/ORGANISME CONCRET)
# Substituir aquest diccionari pels marcadors rellevants del teu experiment.
# ==============================================================================
marker_genes = {
    "T cells":        ["CD3D", "CD3E", "CD3G"],
    "B cells":        ["MS4A1", "CD79A", "CD79B"],
    "NK cells":       ["NKG7", "GNLY", "KLRD1"],
    "Monocytes":      ["CD14", "LYZ", "FCGR3A"],
    "Dendritic cells": ["FCER1A", "CST3"],
    "Platelets":      ["PPBP"],
}

# Filtrar nomes els marcadors presents a l'objecte
marker_genes = {
    ct: [g for g in genes if g in adata.var_names]
    for ct, genes in marker_genes.items()
}
marker_genes = {ct: genes for ct, genes in marker_genes.items() if len(genes) > 0}

if marker_genes:
    sc.tl.score_genes(adata, gene_list=[g for gl in marker_genes.values() for g in gl],
                       score_name="_tmp_all_markers")
    for ct, genes in marker_genes.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{ct}")

    sc.pl.dotplot(
        adata, marker_genes, groupby="leiden", show=False, save="_known_markers.png"
    )
    os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

    print("Scores mitjans per cluster (usar per assignar manualment el tipus cel·lular):")
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    print(adata.obs.groupby("leiden", observed=True)[score_cols].mean())

# ---- Placeholder: assignacio manual final (editar segons els scores anteriors) ----
# cluster_to_celltype = {"0": "T cells", "1": "Monocytes", ...}
# adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype)

out_h5ad = os.path.join(IN_DIR, "adata_annotated.h5ad")
adata.write(out_h5ad)

print(f"\nFet. Taula de marcadors: {os.path.join(IN_DIR, 'marker_genes_per_cluster.csv')}")
print(f"AnnData amb scores/marcadors guardat a: {out_h5ad}")
print("Pas manual següent: assignar 'cell_type' a cada cluster segons els marcadors/scores.")
