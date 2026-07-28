#!/usr/bin/env python3
# ==============================================================================
# 02_filtering.py
# Aplica els llindars de filtratge decidits despres d'inspeccionar els plots
# de 00_qc_import.py i 01_doublet_detection.py, i elimina els doublets predits.
#
# IMPORTANT: els llindars de sota son valors d'exemple habituals per teixit
# huma/raton estandard. S'HAN D'AJUSTAR segons la distribucio observada als
# violin plots de cada projecte concret -- no hi ha llindars universals.
#
# Input : results/python/qc_filtering/adata_doublets.h5ad
# Output: results/python/qc_filtering/adata_filtered.h5ad
#
# Us: python 2_qc_filtering/02_filtering.py
# ==============================================================================

import os
import scanpy as sc

PROJECT_DIR = os.getcwd()
IN_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")

adata = sc.read_h5ad(os.path.join(IN_DIR, "adata_doublets.h5ad"))
n_before = adata.n_obs

# ---- Llindars a ajustar segons el projecte ----
MIN_GENES = 200
MAX_GENES = 6000       # cellules amb massa gens sovint son doublets no detectats
MIN_COUNTS = 500
MAX_PCT_MT = 15         # 15% per teixit estandard; pujar per teixits amb metabolisme alt

# ---- Filtratge de gens (eliminar gens no expressats en cap cellula) ----
sc.pp.filter_genes(adata, min_cells=3)

# ---- Filtratge de cellules per metriques de QC ----
mask = (
    (adata.obs["n_genes_by_counts"] >= MIN_GENES) &
    (adata.obs["n_genes_by_counts"] <= MAX_GENES) &
    (adata.obs["total_counts"] >= MIN_COUNTS) &
    (adata.obs["pct_counts_mt"] <= MAX_PCT_MT) &
    (~adata.obs["predicted_doublet"])
)
adata = adata[mask].copy()

n_after = adata.n_obs
print(f"Cellules abans del filtratge : {n_before}")
print(f"Cellules despres del filtratge: {n_after} ({100*n_after/n_before:.1f}%)")
print(f"Cellules eliminades           : {n_before - n_after}")

print("\nResum per mostra despres del filtratge:")
print(adata.obs.groupby("sample", observed=True).size())

out_h5ad = os.path.join(IN_DIR, "adata_filtered.h5ad")
adata.write(out_h5ad)

print(f"\nFet. AnnData filtrat guardat a: {out_h5ad}")
