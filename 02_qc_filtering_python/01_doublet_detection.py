#!/usr/bin/env python3
# ==============================================================================
# 01_doublet_detection.py
# Deteccio de doublets amb Scrublet (integrat nativament a Scanpy des de la
# v1.9+), calculat per separat per mostra (batch_key="sample").
#
# Input : results/python/qc_filtering/adata_raw.h5ad
# Output: results/python/qc_filtering/adata_doublets.h5ad
#         results/python/qc_filtering/plots/doublet_scores_*.png
#
# Us: python 2_qc_filtering/01_doublet_detection.py
# Requereix: scanpy>=1.9, scikit-image (dependencia de scrublet)
# ==============================================================================

import os
import scanpy as sc

PROJECT_DIR = os.getcwd()
IN_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")
PLOTS_DIR = os.path.join(IN_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

adata = sc.read_h5ad(os.path.join(IN_DIR, "adata_raw.h5ad"))

# ---- Scrublet per mostra (expected_doublet_rate tipic per 10x: 0.05-0.08) ----
sc.pp.scrublet(adata, batch_key="sample", expected_doublet_rate=0.06)

# adata.obs ara conte:
#   doublet_score       -> score continu (0-1)
#   predicted_doublet    -> booleà

n_doublets = int(adata.obs["predicted_doublet"].sum())
print(f"Doublets predits: {n_doublets} / {adata.n_obs} cellules "
      f"({100 * n_doublets / adata.n_obs:.2f}%)")

# ---- Resum per mostra ----
summary = adata.obs.groupby("sample", observed=True)["predicted_doublet"].agg(["sum", "count"])
summary["pct"] = 100 * summary["sum"] / summary["count"]
print(summary)
summary.to_csv(os.path.join(IN_DIR, "doublet_summary_per_sample.csv"))

# ---- Histograma dels scores ----
sc.pl.scrublet_score_distribution(adata, show=False, save="_doublet_scores.png")
os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

out_h5ad = os.path.join(IN_DIR, "adata_doublets.h5ad")
adata.write(out_h5ad)

print(f"\nFet. AnnData amb scores de doublets guardat a: {out_h5ad}")
print("Els doublets s'eliminaran al pas de filtratge (02_filtering.py).")
