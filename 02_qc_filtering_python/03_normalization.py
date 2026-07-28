#!/usr/bin/env python3
# ==============================================================================
# 03_normalization.py
# Normalitza la profunditat de sequenciacio, aplica log1p, selecciona gens
# altament variables (HVG) i escala les dades. Deixa l'AnnData llest per a
# PCA/clustering al modul 3_data_analysis.
#
# Input : results/python/qc_filtering/adata_filtered.h5ad
# Output: results/python/qc_filtering/adata_normalized.h5ad
#
# Us: python 2_qc_filtering/03_normalization.py
# ==============================================================================

import os
import scanpy as sc

PROJECT_DIR = os.getcwd()
IN_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")

adata = sc.read_h5ad(os.path.join(IN_DIR, "adata_filtered.h5ad"))

# ---- Guardar comptatges crus (utils per a analisis posteriors, p.ex. pseudobulk) ----
adata.layers["counts"] = adata.X.copy()

# ---- Normalitzacio de profunditat + log-transformacio ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- Gens altament variables (HVG) ----
# Es calculen per mostra (batch_key) i despres es combinen, per evitar que
# diferencies tecniques entre mostres dominin la seleccio de HVGs.
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, batch_key="sample", flavor="seurat"
)

adata.raw = adata   # guarda l'expressio log-normalitzada completa (tots els gens)

# ---- Escalar nomes els HVG per a PCA (evita que gens d'expressio alta dominin) ----
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)

# Es guarden dos objectes: el complet (per consultar qualsevol gen despres)
# i el reduit a HVG+escalat (input directe per PCA)
out_full = os.path.join(IN_DIR, "adata_normalized.h5ad")
out_hvg = os.path.join(IN_DIR, "adata_normalized_hvg_scaled.h5ad")
adata.write(out_full)
adata_hvg.write(out_hvg)

print(f"Fet.")
print(f"  Objecte complet (log-normalitzat, tots els gens): {out_full}")
print(f"  Objecte HVG escalat (input per PCA)              : {out_hvg}")
print(f"  Nombre de HVG seleccionats: {adata.var['highly_variable'].sum()}")
