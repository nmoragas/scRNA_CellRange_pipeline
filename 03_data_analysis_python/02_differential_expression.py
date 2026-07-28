#!/usr/bin/env python3
# ==============================================================================
# 02_differential_expression.py
# Expressio diferencial entre dues condicions (p.ex. tractat vs control),
# tant a nivell global com per cluster/tipus cel·lular.
#
# IMPORTANT: aquest script assumeix que adata.obs conte una columna "condition"
# (p.ex. "control"/"treated"). Cal afegir-la manualment (o derivar-la del nom
# de mostra) si encara no existeix.
#
# Input : results/python/data_analysis/adata_annotated.h5ad
# Output: results/python/data_analysis/de_global.csv
#         results/python/data_analysis/de_per_cluster/<cluster>.csv
#
# Us: python 3_data_analysis/02_differential_expression.py
# ==============================================================================

import os
import scanpy as sc

PROJECT_DIR = os.getcwd()
IN_DIR = os.path.join(PROJECT_DIR, "results", "python", "data_analysis")
DE_DIR = os.path.join(IN_DIR, "de_per_cluster")
os.makedirs(DE_DIR, exist_ok=True)

adata = sc.read_h5ad(os.path.join(IN_DIR, "adata_annotated.h5ad"))

if "condition" not in adata.obs.columns:
    print("[AVIS] No existeix adata.obs['condition']. Exemple per crear-la a partir")
    print("       del nom de mostra:")
    print("       adata.obs['condition'] = adata.obs['sample'].map({'sample1':'control', 'sample2':'treated'})")
    raise SystemExit("Afegeix la columna 'condition' abans de continuar.")

conditions = adata.obs["condition"].unique().tolist()
if len(conditions) != 2:
    raise SystemExit(f"Aquest script espera exactament 2 condicions, trobades: {conditions}")

groupA, groupB = conditions

# ---- DE global (totes les cellules juntes, sense distingir cluster/tipus) ----
# NOTA: la DE "pseudobulk" (agregant comptatges per mostra abans de fer el test,
# p.ex. amb pydeseq2) es mes robusta estadisticament que comparar cellula a
# cellula, ja que evita pseudo-repliques. Aquest bloc fa servir Wilcoxon per
# simplicitat; per a resultats publicables, considerar pseudobulk + DESeq2.
sc.tl.rank_genes_groups(
    adata, groupby="condition", groups=[groupB], reference=groupA, method="wilcoxon"
)
de_global = sc.get.rank_genes_groups_df(adata, group=groupB)
de_global.to_csv(os.path.join(IN_DIR, "de_global.csv"), index=False)
print(f"DE global ({groupB} vs {groupA}) guardada: {len(de_global)} gens testats")

# ---- DE per cluster/tipus cel·lular ----
cluster_col = "cell_type" if "cell_type" in adata.obs.columns else "leiden"

for cl in sorted(adata.obs[cluster_col].unique(), key=str):
    sub = adata[adata.obs[cluster_col] == cl].copy()
    if sub.obs["condition"].nunique() < 2:
        print(f"  [AVIS] Cluster {cl}: nomes una condicio present, saltant")
        continue

    sc.tl.rank_genes_groups(
        sub, groupby="condition", groups=[groupB], reference=groupA, method="wilcoxon"
    )
    de_cl = sc.get.rank_genes_groups_df(sub, group=groupB)
    de_cl.to_csv(os.path.join(DE_DIR, f"{cl}.csv"), index=False)
    print(f"  Cluster {cl}: {len(de_cl)} gens testats")

print(f"\nFet. Resultats a: {IN_DIR}/de_global.csv i {DE_DIR}/")
