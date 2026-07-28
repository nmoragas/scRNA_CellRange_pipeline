#!/usr/bin/env python3
# ==============================================================================
# 00_qc_import.py
# Carrega les matrius filtrades de STARsolo (una per mostra), les combina en
# un unic AnnData i calcula metriques de QC per cellula.
#
# Input : results/02_starsolo/<SAMPLE_ID>/Solo.out/Gene/filtered/
# Output: results/python/qc_filtering/adata_raw.h5ad
#         results/python/qc_filtering/plots/qc_violin.png, qc_scatter.png
#
# Us: python 2_qc_filtering/00_qc_import.py
# Requereix: scanpy, anndata, matplotlib
# ==============================================================================

import os
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1

PROJECT_DIR = os.getcwd()
STARSOLO_DIR = os.path.join(PROJECT_DIR, "results", "02_starsolo")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "results", "python", "qc_filtering")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

samples = sorted(
    d for d in os.listdir(STARSOLO_DIR)
    if os.path.isdir(os.path.join(STARSOLO_DIR, d))
)
print(f"Mostres detectades: {samples}")

adatas = {}
for sample in samples:
    mat_dir = os.path.join(STARSOLO_DIR, sample, "Solo.out", "Gene", "filtered")
    if not os.path.isdir(mat_dir):
        print(f"  [AVIS] No trobat {mat_dir}, saltant {sample}")
        continue

    a = sc.read_mtx(os.path.join(mat_dir, "matrix.mtx")).T
    barcodes = [l.strip() for l in open(os.path.join(mat_dir, "barcodes.tsv"))]
    features = [l.strip().split("\t") for l in open(os.path.join(mat_dir, "features.tsv"))]

    a.obs_names = barcodes
    a.var_names = [f[1] for f in features]   # nom del gen (columna 2 de features.tsv)
    a.var["gene_id"] = [f[0] for f in features]
    a.var_names_make_unique()
    a.obs["sample"] = sample

    print(f"  {sample}: {a.n_obs} cellules, {a.n_vars} gens")
    adatas[sample] = a

# ---- Combinar totes les mostres ----
adata = ad.concat(adatas, label="sample_batch", join="outer", index_unique="_")
adata.obs["sample"] = adata.obs["sample"].astype("category")

# ---- Metriques de QC ----
# Ajustar el patró mitocondrial segons l'organisme: "MT-" (huma), "mt-" (ratoli)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
)

# ---- Plots de QC ----
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample",
    jitter=0.4,
    multi_panel=True,
    show=False,
    save="_qc_violin.png",
)
sc.pl.scatter(
    adata, x="total_counts", y="pct_counts_mt", color="sample",
    show=False, save="_qc_scatter_mt.png",
)
sc.pl.scatter(
    adata, x="total_counts", y="n_genes_by_counts", color="sample",
    show=False, save="_qc_scatter_genes.png",
)
# scanpy guarda les figures a ./figures/ per defecte; les movem a la carpeta de resultats
os.system(f"mv ./figures/*.png {PLOTS_DIR}/ 2>/dev/null || true")

# ---- Guardar AnnData sense filtrar ----
out_h5ad = os.path.join(OUTPUT_DIR, "adata_raw.h5ad")
adata.write(out_h5ad)

print(f"\nFet. AnnData combinat (sense filtrar) guardat a: {out_h5ad}")
print(f"Plots de QC a: {PLOTS_DIR}")
print("Revisa els plots abans d'aplicar llindars al script 02_filtering.py")
