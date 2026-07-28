# scRNA-seq Pipeline: FASTQ → Count Matrix → QC → Analysis

Preprocessing and analysis pipeline for scRNA-seq data (10x Genomics), from
raw fastq files through to clustering, cell type annotation, and differential
expression.

## Directory structure

```
.
├── 01_preprocessing_scripts/     # Raw fastq QC + quantification with CellRanger
│   ├── 01_pre_QC.sh                 FastQC + MultiQC on raw fastq files
│   ├── 02_CellRanger.sh              cellranger count (alignment + matrix)
│   └── logs/                        SLURM logs (not version-controlled)
│
├── 02_qc_filtering_python/       # Technical cleanup (Python/Scanpy)
│   ├── 00_qc_import.py              import matrix → AnnData + QC metrics
│   ├── 01_doublet_detection.py      Scrublet
│   ├── 02_filtering.py              apply thresholds + remove doublets
│   └── 03_normalization.py          normalization + HVG selection + scaling
│
├── 03_data_analysis_python/      # Biological analysis (Python/Scanpy)
│   ├── 00_pca_clustering.py         PCA + (optional Harmony) + Leiden + UMAP
│   ├── 01_marker_genes_annotation.py  marker genes + cell type annotation
│   ├── 02_differential_expression.py  DE between conditions
│   └── 03_integration_batch_correction.py  Harmony, before/after comparison
│
├── data/
│   ├── raw/                         input fastq.gz files (not version-controlled)
│   ├── cellranger-10.1.0/           CellRanger installation (not version-controlled)
│   └── refdata-gex-GRCh38-2024-A/   10x reference package (not version-controlled)
│
└── results/
    └── 01_preprocessing/
        ├── 1_preqc/                  FastQC + MultiQC output
        └── 02_cellranger/             cellranger count output per sample
```

## Execution order

```bash
# 1. QC on raw fastq files
sbatch 01_preprocessing_scripts/01_pre_QC.sh

# 2. Quantification with CellRanger
sbatch 01_preprocessing_scripts/02_CellRanger.sh

# 3. Per-cell QC and filtering
python 02_qc_filtering_python/00_qc_import.py
python 02_qc_filtering_python/01_doublet_detection.py
python 02_qc_filtering_python/02_filtering.py
python 02_qc_filtering_python/03_normalization.py

# 4. Analysis
python 03_data_analysis_python/00_pca_clustering.py
python 03_data_analysis_python/01_marker_genes_annotation.py
python 03_data_analysis_python/02_differential_expression.py   # requires obs["condition"]
python 03_data_analysis_python/03_integration_batch_correction.py  # only if needed
```

## Requirements

- CellRanger 10.1.0 (installed at `data/cellranger-10.1.0/`)
- 10x reference: `refdata-gex-GRCh38-2024-A`
- Python: `scanpy>=1.9`, `anndata`, `harmonypy`, `leidenalg`, `pandas`

## Points to adjust per project

- Mitochondrial gene pattern in `00_qc_import.py` ("MT-" human vs "mt-" mouse)
- Filtering thresholds in `02_filtering.py` (no universal values — set based on
  the QC distributions observed in your own data)
- Known marker genes in `01_marker_genes_annotation.py`
- `condition` column required before running `02_differential_expression.py`

## License / authorship

Developed by Núria Moragas
