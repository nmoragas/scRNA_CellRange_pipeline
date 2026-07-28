#!/bin/bash
#SBATCH --job-name=fastqc_multiqc
#SBATCH --output=01_preprocessing_scripts/logs/01_pre_fastqc_multiqc_%j.out
#SBATCH --error=01_preprocessing_scripts/logs/01_pre_fastqc_multiqc_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

 
# ---------------------------------------------------------------------------
# 01_run_fastqc_multiqc.sh
# QC de les FASTQ (només R2, lectures cDNA) amb FastQC, agregat amb MultiQC
#
# Ús:
#   sbatch 01_preprocessing/scripts/01_run_fastqc_multiqc.sh
# ---------------------------------------------------------------------------
 
set -euo pipefail


# --- Paths segons l'estructura del projecte -------------------------------
RAW_DATA_DIR="data/raw"
QC_OUTDIR="results/01_preprocessing/1_preqc"
MULTIQC_TITLE="GSE174609 Quality Control"
 
mkdir -p "${QC_OUTDIR}"
mkdir -p 01_preprocessing_scripts/logs
 
# --- Càrrega de mòduls (ajusta noms/versions al teu HPC) -------------------
module load apps/multiqc/1.25.1
module load apps/fastqc/0.12.1 
 
# --- FastQC sobre R2 (lectures cDNA) ----------------------------------------
echo "[$(date)] Corrent FastQC..."
fastqc "${RAW_DATA_DIR}"/*_R2_001.fastq.gz \
    --outdir "${QC_OUTDIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --quiet
 
# --- Agregació amb MultiQC --------------------------------------------------
echo "[$(date)] Corrent MultiQC..."
multiqc "${QC_OUTDIR}" \
    --outdir "${QC_OUTDIR}" \
    --filename multiqc_report \
    --title "${MULTIQC_TITLE}"
 
echo "[$(date)] Fet. Resultats a: ${QC_OUTDIR}"