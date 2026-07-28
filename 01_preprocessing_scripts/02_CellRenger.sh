#!/bin/bash
#SBATCH --job-name=cellranger_array
#SBATCH --output=01_preprocessing_scripts/logs/02_cellranger_%j.out
#SBATCH --error=01_preprocessing_scripts/logs/02_cellranger_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --array=1-1    

 
#-----------------------------------------------
# Cell Ranger Count - Array Job for Multiple Samples
#-----------------------------------------------
 
# Cell Ranger al PATH
export PATH=/mnt/hydra/ubs/shared/users/Nuria/courses/scRNA/CellRenger/data/cellranger-10.1.0:$PATH 
 
# Definim rutes del projecte
PROJECT_DIR=$(pwd)
FASTQ_DIR=${PROJECT_DIR}/data/raw
OUTPUT_DIR=${PROJECT_DIR}/results/01_preprocessing/02_cellranger
REFERENCE=${PROJECT_DIR}/data/refdata-gex-GRCh38-2024-A

# Creem directoris necessaris
mkdir -p ${OUTPUT_DIR}

# ---- Detecció automàtica de mostres a partir dels FASTQ ----
# Cell Ranger espera noms tipus: SAMPLE_S1_L001_R1_001.fastq.gz
# Extraiem el nom de mostra (tot el que hi ha abans de "_S<num>_")
mapfile -t SAMPLES < <(ls ${FASTQ_DIR}/*.fastq.gz | \
    sed -E 's/.*\///; s/_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz//' | \
    sort -u)

# Mapegem SLURM_ARRAY_TASK_ID (1-indexat) a la mostra corresponent (0-indexat)
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Paràmetres de Cell Ranger
EXPECTED_CELLS=5000
LOCALCORES=16
LOCALMEM=110

# Anem al directori de sortida (Cell Ranger crea subdirectoris aquí)
cd ${OUTPUT_DIR}

# Executem Cell Ranger Count
cellranger count \
    --id=${SAMPLE_ID} \
    --fastqs=${FASTQ_DIR} \
    --sample=${SAMPLE_ID} \
    --transcriptome=${REFERENCE} \
    --expect-cells=${EXPECTED_CELLS} \
    --localcores=${LOCALCORES} \
    --localmem=${LOCALMEM} \
    --chemistry=auto \
    --create-bam=false