#!/bin/bash
#SBATCH --job-name=repeatsdb-align
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --array=0-16
#SBATCH --output=logs/fetch_alignments_%A_%a.out
#SBATCH --error=logs/fetch_alignments_%A_%a.err

BASEDIR="/groups/banfield/projects/human/nano/scratch/int/kinsey/repeatdetection/solenoid_classifier/repeatsdb-scrape"
ANNOTATIONS_DIR="${BASEDIR}/result-annotations"
OUTPUT_DIR="${BASEDIR}/result-alignments"
SCRIPT="${BASEDIR}/scripts/fetch_alignments.py"

CSV_FILES=(
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_3_1.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_3_2.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_3_3.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_3_4.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_1.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_2.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_3.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_4.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_5.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_6.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_7.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_4_8.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_5_1.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_5_2.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_5_3.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_5_4.csv"
    "${ANNOTATIONS_DIR}/repeatsDB_annotations_5_5.csv"
)

INPUT_CSV="${CSV_FILES[$SLURM_ARRAY_TASK_ID]}"

echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${INPUT_CSV}"
echo "Output directory: ${OUTPUT_DIR}"

python3 "${SCRIPT}" "${INPUT_CSV}" --output-dir "${OUTPUT_DIR}"
