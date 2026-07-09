#!/bin/bash --login
#SBATCH --job-name=fastqc_multiqc
#SBATCH --output=fastqc_multiqc_%j.out
#SBATCH --error=fastqc_multiqc_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

set -euo pipefail

###############################################################################
# EDIT THIS BLOCK FOR YOUR DATA AND CLUSTER
###############################################################################

FASTQ_DIR="FASTQ"
OUTPUT_DIR="fastqc_multiqc_output"

FASTQC_MODULE="FastQC/0.12.1-Java-11"
MULTIQC_MODULE="MultiQC/1.14-foss-2022b"
PARALLEL_MODULE="parallel/20230722-GCCcore-12.2.0"

# Run this many FastQC processes at once. Each process handles one FASTQ file.
PARALLEL_JOBS=16
FASTQC_THREADS=1

###############################################################################
# SETUP
###############################################################################

module purge
module load "$FASTQC_MODULE"
module load "$MULTIQC_MODULE"
module load "$PARALLEL_MODULE"

command -v fastqc >/dev/null 2>&1 || { echo "ERROR: fastqc not found after loading modules." >&2; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc not found after loading modules." >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: parallel not found after loading modules." >&2; exit 1; }

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

FASTQC_OUTDIR="$OUTPUT_DIR/fastqc_reports"
MULTIQC_OUTDIR="$OUTPUT_DIR/multiqc_report"
FASTQC_JOBLOG="$OUTPUT_DIR/fastqc_parallel_joblog.tsv"

mkdir -p "$FASTQC_OUTDIR" "$MULTIQC_OUTDIR"

###############################################################################
# FIND FASTQ FILES
###############################################################################

mapfile -t FASTQ_FILES < <(
    find "$FASTQ_DIR" -maxdepth 1 -type f \( \
        -name "*.fq.gz" -o \
        -name "*.fastq.gz" \
    \) | sort
)

if [[ "${#FASTQ_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No .fq.gz or .fastq.gz files found in: $FASTQ_DIR" >&2
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ files."
echo "Writing FastQC reports to: $FASTQC_OUTDIR"

###############################################################################
# RUN FASTQC, THEN MULTIQC
###############################################################################

parallel --jobs "$PARALLEL_JOBS" \
         --joblog "$FASTQC_JOBLOG" \
         --halt soon,fail=1 \
         fastqc --quiet --threads "$FASTQC_THREADS" --outdir "$FASTQC_OUTDIR" ::: "${FASTQ_FILES[@]}"

multiqc "$FASTQC_OUTDIR" \
    --outdir "$MULTIQC_OUTDIR" \
    --filename "multiqc_FASTQ_report.html" \
    --force

if [[ ! -s "$MULTIQC_OUTDIR/multiqc_FASTQ_report.html" ]]; then
    echo "ERROR: MultiQC report was not created." >&2
    exit 1
fi

echo "Done."
echo "FastQC reports:  $FASTQC_OUTDIR"
echo "MultiQC report:  $MULTIQC_OUTDIR/multiqc_FASTQ_report.html"
