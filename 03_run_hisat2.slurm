#!/bin/bash --login
#SBATCH --job-name=hisat2_cds_bam
#SBATCH --output=hisat2_cds_bam_%j.out
#SBATCH --error=hisat2_cds_bam_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=120G

set -euo pipefail

###############################################################################
# EDIT THIS BLOCK FOR YOUR DATA AND CLUSTER
###############################################################################

FASTQ_DIR="FASTQ"
REFERENCE_FASTA="Tidestromia_oblongifolia_DP_CDS.longest_isoform_per_locus.fasta"
OUTPUT_DIR="hisat2_cds_mapping"

HISAT2_MODULE="HISAT2/2.2.1-gompi-2023a"
SAMTOOLS_MODULE="SAMtools/1.18-GCC-12.3.0"
PARALLEL_MODULE="parallel/20230722-GCCcore-12.3.0"

# Run 10 paired-end samples at a time, using 10 HISAT2 threads per sample.
PARALLEL_JOBS=10
HISAT2_THREADS=10
SAMTOOLS_THREADS=1
INDEX_THREADS=10

# For a CDS/transcript reference, spliced alignment is usually not needed.
HISAT2_EXTRA_ARGS="--no-spliced-alignment"

###############################################################################
# SETUP
###############################################################################

module purge
module load "$HISAT2_MODULE"
module load "$SAMTOOLS_MODULE"
module load "$PARALLEL_MODULE"

command -v hisat2 >/dev/null 2>&1 || { echo "ERROR: hisat2 not found after loading modules." >&2; exit 1; }
command -v hisat2-build >/dev/null 2>&1 || { echo "ERROR: hisat2-build not found after loading modules." >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found after loading modules." >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: parallel not found after loading modules." >&2; exit 1; }

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

if [[ ! -s "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found or empty: $REFERENCE_FASTA" >&2
    exit 1
fi

INDEX_DIR="$OUTPUT_DIR/hisat2_index"
BAM_DIR="$OUTPUT_DIR/bam"
LOG_DIR="$OUTPUT_DIR/logs"
PAIR_TSV="$OUTPUT_DIR/sample_pairs.tsv"
MAP_JOBLOG="$OUTPUT_DIR/hisat2_parallel_joblog.tsv"
INDEX_PREFIX="$INDEX_DIR/reference_cds"

mkdir -p "$INDEX_DIR" "$BAM_DIR" "$LOG_DIR"

###############################################################################
# BUILD HISAT2 INDEX
###############################################################################

echo "Building HISAT2 index from: $REFERENCE_FASTA"
hisat2-build -p "$INDEX_THREADS" "$REFERENCE_FASTA" "$INDEX_PREFIX"

###############################################################################
# FIND PAIRED FASTQ FILES
###############################################################################

mapfile -t R1_FILES < <(
    find "$FASTQ_DIR" -maxdepth 1 -type f \( \
        -name "*_R1.fq.gz" -o \
        -name "*_R1.fastq.gz" \
    \) | sort
)

if [[ "${#R1_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No R1 FASTQ files found in: $FASTQ_DIR" >&2
    exit 1
fi

: > "$PAIR_TSV"

for R1 in "${R1_FILES[@]}"; do
    if [[ "$R1" == *_R1.fq.gz ]]; then
        SAMPLE=$(basename "${R1%_R1.fq.gz}")
        R2="${R1%_R1.fq.gz}_R2.fq.gz"
    elif [[ "$R1" == *_R1.fastq.gz ]]; then
        SAMPLE=$(basename "${R1%_R1.fastq.gz}")
        R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    else
        echo "ERROR: Unexpected R1 filename: $R1" >&2
        exit 1
    fi

    if [[ ! -s "$R2" ]]; then
        echo "ERROR: Missing R2 file for sample $SAMPLE: $R2" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE" \
        "$R1" \
        "$R2" \
        "$BAM_DIR/${SAMPLE}.bam" \
        "$LOG_DIR/${SAMPLE}.hisat2.log" >> "$PAIR_TSV"
done

echo "Found ${#R1_FILES[@]} paired-end samples."
echo "Writing BAM files to: $BAM_DIR"

###############################################################################
# MAP READS AND WRITE SORTED BAM FILES
###############################################################################

# HISAT2 writes SAM to stdout. SAMtools converts that stream to a sorted BAM.
# --no-unal and samtools -F 4 exclude unmapped reads from the final BAM file.
export SHELL=/bin/bash
export INDEX_PREFIX HISAT2_THREADS SAMTOOLS_THREADS HISAT2_EXTRA_ARGS

parallel --colsep '\t' \
         --jobs "$PARALLEL_JOBS" \
         --joblog "$MAP_JOBLOG" \
         --halt soon,fail=1 \
         'set -euo pipefail;
          echo "Mapping sample: {1}";
          hisat2 -x "$INDEX_PREFIX" \
                 -1 "{2}" \
                 -2 "{3}" \
                 -p "$HISAT2_THREADS" \
                 $HISAT2_EXTRA_ARGS \
                 --no-unal \
                 2> "{5}" | \
              samtools view -@ "$SAMTOOLS_THREADS" -b -F 4 - | \
              samtools sort -@ "$SAMTOOLS_THREADS" -o "{4}" -;
          samtools index -@ "$SAMTOOLS_THREADS" "{4}"' \
         :::: "$PAIR_TSV"

echo "Done."
echo "BAM files:    $BAM_DIR"
echo "HISAT2 logs:  $LOG_DIR"
echo "Sample table: $PAIR_TSV"
echo "Job log:      $MAP_JOBLOG"
