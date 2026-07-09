#!/bin/bash --login
#SBATCH --job-name=trimmomatic_fastq
#SBATCH --output=trimmomatic_fastq_%j.out
#SBATCH --error=trimmomatic_fastq_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

set -euo pipefail

###############################################################################
# Edit this block for your own project
###############################################################################

TRIMMOMATIC_MODULE="Trimmomatic/0.39-Java-17"
PARALLEL_MODULE="parallel/20230722-GCCcore-12.2.0"

FASTQ_DIR="FASTQ"
OUTPUT_DIR="trimmomatic_output"

# Leave ADAPTER_FASTA blank to use a bundled Trimmomatic adapter file.
# Common bundled choices include TruSeq3-PE.fa, TruSeq2-PE.fa, and NexteraPE-PE.fa.
ADAPTER_FILE_NAME="TruSeq3-PE.fa"
ADAPTER_FASTA=""

# Parallelization: PARALLEL_JOBS samples at a time, each using TRIMMOMATIC_THREADS.
# Keep PARALLEL_JOBS * TRIMMOMATIC_THREADS <= SLURM_CPUS_PER_TASK.
PARALLEL_JOBS=8
TRIMMOMATIC_THREADS=4
JAVA_XMX="6g"

# Conservative starting settings for Illumina paired-end RNA-seq/short-read data.
SEED_MISMATCHES=2
PALINDROME_CLIP_THRESHOLD=30
SIMPLE_CLIP_THRESHOLD=10
MIN_ADAPTER_LENGTH=2
KEEP_BOTH_READS=true
TRIM_STEPS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"

###############################################################################
# End edit block
###############################################################################

module purge
module load "$TRIMMOMATIC_MODULE"
module load "$PARALLEL_MODULE"

TRIMMOMATIC_JAR="${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar"

if [[ -z "$ADAPTER_FASTA" ]]; then
    ADAPTER_FASTA="${EBROOTTRIMMOMATIC}/adapters/${ADAPTER_FILE_NAME}"
fi

if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    echo "ERROR: Trimmomatic jar not found: $TRIMMOMATIC_JAR" >&2
    exit 1
fi

if [[ ! -f "$ADAPTER_FASTA" ]]; then
    echo "ERROR: Adapter FASTA not found: $ADAPTER_FASTA" >&2
    echo "Bundled adapter files available here:" >&2
    ls -1 "${EBROOTTRIMMOMATIC}/adapters" >&2 || true
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

PAIRED_DIR="${OUTPUT_DIR}/paired"
UNPAIRED_DIR="${OUTPUT_DIR}/unpaired"
LOG_DIR="${OUTPUT_DIR}/logs"
PAIR_LIST="${OUTPUT_DIR}/fastq_pairs.tsv"
JOBLOG="${OUTPUT_DIR}/parallel_joblog.tsv"

mkdir -p "$PAIRED_DIR" "$UNPAIRED_DIR" "$LOG_DIR"
: > "$PAIR_LIST"

shopt -s nullglob
R1_FILES=("${FASTQ_DIR}"/*_R1.fq.gz "${FASTQ_DIR}"/*_R1.fastq.gz)
shopt -u nullglob

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 files found in $FASTQ_DIR" >&2
    echo "Expected names like sample_R1.fq.gz or sample_R1.fastq.gz" >&2
    exit 1
fi

missing_pairs=0
for r1 in "${R1_FILES[@]}"; do
    if [[ "$r1" == *_R1.fq.gz ]]; then
        stem="${r1%_R1.fq.gz}"
        r2="${stem}_R2.fq.gz"
    elif [[ "$r1" == *_R1.fastq.gz ]]; then
        stem="${r1%_R1.fastq.gz}"
        r2="${stem}_R2.fastq.gz"
    else
        continue
    fi

    sample="$(basename "$stem")"

    if [[ -f "$r2" ]]; then
        printf '%s\t%s\t%s\n' "$r1" "$r2" "$sample" >> "$PAIR_LIST"
    else
        echo "Missing R2 for: $r1" >&2
        missing_pairs=$((missing_pairs + 1))
    fi
done

if [[ "$missing_pairs" -gt 0 ]]; then
    echo "ERROR: Found $missing_pairs R1 file(s) without matching R2 files." >&2
    exit 1
fi

N_PAIRS="$(wc -l < "$PAIR_LIST")"

if [[ "$N_PAIRS" -eq 0 ]]; then
    echo "ERROR: No complete FASTQ pairs found." >&2
    exit 1
fi

if (( PARALLEL_JOBS * TRIMMOMATIC_THREADS > SLURM_CPUS_PER_TASK )); then
    echo "ERROR: PARALLEL_JOBS * TRIMMOMATIC_THREADS exceeds requested CPUs." >&2
    echo "PARALLEL_JOBS=$PARALLEL_JOBS" >&2
    echo "TRIMMOMATIC_THREADS=$TRIMMOMATIC_THREADS" >&2
    echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK" >&2
    exit 1
fi

echo "FASTQ directory: $FASTQ_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Adapter FASTA: $ADAPTER_FASTA"
echo "FASTQ pairs found: $N_PAIRS"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Threads per Trimmomatic job: $TRIMMOMATIC_THREADS"
echo "Started: $(date)"

export TRIMMOMATIC_JAR ADAPTER_FASTA PAIRED_DIR UNPAIRED_DIR LOG_DIR
export TRIMMOMATIC_THREADS JAVA_XMX
export SEED_MISMATCHES PALINDROME_CLIP_THRESHOLD SIMPLE_CLIP_THRESHOLD
export MIN_ADAPTER_LENGTH KEEP_BOTH_READS TRIM_STEPS
export SHELL=/bin/bash

parallel --colsep '\t' \
    --jobs "$PARALLEL_JOBS" \
    --joblog "$JOBLOG" \
    --halt soon,fail=1 \
    'java -Xmx"${JAVA_XMX}" -jar "${TRIMMOMATIC_JAR}" PE \
        -threads "${TRIMMOMATIC_THREADS}" \
        -phred33 \
        "{1}" "{2}" \
        "${PAIRED_DIR}/{3}_R1.trimmed.paired.fq.gz" \
        "${UNPAIRED_DIR}/{3}_R1.trimmed.unpaired.fq.gz" \
        "${PAIRED_DIR}/{3}_R2.trimmed.paired.fq.gz" \
        "${UNPAIRED_DIR}/{3}_R2.trimmed.unpaired.fq.gz" \
        "ILLUMINACLIP:${ADAPTER_FASTA}:${SEED_MISMATCHES}:${PALINDROME_CLIP_THRESHOLD}:${SIMPLE_CLIP_THRESHOLD}:${MIN_ADAPTER_LENGTH}:${KEEP_BOTH_READS}" \
        ${TRIM_STEPS} \
        > "${LOG_DIR}/{3}.trimmomatic.log" 2>&1' \
    :::: "$PAIR_LIST"

N_PAIRED_R1="$(find "$PAIRED_DIR" -name '*_R1.trimmed.paired.fq.gz' | wc -l)"
N_PAIRED_R2="$(find "$PAIRED_DIR" -name '*_R2.trimmed.paired.fq.gz' | wc -l)"

echo "Finished: $(date)"
echo "Trimmed paired R1 files: $N_PAIRED_R1"
echo "Trimmed paired R2 files: $N_PAIRED_R2"
echo "Paired output: $PAIRED_DIR"
echo "Unpaired output: $UNPAIRED_DIR"
echo "Logs: $LOG_DIR"
echo "Parallel job log: $JOBLOG"
