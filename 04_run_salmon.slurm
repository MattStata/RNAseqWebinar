#!/bin/bash --login
#SBATCH --job-name=salmon_gff_quant
#SBATCH --output=salmon_gff_quant_%j.out
#SBATCH --error=salmon_gff_quant_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=256G

set -euo pipefail

###############################################################################
# Edit this block for your own data and cluster environment
###############################################################################

SALMON_MODULE="Salmon/1.10.1-GCC-12.3.0"
PARALLEL_MODULE="parallel/20230722-GCCcore-12.3.0"

FASTQ_DIR="FASTQ"
GENOME_FASTA="Tidestromia_oblongifolia_DV.genome.assembly.softmasked.fasta"
GFF3_FILE="Tidestromia_oblongifolia_DV.gene_models.gff3"

OUTDIR="salmon_gff_quant_output"

PARALLEL_JOBS=2
SALMON_THREADS=16
INDEX_THREADS=16

###############################################################################
# Setup
###############################################################################

module purge
module load "$SALMON_MODULE"
module load "$PARALLEL_MODULE"

REFERENCE_DIR="${OUTDIR}/reference"
INDEX_DIR="${OUTDIR}/salmon_index"
QUANT_DIR="${OUTDIR}/quants"
LOG_DIR="${OUTDIR}/logs"
SCRIPT_DIR="${OUTDIR}/sample_scripts"

TRANSCRIPTS_FASTA="${REFERENCE_DIR}/transcripts.fa"
TX2GENE="${REFERENCE_DIR}/tx2gene.tsv"
SAMPLE_TABLE="${OUTDIR}/sample_pairs.tsv"
PARALLEL_JOBLOG="${OUTDIR}/parallel_joblog.tsv"

mkdir -p "$REFERENCE_DIR" "$INDEX_DIR" "$QUANT_DIR" "$LOG_DIR" "$SCRIPT_DIR"

command -v salmon >/dev/null 2>&1 || { echo "ERROR: salmon not found after loading modules." >&2; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: parallel not found after loading modules." >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found." >&2; exit 1; }

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

if [[ ! -s "$GENOME_FASTA" ]]; then
    echo "ERROR: genome FASTA not found or empty: $GENOME_FASTA" >&2
    exit 1
fi

if [[ ! -s "$GFF3_FILE" ]]; then
    echo "ERROR: GFF3 file not found or empty: $GFF3_FILE" >&2
    exit 1
fi

###############################################################################
# Build transcript FASTA and transcript-to-gene map from the genome and GFF3
###############################################################################

export GENOME_FASTA GFF3_FILE TRANSCRIPTS_FASTA TX2GENE

python3 - <<'PY'
from pathlib import Path
from collections import defaultdict
import os
import sys

GENOME_FASTA = Path(os.environ["GENOME_FASTA"])
GFF3_FILE = Path(os.environ["GFF3_FILE"])
TRANSCRIPTS_FASTA = Path(os.environ["TRANSCRIPTS_FASTA"])
TX2GENE = Path(os.environ["TX2GENE"])

def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs

def parse_attrs(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs

def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

print(f"Reading genome FASTA: {GENOME_FASTA}")
genome = read_fasta(GENOME_FASTA)
print(f"Loaded {len(genome)} genome sequences")

transcript_to_gene = {}
features_by_tx = defaultdict(list)
fallback_cds_by_tx = defaultdict(list)

with GFF3_FILE.open() as handle:
    for line in handle:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            if line == "##FASTA":
                break
            continue

        fields = line.split("\t")
        if len(fields) != 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attrs_raw = fields
        attrs = parse_attrs(attrs_raw)
        start = int(start)
        end = int(end)

        if feature in {"mRNA", "transcript"}:
            tx_id = attrs.get("ID")
            gene_id = attrs.get("Parent", attrs.get("gene", tx_id))
            if tx_id:
                gene_id = gene_id.split(",")[0]
                transcript_to_gene[tx_id] = gene_id

        elif feature in {"exon", "CDS"}:
            parents = attrs.get("Parent")
            if not parents:
                continue
            for tx_id in parents.split(","):
                if feature == "exon":
                    features_by_tx[tx_id].append((seqid, start, end, strand))
                else:
                    fallback_cds_by_tx[tx_id].append((seqid, start, end, strand))

if not features_by_tx and fallback_cds_by_tx:
    print("No exon features were found; using CDS features instead.")
    features_by_tx = fallback_cds_by_tx

written = 0
missing_seqids = set()
missing_gene_ids = 0

with TRANSCRIPTS_FASTA.open("w") as fasta_out, TX2GENE.open("w") as map_out:
    for tx_id in sorted(features_by_tx):
        parts = features_by_tx[tx_id]
        if not parts:
            continue

        seqids = {p[0] for p in parts}
        strands = {p[3] for p in parts}
        if len(seqids) != 1 or len(strands) != 1:
            continue

        seqid = next(iter(seqids))
        strand = next(iter(strands))
        if seqid not in genome:
            missing_seqids.add(seqid)
            continue

        if strand == "+":
            parts = sorted(parts, key=lambda x: x[1])
        elif strand == "-":
            parts = sorted(parts, key=lambda x: x[1], reverse=True)
        else:
            parts = sorted(parts, key=lambda x: x[1])

        seq = "".join(genome[seqid][start - 1:end] for _, start, end, _ in parts)
        if strand == "-":
            seq = revcomp(seq)

        if not seq:
            continue

        gene_id = transcript_to_gene.get(tx_id)
        if gene_id is None:
            gene_id = tx_id.split(".")[0]
            missing_gene_ids += 1

        fasta_out.write(f">{tx_id} gene={gene_id}\n")
        for i in range(0, len(seq), 80):
            fasta_out.write(seq[i:i+80] + "\n")

        map_out.write(f"{tx_id}\t{gene_id}\n")
        written += 1

if missing_seqids:
    example = ", ".join(sorted(missing_seqids)[:5])
    print(f"ERROR: {len(missing_seqids)} GFF3 sequence IDs were not found in the genome FASTA.", file=sys.stderr)
    print(f"Examples: {example}", file=sys.stderr)
    sys.exit(1)

if written == 0:
    print("ERROR: no transcript sequences were written.", file=sys.stderr)
    sys.exit(1)

print(f"Wrote transcript FASTA: {TRANSCRIPTS_FASTA}")
print(f"Wrote transcript-to-gene map: {TX2GENE}")
print(f"Transcript sequences written: {written}")
if missing_gene_ids:
    print(f"Transcripts without explicit gene parent, using ID prefix fallback: {missing_gene_ids}")
PY

###############################################################################
# Build Salmon index
###############################################################################

salmon index \
    -t "$TRANSCRIPTS_FASTA" \
    -i "$INDEX_DIR" \
    -p "$INDEX_THREADS"

###############################################################################
# Find paired FASTQ files
###############################################################################

: > "$SAMPLE_TABLE"

find "$FASTQ_DIR" -maxdepth 1 -type f \( -name "*_R1.fq.gz" -o -name "*_R1.fastq.gz" \) | sort | while read -r R1; do
    R2="${R1/_R1.fq.gz/_R2.fq.gz}"
    R2="${R2/_R1.fastq.gz/_R2.fastq.gz}"

    if [[ ! -s "$R2" ]]; then
        echo "ERROR: missing R2 file for $R1" >&2
        exit 1
    fi

    SAMPLE=$(basename "$R1")
    SAMPLE="${SAMPLE%_R1.fq.gz}"
    SAMPLE="${SAMPLE%_R1.fastq.gz}"

    printf "%s\t%s\t%s\n" "$SAMPLE" "$R1" "$R2" >> "$SAMPLE_TABLE"
done

N_SAMPLES=$(wc -l < "$SAMPLE_TABLE")
if [[ "$N_SAMPLES" -eq 0 ]]; then
    echo "ERROR: no paired FASTQ files found in $FASTQ_DIR" >&2
    exit 1
fi

echo "Found $N_SAMPLES paired-end samples."

###############################################################################
# Create one small Salmon command script per sample
###############################################################################

rm -f "${SCRIPT_DIR}"/*.sh

while IFS=$'\t' read -r SAMPLE R1 R2; do
    SAMPLE_SCRIPT="${SCRIPT_DIR}/${SAMPLE}.salmon.sh"
    SAMPLE_OUTDIR="${QUANT_DIR}/${SAMPLE}"
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE}.salmon.log"

    cat > "$SAMPLE_SCRIPT" <<EOS
#!/bin/bash --login
set -euo pipefail

salmon quant \\
    -i "$INDEX_DIR" \\
    -l A \\
    -1 "$R1" \\
    -2 "$R2" \\
    -p "$SALMON_THREADS" \\
    --validateMappings \\
    --seqBias \\
    --gcBias \\
    --geneMap "$TX2GENE" \\
    -o "$SAMPLE_OUTDIR" \\
    > "$SAMPLE_LOG" 2>&1
EOS

    chmod +x "$SAMPLE_SCRIPT"
done < "$SAMPLE_TABLE"

###############################################################################
# Quantify samples
###############################################################################

find "$SCRIPT_DIR" -maxdepth 1 -type f -name "*.salmon.sh" | sort | \
    parallel \
        --jobs "$PARALLEL_JOBS" \
        --joblog "$PARALLEL_JOBLOG" \
        --halt soon,fail=1 \
        bash {}

###############################################################################
# Check outputs
###############################################################################

MISSING=0
while IFS=$'\t' read -r SAMPLE R1 R2; do
    if [[ ! -s "${QUANT_DIR}/${SAMPLE}/quant.sf" ]]; then
        echo "ERROR: missing quant.sf for $SAMPLE" >&2
        MISSING=1
    fi
    if [[ ! -s "${QUANT_DIR}/${SAMPLE}/quant.genes.sf" ]]; then
        echo "ERROR: missing quant.genes.sf for $SAMPLE" >&2
        MISSING=1
    fi
done < "$SAMPLE_TABLE"

if [[ "$MISSING" -ne 0 ]]; then
    exit 1
fi

echo "Salmon quantification complete."
echo "Transcript-level quantification: ${QUANT_DIR}/*/quant.sf"
echo "Gene-level quantification:       ${QUANT_DIR}/*/quant.genes.sf"
echo "Transcript-to-gene map:          $TX2GENE"
