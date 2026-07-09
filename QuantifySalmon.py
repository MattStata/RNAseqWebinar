#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set


def read_salmon_sf(path: Path, metric: str) -> Dict[str, str]:
    """
    Read a Salmon quant.sf / quant.genes.sf file and return:
        {Name: metric_value_as_string}
    """
    out: Dict[str, str] = {}

    with path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n\r").split("\t")

        try:
            name_i = header.index("Name")
        except ValueError:
            raise SystemExit(f"ERROR: 'Name' column not found in {path}")

        try:
            metric_i = header.index(metric)
        except ValueError:
            raise SystemExit(
                f"ERROR: metric column '{metric}' not found in {path}\n"
                f"       Available columns: {', '.join(header)}"
            )

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n\r").split("\t")
            if len(fields) <= max(name_i, metric_i):
                continue

            out[fields[name_i]] = fields[metric_i]

    return out


def detect_samples(quant_root: Path) -> List[str]:
    """
    Auto-detect sample directories under the quant root.
    """
    return sorted([p.name for p in quant_root.iterdir() if p.is_dir()])


def derive_prefix_from_quant_root(quant_root: Path) -> str:
    """
    Derive a sensible output prefix from the quant root directory name.
    For example:
        Tidestromia_oblongifolia_Salmon_Quantification_salmon_quant
    becomes:
        Tidestromia_oblongifolia_Salmon_Quantification
    """
    name = quant_root.name
    suffix = "_salmon_quant"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def write_matrix(
    out_csv: Path,
    id_header: str,
    samples: List[str],
    per_sample: Dict[str, Dict[str, str]],
    all_ids: Set[str],
    missing_value: str,
) -> None:
    """
    Write a wide CSV matrix with IDs as rows and samples as columns.
    """
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([id_header] + samples)

        for feature_id in sorted(all_ids):
            row = [feature_id] + [per_sample[s].get(feature_id, missing_value) for s in samples]
            w.writerow(row)


def build_one_matrix(
    quant_root: Path,
    samples: List[str],
    metric: str,
    which: str,
    outdir: Path,
    prefix: str,
    missing_value: str,
) -> None:
    """
    Build one or both matrix types for a single Salmon metric.
    """
    # Transcript-level matrix
    if which in ("both", "tx"):
        per_sample_tx: Dict[str, Dict[str, str]] = {}
        all_tx_ids: Set[str] = set()

        for sample in samples:
            sf = quant_root / sample / "quant.sf"
            if not sf.exists():
                raise SystemExit(f"ERROR: missing transcript-level quant file: {sf}")

            d = read_salmon_sf(sf, metric)
            per_sample_tx[sample] = d
            all_tx_ids.update(d.keys())

        out_csv = outdir / f"{prefix}.transcripts.{metric}.csv"

        write_matrix(
            out_csv=out_csv,
            id_header="transcript_id",
            samples=samples,
            per_sample=per_sample_tx,
            all_ids=all_tx_ids,
            missing_value=missing_value,
        )

        print(f"Wrote transcript matrix: {out_csv}")

    # Gene-level matrix
    if which in ("both", "gene"):
        per_sample_g: Dict[str, Dict[str, str]] = {}
        all_g_ids: Set[str] = set()

        for sample in samples:
            sf = quant_root / sample / "quant.genes.sf"
            if not sf.exists():
                raise SystemExit(f"ERROR: missing gene-level quant file: {sf}")

            d = read_salmon_sf(sf, metric)
            per_sample_g[sample] = d
            all_g_ids.update(d.keys())

        out_csv = outdir / f"{prefix}.genes.{metric}.csv"

        write_matrix(
            out_csv=out_csv,
            id_header="gene_id",
            samples=samples,
            per_sample=per_sample_g,
            all_ids=all_g_ids,
            missing_value=missing_value,
        )

        print(f"Wrote gene matrix: {out_csv}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Build wide CSV expression matrices from Salmon outputs.\n"
            "By default writes both transcript-level (quant.sf) and gene-level "
            "(quant.genes.sf) matrices for both TPM and NumReads."
        )
    )

    ap.add_argument(
        "--quant-root",
        required=True,
        help="Directory containing per-sample Salmon output subdirectories",
    )
    ap.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help=(
            "Optional explicit sample order as subdirectory names under --quant-root. "
            "If omitted, sample directories are auto-detected and sorted alphabetically."
        ),
    )
    ap.add_argument(
        "--metric",
        choices=["Both", "TPM", "NumReads"],
        default="Both",
        help=(
            "Which Salmon metric to extract: Both, TPM, or NumReads "
            "(default: Both)"
        ),
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help=(
            "Output directory for CSV matrices. "
            "Default: derived from the quant root as <prefix>_salmon_matrices"
        ),
    )
    ap.add_argument(
        "--prefix",
        default=None,
        help=(
            "Prefix for output filenames. "
            "Default: derived from the quant root directory name"
        ),
    )
    ap.add_argument(
        "--which",
        choices=["both", "tx", "gene"],
        default="both",
        help="Which matrices to write: both, tx, or gene (default: both)",
    )
    ap.add_argument(
        "--missing",
        default="0",
        help="Value to write for IDs missing from a sample (default: 0)",
    )

    args = ap.parse_args()

    quant_root = Path(args.quant_root)

    if not quant_root.exists():
        raise SystemExit(f"ERROR: quant root does not exist: {quant_root}")

    if not quant_root.is_dir():
        raise SystemExit(f"ERROR: quant root is not a directory: {quant_root}")

    prefix = args.prefix if args.prefix is not None else derive_prefix_from_quant_root(quant_root)
    outdir = Path(args.outdir) if args.outdir is not None else Path(f"{prefix}_salmon_matrices")
    outdir.mkdir(parents=True, exist_ok=True)

    samples = args.samples if args.samples is not None else detect_samples(quant_root)

    if not samples:
        raise SystemExit(f"ERROR: no sample subdirectories found under {quant_root}")

    print(f"Quant root: {quant_root}")
    print(f"Detected samples ({len(samples)}): {', '.join(samples)}")
    print(f"Metric mode: {args.metric}")
    print(f"Output directory: {outdir}")
    print(f"Output prefix: {prefix}")
    print(f"Writing: {args.which}")

    metrics_to_write = ["TPM", "NumReads"] if args.metric == "Both" else [args.metric]

    for metric in metrics_to_write:
        print(f"\nProcessing metric: {metric}")
        build_one_matrix(
            quant_root=quant_root,
            samples=samples,
            metric=metric,
            which=args.which,
            outdir=outdir,
            prefix=prefix,
            missing_value=args.missing,
        )

    print("Done.")


if __name__ == "__main__":
    main()