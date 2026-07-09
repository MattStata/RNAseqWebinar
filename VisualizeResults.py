#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import warnings
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# ---------------- CLI ----------------

def get_parser():
    p = argparse.ArgumentParser(
        description=("Read an expression matrix (rows=features, columns=samples) and write PCA, optional t-SNE, "
                     "and correlation heatmaps (Pearson & Spearman) as JPEGs. "
                     "Optionally color samples by category strings (comma-separated)."))
    p.add_argument("--data", required=True, help="Expression CSV matrix")
    p.add_argument(
        "--output",
        required=True,
        help="Output prefix (writes <prefix>_PCA.jpg, optional <prefix>_TSNE.jpg, and correlation plots)"
    )
    p.add_argument(
        "--categories",
        default=None,
        help=("Comma-separated list of category strings (e.g., 'GBF1_01,GBF1_12,GFP_01,GFP_12'). "
              "A sample matches a category if that category occurs anywhere in the sample name as a whole token "
              "(i.e., delimited by non-alphanumeric characters or string boundaries). "
              "Each sample is assigned to the longest matching category.")
    )
    p.add_argument(
        "--cat_names",
        default=None,
        help=("Optional: comma-separated list of display names for categories in the legend, "
              "in the SAME ORDER as --categories. "
              "If provided, --categories MUST be provided and lengths MUST match.")
    )
    p.add_argument(
        "--min",
        type=float,
        default=None,
        help=("Optional minimum mean expression cutoff across samples. "
              "Features with mean expression < --min are removed before plotting.")
    )
    p.add_argument(
        "--top_var_pct",
        type=float,
        default=None,
        help=("Optional: keep only the top X%% most variable features after --min filtering. "
              "Variability is computed as variance across samples on log1p(expression). "
              "Example: --top_var_pct 10 keeps the top 10%% most variable features.")
    )
    p.add_argument(
        "--tsne",
        action="store_true",
        help="Enable t-SNE plot (default: off)."
    )
    p.add_argument(
        "--pcs",
        default="1,2",
        help=("Two-integer comma-separated list of which principal components to plot for PCA (1-based). "
              "Example: --pcs 1,3. Default: 1,2.")
    )
    return p

# ---------------- Data loading ----------------

def load_expression_matrix(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)

    possible_id_cols = ["GeneID", "gene_id", "transcript_id"]

    id_col = None
    for c in possible_id_cols:
        if c in df.columns:
            id_col = c
            break

    if id_col is None:
        print(
            "Error: CSV must have one of these feature ID columns: "
            + ", ".join(possible_id_cols),
            file=sys.stderr
        )
        sys.exit(1)

    sample_cols = [c for c in df.columns if c != id_col]
    if not sample_cols:
        print("Error: CSV has no sample columns.", file=sys.stderr)
        sys.exit(1)

    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df[sample_cols] = df[sample_cols].fillna(0.0)

    return df.set_index(id_col)[sample_cols].T  # Samples x features

# ---------------- Categories & colors ----------------

def _compile_category_regex(cat: str) -> re.Pattern:
    """
    Match `cat` as a whole token in a sample name.
    Token boundaries are non-alphanumeric characters (anything other than [A-Za-z0-9])
    or string boundaries.

    Prevents e.g. 'Condition_1' matching 'Condition_11' because the next char is alphanumeric.
    """
    return re.compile(rf"(?<![A-Za-z0-9]){re.escape(cat)}(?![A-Za-z0-9])")

def parse_categories(arg, sample_names):
    if not arg:
        return None, {}
    cats = [x.strip() for x in arg.split(",") if x.strip()]
    if not cats:
        return None, {}

    cat_re = {c: _compile_category_regex(c) for c in cats}

    mapping, multi_hit, no_hit = {}, [], []
    for s in sample_names:
        hits = [c for c in cats if cat_re[c].search(s)]
        if not hits:
            no_hit.append(s)
        else:
            best = max(hits, key=len)  # most specific (longest) match
            mapping[s] = best
            if len(hits) > 1:
                multi_hit.append((s, hits))

    empty_cats = [c for c in cats if c not in mapping.values()]
    if multi_hit:
        for s, hits in multi_hit:
            best = max(hits, key=len)
            print(f"Warning: sample '{s}' matches multiple categories {hits}; using '{best}'.", file=sys.stderr)
    if no_hit:
        print(f"Warning: {len(no_hit)} sample(s) did not match any category: {no_hit}", file=sys.stderr)
    if empty_cats:
        print(f"Warning: {len(empty_cats)} specified categories have no samples: {empty_cats}", file=sys.stderr)

    return cats, mapping

def colors_from_distinct(labels, order=None):
    """
    Distinct categorical colors:
    - <=10: tab10
    - <=20: tab20
    - >20: evenly spaced HSV hues

    If `order` is provided, colors and legend order follow that order for any labels present.
    """
    from matplotlib import colormaps
    from matplotlib.colors import hsv_to_rgb

    if order is not None:
        present = set(labels)
        uniq = [x for x in order if x in present]
        for x in labels:
            if x in present and x not in uniq:
                uniq.append(x)
    else:
        uniq = list(dict.fromkeys(labels))  # preserve first-seen order

    n = len(uniq)
    if n == 0:
        return [], {}

    if n <= 10:
        palette = list(colormaps["tab10"].colors[:n])
    elif n <= 20:
        palette = list(colormaps["tab20"].colors[:n])
    else:
        hues = np.linspace(0.0, 1.0, n, endpoint=False)
        palette = [tuple(hsv_to_rgb((h, 0.75, 0.90))) for h in hues]

    cmap_dict = {cat: palette[i] for i, cat in enumerate(uniq)}
    return [cmap_dict[l] for l in labels], cmap_dict

# ---------------- Feature filtering ----------------

def filter_features(expr_df: pd.DataFrame, min_mean=None, top_var_pct=None) -> pd.DataFrame:
    """
    expr_df: samples x features
    - min_mean filters features by mean expression across samples
    - top_var_pct keeps top X% features by variance across samples on log1p(expression)
    """
    if expr_df.shape[1] == 0:
        return expr_df

    if min_mean is not None:
        if min_mean < 0:
            print("Error: --min must be >= 0.", file=sys.stderr)
            sys.exit(1)
        means = expr_df.mean(axis=0)
        expr_df = expr_df.loc[:, means >= float(min_mean)]

    if top_var_pct is not None:
        pct = float(top_var_pct)
        if not (0.0 < pct <= 100.0):
            print("Error: --top_var_pct must be in (0, 100].", file=sys.stderr)
            sys.exit(1)

        if expr_df.shape[1] > 0:
            X = np.log1p(expr_df.values)
            vars_ = X.var(axis=0, ddof=0)

            k = int(np.ceil((pct / 100.0) * expr_df.shape[1]))
            if expr_df.shape[1] >= 2:
                k = max(2, min(expr_df.shape[1], k))
            else:
                k = 1

            top_idx = np.argsort(vars_)[::-1][:k]
            expr_df = expr_df.iloc[:, top_idx]

    return expr_df

# ---------------- Plot helpers ----------------

def annotate_points(ax, xs, ys, labels, fontsize=8, alpha=0.9):
    for x, y, lab in zip(xs, ys, labels):
        ax.text(x, y, lab, fontsize=fontsize, alpha=alpha)

def plot_scatter(
    xy,
    xlab,
    ylab,
    outpath,
    labels,
    category_names=None,
    legend_order=None,
    dpi=500,
    annotate_if_single_color=True
):
    n = xy.shape[0]
    rng = np.random.default_rng()
    order = rng.permutation(n)
    xs, ys = xy[order, 0], xy[order, 1]
    labels = [labels[i] for i in order]
    if category_names is not None:
        category_names = [category_names[i] for i in order]

    if category_names is None:
        colors, legend_handles = ["C0"] * n, None
    else:
        import matplotlib.patches as mpatches
        colors, cmap = colors_from_distinct(category_names, order=legend_order)
        keys = list(cmap.keys())
        if legend_order is not None:
            ordered = [k for k in legend_order if k in cmap]
            ordered += [k for k in keys if k not in ordered]
        else:
            ordered = keys
        legend_handles = [mpatches.Patch(color=cmap[k], label=k) for k in ordered]

    fig = plt.figure(figsize=(6.6, 5.2))
    ax = fig.add_subplot(111)
    ax.scatter(xs, ys, c=colors, s=40, linewidths=0.5, edgecolors="k", alpha=0.95)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.grid(True, linewidth=0.3, alpha=0.4)

    if legend_handles:
        ax.legend(handles=legend_handles, loc="best", fontsize=8, frameon=True)
    elif annotate_if_single_color:
        annotate_points(ax, xs, ys, labels, fontsize=7, alpha=0.85)

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, format="jpg", bbox_inches="tight")
    plt.close(fig)

# ---------------- Correlation clustermaps ----------------

def _tick_fontsize(n):
    if n <= 15: return 11
    if n <= 30: return 9
    if n <= 60: return 7
    if n <= 100: return 6
    return 5

def _safe_range_from_corr_mats(pear: np.ndarray, spea: np.ndarray):
    """
    Compute vmin/vmax from off-diagonal values across both matrices.
    Diagonal values are excluded because they are always 1.0.
    """
    if pear.shape[0] < 2:
        return -1.0, 1.0, 0.0

    n = pear.shape[0]
    mask = ~np.eye(n, dtype=bool)

    vals = np.concatenate([pear[mask].ravel(), spea[mask].ravel()])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return -1.0, 1.0, 0.0

    vmin = float(np.min(vals))
    vmax = float(np.max(vals))

    vmin = max(-1.0, min(1.0, vmin))
    vmax = max(-1.0, min(1.0, vmax))

    if np.isclose(vmin, vmax):
        eps = 1e-3
        vmin = max(-1.0, vmin - eps)
        vmax = min(1.0, vmax + eps)

    center = 0.0 if (vmin < 0.0 < vmax) else None
    return vmin, vmax, center

def plot_corr_heatmaps(expr_df, sample_names, outpath, dpi=500, mode="fixed"):
    """
    Pearson + Spearman clustermaps stitched into a 1x2 JPEG.

    mode:
      - "fixed": colorbar fixed to [-1, 1]
      - "fit"  : colorbar fit to min/max off-diagonal values present
    """
    try:
        from scipy.stats import rankdata
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        print("Error: scipy is required for correlation heatmaps with dendrograms (pip/conda install scipy).",
              file=sys.stderr)
        return

    X = expr_df.values  # samples x features
    n = X.shape[0]
    fs = _tick_fontsize(n)

    keep = X.var(axis=0) > 0
    X = X[:, keep]

    if X.shape[1] < 2:
        pear = np.eye(n)
        spea = np.eye(n)
    else:
        pear = np.corrcoef(X)
        ranks = np.vstack([rankdata(row, method="average") for row in X])
        spea = np.corrcoef(ranks)

    for M in (pear, spea):
        np.nan_to_num(M, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        np.fill_diagonal(M, 1.0)

    pear_df = pd.DataFrame(pear, index=sample_names, columns=sample_names)
    spea_df = pd.DataFrame(spea, index=sample_names, columns=sample_names)

    def make_linkage(corr_df):
        D = 1.0 - corr_df.values
        D = np.clip(D, 0.0, 2.0)
        np.fill_diagonal(D, 0.0)
        if np.allclose(D, 0.0):
            return None
        return linkage(squareform(D, checks=False), method="average")

    row_link_p = col_link_p = make_linkage(pear_df)
    row_link_s = col_link_s = make_linkage(spea_df)

    tree_kws = dict(color="k", linewidths=2.0)

    def _style_clustermap(cg, fs):
        for lbl in cg.ax_heatmap.get_xticklabels():
            lbl.set_rotation(90)
            lbl.set_fontsize(fs)
        for lbl in cg.ax_heatmap.get_yticklabels():
            lbl.set_fontsize(fs)

        for spine in cg.ax_heatmap.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.5)
            spine.set_edgecolor("black")

        if cg.cax is not None:
            for spine in cg.cax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.2)
                spine.set_edgecolor("black")

    cmap = "RdBu_r"

    if mode == "fixed":
        vmin, vmax, center = -1.0, 1.0, 0.0
    elif mode == "fit":
        vmin, vmax, center = _safe_range_from_corr_mats(pear, spea)
    else:
        print(f"Error: plot_corr_heatmaps mode must be 'fixed' or 'fit' (got {mode!r}).", file=sys.stderr)
        return

    cg1 = sns.clustermap(
        pear_df, cmap=cmap, vmin=vmin, vmax=vmax, center=center,
        row_cluster=True, col_cluster=True,
        row_linkage=row_link_p, col_linkage=col_link_p,
        xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.18, 0.18), colors_ratio=0.01,
        cbar_pos=(0.02, 0.80, 0.02, 0.16),
        tree_kws=tree_kws
    )
    _style_clustermap(cg1, fs)

    cg2 = sns.clustermap(
        spea_df, cmap=cmap, vmin=vmin, vmax=vmax, center=center,
        row_cluster=True, col_cluster=True,
        row_linkage=row_link_s, col_linkage=col_link_s,
        xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.18, 0.18), colors_ratio=0.01,
        cbar_pos=(0.02, 0.80, 0.02, 0.16),
        tree_kws=tree_kws
    )
    _style_clustermap(cg2, fs)

    import os
    import tempfile

    tmp_files = []
    for cg in (cg1, cg2):
        f = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        f.close()
        tmp_files.append(f.name)
        cg.savefig(f.name, dpi=dpi, bbox_inches="tight")
        plt.close(cg.fig)

    img1 = plt.imread(tmp_files[0])
    img2 = plt.imread(tmp_files[1])

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))
    axes[0].imshow(img1)
    axes[0].axis("off")
    axes[0].set_title("Pearson correlation", fontsize=13, pad=10)

    axes[1].imshow(img2)
    axes[1].axis("off")
    axes[1].set_title("Spearman (rank) correlation", fontsize=13, pad=10)

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, format="jpg", bbox_inches="tight")
    plt.close(fig)

    for f in tmp_files:
        try:
            os.remove(f)
        except OSError:
            pass

# ---------------- t-SNE perplexity ----------------

def compute_tsne_perplexity(n_samples: int) -> int:
    if n_samples <= 3:
        return 2
    p = int(round(n_samples / 3))
    p = max(5, min(30, p))
    if p >= n_samples:
        p = n_samples - 1
    return p

# ---------------- PCA pcs parsing ----------------

def parse_pcs(pcs_arg: str):
    raw = pcs_arg.strip()
    parts = [x.strip() for x in raw.split(",") if x.strip()]
    if len(parts) != 2:
        print("Error: --pcs must be exactly two comma-separated integers (e.g., 1,2 or 1,3).", file=sys.stderr)
        sys.exit(1)
    try:
        a = int(parts[0])
        b = int(parts[1])
    except ValueError:
        print("Error: --pcs must be integers (e.g., 1,2).", file=sys.stderr)
        sys.exit(1)
    if a < 1 or b < 1:
        print("Error: --pcs values must be >= 1 (1-based component numbers).", file=sys.stderr)
        sys.exit(1)
    if a == b:
        print("Error: --pcs must specify two different component numbers.", file=sys.stderr)
        sys.exit(1)
    return a, b

# ---------------- Main ----------------

def main():
    parser = get_parser()
    args = parser.parse_args()

    # Validate category display names
    if args.cat_names is not None and not args.categories:
        print("Error: --cat_names was provided but --categories was not provided.", file=sys.stderr)
        sys.exit(1)

    pc_a, pc_b = parse_pcs(args.pcs)
    n_pca = max(pc_a, pc_b)

    # Load data as samples x features
    expr = load_expression_matrix(args.data)
    sample_names = list(expr.index)

    # Feature filtering
    expr = filter_features(expr, min_mean=args.min, top_var_pct=args.top_var_pct)

    print(f"Features_used_for_plotting\t{expr.shape[1]}")

    if expr.shape[1] < 2:
        print("Error: after filtering, fewer than 2 features remain; cannot run PCA.", file=sys.stderr)
        sys.exit(1)

    # Transform and scale
    X = np.log1p(expr.values)  # samples x features
    Xz = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    if Xz.shape[0] < 2:
        print("Error: need at least 2 samples for PCA.", file=sys.stderr)
        sys.exit(1)

    max_possible_pcs = min(Xz.shape[0], Xz.shape[1])
    if n_pca > max_possible_pcs:
        print(
            f"Error: requested --pcs {pc_a},{pc_b} requires computing up to PC{n_pca}, but only "
            f"{max_possible_pcs} PCs are available (min(n_samples, n_features)).",
            file=sys.stderr
        )
        sys.exit(1)

    # Categories and optional display names
    cats, mapping = parse_categories(args.categories, sample_names) if args.categories else (None, {})

    cat_display = None
    legend_order = None
    if cats is not None:
        if args.cat_names is not None:
            cat_names = [x.strip() for x in args.cat_names.split(",") if x.strip()]
            if len(cat_names) != len(cats):
                print(
                    f"Error: --cat_names length ({len(cat_names)}) does not match "
                    f"--categories length ({len(cats)}).",
                    file=sys.stderr
                )
                sys.exit(1)
            cat_display = dict(zip(cats, cat_names))
            legend_order = cat_names[:]
        else:
            cat_display = {c: c for c in cats}
            legend_order = cats[:]

    cat_labels = None
    if cats is not None:
        cat_labels = []
        for s in sample_names:
            c = mapping.get(s, None)
            if c is None:
                cat_labels.append("Uncategorized")
            else:
                cat_labels.append(cat_display[c])

        if "Uncategorized" in cat_labels:
            print("Warning: Some samples were not assigned; labeled as 'Uncategorized'.", file=sys.stderr)

        if legend_order is not None and "Uncategorized" in set(cat_labels) and "Uncategorized" not in legend_order:
            legend_order = legend_order + ["Uncategorized"]

    # PCA
    pca = PCA(n_components=n_pca, random_state=42)
    pca_scores = pca.fit_transform(Xz)
    pca_xy = pca_scores[:, [pc_a - 1, pc_b - 1]]

    pc_a_var = pca.explained_variance_ratio_[pc_a - 1] * 100.0
    pc_b_var = pca.explained_variance_ratio_[pc_b - 1] * 100.0

    pca_out = f"{args.output}_PCA.jpg"
    plot_scatter(
        pca_xy,
        xlab=f"PC{pc_a} ({pc_a_var:.1f}%)",
        ylab=f"PC{pc_b} ({pc_b_var:.1f}%)",
        outpath=pca_out,
        labels=sample_names,
        category_names=cat_labels,
        legend_order=legend_order
    )
    print(f"Wrote {pca_out}", file=sys.stderr)

    # t-SNE
    if args.tsne:
        perplexity = compute_tsne_perplexity(Xz.shape[0])
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            init="pca",
            learning_rate="auto",
            random_state=42
        )
        tsne_xy = tsne.fit_transform(Xz)
        tsne_out = f"{args.output}_TSNE.jpg"
        plot_scatter(
            tsne_xy,
            xlab="TSNE1",
            ylab="TSNE2",
            outpath=tsne_out,
            labels=sample_names,
            category_names=cat_labels,
            legend_order=legend_order
        )
        print(f"Wrote {tsne_out}", file=sys.stderr)

    # Correlation clustermaps
    corr_out_fixed = f"{args.output}_CORR.jpg"
    plot_corr_heatmaps(
        expr_df=pd.DataFrame(X, index=sample_names),
        sample_names=sample_names,
        outpath=corr_out_fixed,
        mode="fixed"
    )
    print(f"Wrote {corr_out_fixed}", file=sys.stderr)

    corr_out_fit = f"{args.output}_CORR_FIT.jpg"
    plot_corr_heatmaps(
        expr_df=pd.DataFrame(X, index=sample_names),
        sample_names=sample_names,
        outpath=corr_out_fit,
        mode="fit"
    )
    print(f"Wrote {corr_out_fit}", file=sys.stderr)

if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()
