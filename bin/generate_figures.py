#!/usr/bin/env python3
"""
bin/generate_figures.py

Author: Nadeem Khan, Bioinformatician, INRS-CAFSB
"""

import argparse, gzip, csv, subprocess, sys, warnings
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
import matplotlib.patheffects as pe
from scipy import stats
warnings.filterwarnings("ignore")

# ── Nature-style rcParams ──────────────────────────────────────────
plt.rcParams.update({
    "font.family"        : "sans-serif",
    "font.sans-serif"    : ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size"          : 7,
    "axes.titlesize"     : 8,
    "axes.labelsize"     : 7,
    "xtick.labelsize"    : 6.5,
    "ytick.labelsize"    : 6.5,
    "legend.fontsize"    : 6.5,
    "figure.dpi"         : 300,
    "savefig.dpi"        : 300,
    "savefig.bbox"       : "tight",
    "savefig.pad_inches" : 0.05,
    "axes.linewidth"     : 0.6,
    "xtick.major.width"  : 0.6,
    "ytick.major.width"  : 0.6,
    "xtick.major.size"   : 3,
    "ytick.major.size"   : 3,
    "xtick.minor.size"   : 1.5,
    "ytick.minor.size"   : 1.5,
    "lines.linewidth"    : 0.8,
    "patch.linewidth"    : 0.5,
    "axes.spines.top"    : False,
    "axes.spines.right"  : False,
    "pdf.fonttype"       : 42,   # embeds fonts in PDF
    "ps.fonttype"        : 42,
})

# ── Colour palettes ────────────────────────────────────────────────
CHR_COLORS  = ["#2166ac", "#4d9de0"]            # alternating CHR colours
GW_RED      = "#d62728"                          # GW-sig hits
SUG_AMBER   = "#ff7f0e"                          # suggestive threshold
COHORT_COL  = {"Cohort_A": "#2166ac", "Cohort_B": "#c0392b"}

# 1000G superpopulation colours (publication standard)
SUPERPOP_COL = {
    "AFR": "#d62728", "AMR": "#2ca02c", "EAS": "#1f77b4",
    "EUR": "#9467bd", "SAS": "#ff7f0e"
}
# Fine-grained population colours
POP_COL = {
    # AFR
    "YRI":"#d62728","LWK":"#e74c3c","GWD":"#c0392b","MSL":"#e84393","ESN":"#922b21",
    "ASW":"#f1948a","ACB":"#a93226",
    # AMR
    "MXL":"#2ca02c","PUR":"#27ae60","CLM":"#1e8449","PEL":"#52be80",
    # EAS
    "CHB":"#1f77b4","JPT":"#2e86c1","CHS":"#3498db","CDX":"#1a5276","KHV":"#5dade2",
    # EUR
    "CEU":"#9467bd","TSI":"#7d3c98","FIN":"#a569bd","GBR":"#6c3483","IBS":"#d7bde2",
    # SAS
    "GIH":"#ff7f0e","PJL":"#e67e22","BEB":"#d35400","STU":"#f0b27a","ITU":"#fdebd0",
}

# Fine-grained code → super-population lookup (all 26 1000G populations)
POP_TO_SPOP = {
    # AFR
    "YRI":"AFR","LWK":"AFR","GWD":"AFR","MSL":"AFR","ESN":"AFR","ASW":"AFR","ACB":"AFR",
    # AMR
    "MXL":"AMR","PUR":"AMR","CLM":"AMR","PEL":"AMR",
    # EAS
    "CHB":"EAS","JPT":"EAS","CHS":"EAS","CDX":"EAS","KHV":"EAS",
    # EUR
    "CEU":"EUR","TSI":"EUR","FIN":"EUR","GBR":"EUR","IBS":"EUR",
    # SAS
    "GIH":"SAS","PJL":"SAS","BEB":"SAS","STU":"SAS","ITU":"SAS",
}

TRAIT_COL = {
    "ldl_cholesterol": "#1f77b4",
    "bmi":             "#ff7f0e",
    "crp_log":         "#2ca02c",
    "height_cm":       "#d62728",
    "cad":             "#9467bd",
    "t2d":             "#8c564b",
    "hypertension":    "#e377c2",
}
TRAIT_LABEL = {
    "ldl_cholesterol": "LDL cholesterol",
    "bmi":             "BMI",
    "crp_log":         "CRP (log)",
    "height_cm":       "Height",
    "cad":             "CAD",
    "t2d":             "T2D",
    "hypertension":    "Hypertension",
}

# ── Helpers ────────────────────────────────────────────────────────
def save(fig, out_dir, name, formats=("png","pdf")):
    for fmt in formats:
        path = Path(out_dir) / f"{name}.{fmt}"
        fig.savefig(path, format=fmt)
        print(f"  Saved: {path}")

def load_sumstat(gz_path, pval_col="P", max_pval=0.1):
    """Load a harmonised sumstat gz, returning DataFrame filtered to p < max_pval."""
    chunks = []
    with gzip.open(gz_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                if float(row.get(pval_col,"1")) <= max_pval:
                    chunks.append(row)
            except: pass
    if not chunks:
        return pd.DataFrame()
    df = pd.DataFrame(chunks)
    for c in ["CHR","BP","BETA","SE","P","LOG10P","EAF","N"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["CHR","BP","P"])

def compute_manhattan_positions(df):
    """Add absolute genome position column for Manhattan plot."""
    chrom_sizes = (
        df.groupby("CHR")["BP"].max()
          .sort_index()
          .cumsum()
          .shift(1, fill_value=0)
    )
    df = df.copy()
    df["ABS_POS"] = df["BP"] + df["CHR"].map(chrom_sizes)
    df["CHR_MID"] = df["CHR"].map(
        df.groupby("CHR")["ABS_POS"].mean()
    )
    return df, chrom_sizes

def gc_lambda(pvals):
    """Compute genomic inflation factor λGC."""
    chi2 = stats.chi2.ppf(1 - np.array(pvals), df=1)
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))

def qq_ci(n):
    """95% confidence interval for QQ plot (beta distribution)."""
    ranks = np.arange(1, n+1)
    ci_lo = -np.log10(stats.beta.ppf(0.975, ranks, n-ranks+1))
    ci_hi = -np.log10(stats.beta.ppf(0.025, ranks, n-ranks+1))
    return ci_lo, ci_hi

# ══════════════════════════════════════════════════════════════════
# F1 — Miami Plot
# ══════════════════════════════════════════════════════════════════
def figure_miami(sumstat_dir, traits, cohorts, out_dir, author, institute):
    print("[F1] Miami plot...")
    fig, axes = plt.subplots(2, 1, figsize=(18/2.54, 12/2.54),
                              gridspec_kw={"hspace":0.08})

    GW     = 5e-8
    SUGG   = 1e-5
    n_plotted = {c:0 for c in cohorts}

    for ax_i, cohort in enumerate(cohorts):
        ax = axes[ax_i]
        all_rows = []
        for trait in traits:
            gz = Path(sumstat_dir) / cohort / "sumstats" / f"{cohort}_{trait}_sumstats.tsv.gz"
            if not gz.exists():
                print(f"  [WARN] Missing: {gz}")
                continue
            df = load_sumstat(gz, max_pval=0.05)
            if df.empty: continue
            df["TRAIT"] = trait
            all_rows.append(df)

        if not all_rows:
            ax.text(0.5, 0.5, f"No data for {cohort}",
                    ha="center", va="center", transform=ax.transAxes, fontsize=7)
            continue

        data = pd.concat(all_rows, ignore_index=True)
        data = data.dropna(subset=["CHR","BP","P"])
        data["CHR"] = data["CHR"].astype(int)
        data, chrom_off = compute_manhattan_positions(data)

        # Plot chromosomes alternating colours
        chrs = sorted(data["CHR"].unique())
        for ci, chrom in enumerate(chrs):
            sub = data[data["CHR"] == chrom].copy()
            sub["-log10p"] = -np.log10(sub["P"].clip(lower=1e-300))
            col = CHR_COLORS[ci % 2]

            # Non-significant in small dots
            ns = sub[sub["P"] > SUGG]
            if len(ns):
                ax.scatter(ns["ABS_POS"], ns["-log10p"],
                           c=col, s=0.8, linewidths=0, rasterized=True, alpha=0.6)

            # Suggestive
            sg = sub[(sub["P"] <= SUGG) & (sub["P"] > GW)]
            if len(sg):
                ax.scatter(sg["ABS_POS"], sg["-log10p"],
                           c=SUG_AMBER, s=3, linewidths=0, rasterized=True)

            # GW significant
            gw = sub[sub["P"] <= GW]
            if len(gw):
                ax.scatter(gw["ABS_POS"], gw["-log10p"],
                           c=GW_RED, s=6, linewidths=0.3,
                           edgecolors="white", zorder=5, rasterized=True)
                # Label top SNP per locus (within 500kb)
                gw_sorted = gw.sort_values("P")
                labelled = []
                for _, row in gw_sorted.iterrows():
                    if all(abs(row["ABS_POS"]-lp) > 5e5 for lp in labelled):
                        lbl = row.get("SNP","")
                        if lbl and lbl != "." and not lbl.startswith("chr"):
                            ax.annotate(lbl,
                                xy=(row["ABS_POS"], row["-log10p"]),
                                xytext=(0, 4), textcoords="offset points",
                                fontsize=4.5, ha="center", va="bottom",
                                color="#1a1a1a",
                                path_effects=[pe.withStroke(linewidth=0.8, foreground="white")])
                        labelled.append(row["ABS_POS"])
                        if len(labelled) >= 8: break

        n_plotted[cohort] = len(data)

        # Thresholds
        ymax = max(data["-log10p"].max() if "-log10p" in data else 8, 8) * 1.12
        ax.axhline(-np.log10(GW),   color=GW_RED,    ls="--", lw=0.7, alpha=0.8)
        ax.axhline(-np.log10(SUGG), color=SUG_AMBER, ls=":",  lw=0.6, alpha=0.7)

        # CHR tick labels
        chrom_mids = data.groupby("CHR")["ABS_POS"].median()
        ax.set_xticks(chrom_mids.values)
        ax.set_xticklabels(
            [str(c) if c not in (6,11,16,21) else "" for c in chrom_mids.index],
            fontsize=5.5)
        ax.set_xlim(data["ABS_POS"].min() - 1e7, data["ABS_POS"].max() + 1e7)
        ax.set_ylim(-0.3, ymax)

        # Y label
        lbl = "−log₁₀(P)"
        if ax_i == 0:
            ax.set_ylabel(f"Cohort A  {lbl}", fontsize=7, labelpad=3)
        else:
            ax.set_ylabel(f"Cohort B  {lbl}", fontsize=7, labelpad=3)
            ax.set_xlabel("Chromosome", fontsize=7)

        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))

        # Cohort label inside panel
        ax.text(0.01, 0.97, cohort.replace("_"," "),
                transform=ax.transAxes, fontsize=7, va="top",
                fontweight="bold", color=COHORT_COL[cohort])

    # Legend
    handles = [
        mpatches.Patch(color=GW_RED,    label=f"GW-sig (p<5×10⁻⁸)"),
        mpatches.Patch(color=SUG_AMBER, label=f"Suggestive (p<10⁻⁵)"),
        mpatches.Patch(color=CHR_COLORS[0], label="Alternating chr"),
    ]
    axes[0].legend(handles=handles, loc="upper right", fontsize=5.5,
                   frameon=True, framealpha=0.9, edgecolor="0.8",
                   ncol=3, columnspacing=0.8)

    # Watermark
    fig.text(0.99, 0.01, f"{author} · {institute}",
             ha="right", va="bottom", fontsize=4.5, color="0.6", style="italic")

    fig.suptitle("Genome-Wide Association Study — Miami Plot",
                 fontsize=9, fontweight="bold", y=1.01)

    save(fig, out_dir, "F1_miami_plot")
    plt.close(fig)
    print(f"  [F1] Done. Cohort A: {n_plotted.get('Cohort_A',0):,}  Cohort B: {n_plotted.get('Cohort_B',0):,} variants plotted")

# ══════════════════════════════════════════════════════════════════
# F2 — QQ Plot grid
# ══════════════════════════════════════════════════════════════════
def figure_qq(sumstat_dir, traits, cohorts, out_dir, author, institute):
    print("[F2] QQ plot grid...")
    ncols = len(traits)
    nrows = len(cohorts)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols*3.5/2.54, nrows*3.5/2.54))
    if nrows == 1: axes = axes[np.newaxis,:]
    if ncols == 1: axes = axes[:,np.newaxis]

    for ri, cohort in enumerate(cohorts):
        for ci, trait in enumerate(traits):
            ax = axes[ri, ci]
            gz = Path(sumstat_dir) / cohort / "sumstats" / f"{cohort}_{trait}_sumstats.tsv.gz"

            if not gz.exists():
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=6)
                ax.set_visible(False)
                continue

            # Load all p-values for QQ (no threshold)
            pvals = []
            with gzip.open(gz, "rt") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    try: pvals.append(float(row["P"]))
                    except: pass
            pvals = np.array([p for p in pvals if 0 < p <= 1])
            if len(pvals) < 10:
                ax.set_visible(False); continue

            n = len(pvals)
            pvals_sorted = np.sort(pvals)

            obs  = -np.log10(pvals_sorted)
            exp  = -np.log10(np.arange(1, n+1) / (n+1))
            lgc  = gc_lambda(pvals)

            # 95% CI band
            ci_lo, ci_hi = qq_ci(n)
            # Subsample for speed (keep all extreme points)
            idx_hi  = np.where(obs > -np.log10(1e-4))[0]
            idx_sub = np.unique(np.concatenate([
                np.round(np.linspace(0, n-1, min(5000,n))).astype(int),
                idx_hi
            ]))
            idx_sub.sort()

            ax.fill_between(exp[idx_sub], ci_lo[idx_sub], ci_hi[idx_sub],
                            color="#d0d0d0", alpha=0.4, linewidth=0, label="95% CI")
            ax.plot([0, exp.max()], [0, exp.max()],
                    color="#666666", lw=0.6, ls="--", zorder=1)

            # Points coloured by significance
            col_arr = np.where(obs[idx_sub] > -np.log10(5e-8),
                               GW_RED,
                               np.where(obs[idx_sub] > -np.log10(1e-5),
                                        SUG_AMBER, COHORT_COL[cohort]))
            ax.scatter(exp[idx_sub], obs[idx_sub],
                       c=col_arr, s=1.5, linewidths=0,
                       rasterized=True, zorder=2)

            # λGC
            ax.text(0.05, 0.95, f"λGC={lgc:.3f}",
                    transform=ax.transAxes, fontsize=5.5,
                    va="top", color="#333333",
                    bbox=dict(facecolor="white", edgecolor="0.8",
                              boxstyle="round,pad=0.2", linewidth=0.4))

            # GW threshold
            ax.axhline(-np.log10(5e-8), color=GW_RED, ls="--", lw=0.5, alpha=0.7)

            lim = max(obs.max(), exp.max()) * 1.08
            ax.set_xlim(-0.2, lim); ax.set_ylim(-0.2, lim)

            # Labels
            if ri == nrows-1:
                ax.set_xlabel("Expected −log₁₀(P)", fontsize=6)
            if ci == 0:
                ax.set_ylabel(f"{cohort.replace('_',' ')}\nObserved", fontsize=6)
            if ri == 0:
                ax.set_title(TRAIT_LABEL.get(trait, trait), fontsize=6.5,
                             fontweight="bold", pad=3,
                             color=TRAIT_COL.get(trait,"#333"))

    fig.text(0.99, 0.01, f"{author} · {institute}",
             ha="right", va="bottom", fontsize=4, color="0.6", style="italic")
    fig.suptitle("QQ Plots — Expected vs Observed −log₁₀(P)",
                 fontsize=9, fontweight="bold", y=1.01)

    save(fig, out_dir, "F2_qq_plots")
    plt.close(fig)
    print("  [F2] Done")

# ══════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════
# F3 — THREE PCA figures:
#   F3a  Combined PCA (all samples) — coloured by fine-grained population
#   F3b  Cohort A PCA (from covariates.tsv) — EUR + EAS populations
#   F3c  Cohort B PCA (from covariates.tsv) — AFR + AMR + SAS populations
# ══════════════════════════════════════════════════════════════════

def _load_eigenvec(path):
    """Load a plink eigenvec file, normalise column names."""
    evc = pd.read_csv(path, sep=r"\s+")
    if "#IID" in evc.columns:
        evc = evc.rename(columns={"#IID": "IID"})
    elif evc.columns[0].startswith("#"):
        evc.columns = ["FID","IID"] + [f"PC{i}" for i in range(1, len(evc.columns)-1)]
    if "FID" not in evc.columns:
        evc.insert(0, "FID", evc["IID"])
    return evc

def _load_covar_pcs(covar_path):
    """Load PCs from a covariates.tsv file (columns: FID IID Age Sex PC1..PCn)."""
    cov = pd.read_csv(covar_path, sep=r"\s+")
    pc_cols = [c for c in cov.columns if c.upper().startswith("PC")]
    if not pc_cols:
        return None
    return cov[["IID"] + pc_cols].rename(
        columns={pc_cols[i]: f"PC{i+1}" for i in range(len(pc_cols))}
    )

def _annotate_populations(evc, population_file):
    """Merge eigenvec with population annotation file."""
    if not population_file or not Path(population_file).exists():
        return evc, None, None
    pop_df = pd.read_csv(population_file, sep="\t")
    pop_df.columns = [c.strip() for c in pop_df.columns]
    id_col   = next((c for c in pop_df.columns if c.strip("#").upper() in ("IID","SAMPLE","ID")), None)
    pop_col  = next((c for c in pop_df.columns if "POP" in c.upper() and "SUPER" not in c.upper()), None)
    spop_col = next((c for c in pop_df.columns if "SUPER" in c.upper() or c.upper()=="SUPERPOP"), None)
    if id_col: pop_df = pop_df.rename(columns={id_col:"IID"})
    if pop_col: pop_df = pop_df.rename(columns={pop_col:"Population"})
    if spop_col: pop_df = pop_df.rename(columns={spop_col:"SuperPop"})
    merge_cols = ["IID"] + (["Population"] if pop_col else []) + (["SuperPop"] if spop_col else [])
    evc = evc.merge(pop_df[merge_cols], on="IID", how="left")
    return evc, "Population" if pop_col else None, "SuperPop" if spop_col else None

def _draw_pca_panel(fig, gs_rows, gs_cols, evc, var_pct,
                    colour_by, col_map, pops, title_prefix, show_legend=True):
    """
    Draw PC1×PC2, PC1×PC3, PC2×PC3 scatter + scree into a GridSpec region.
    gs_rows / gs_cols: index slices into the parent GridSpec.
    """
    axs = [fig.add_subplot(gs_rows[0], gs_cols[i]) for i in range(3)]
    ax_scree = fig.add_subplot(gs_rows[1], gs_cols[:2])
    ax_leg   = fig.add_subplot(gs_rows[1], gs_cols[2])

    def _scatter(ax, pcx, pcy):
        px, py = f"PC{pcx}", f"PC{pcy}"
        if px not in evc.columns or py not in evc.columns:
            ax.text(0.5,0.5,"No data",ha="center",va="center",transform=ax.transAxes)
            return
        for pop in pops:
            sub = evc[evc[colour_by]==pop] if colour_by in evc.columns else evc
            if sub.empty: continue
            ax.scatter(sub[px], sub[py],
                       c=col_map.get(pop,"#999"), s=4, alpha=0.72,
                       linewidths=0, label=pop, rasterized=True)
        vx = var_pct[pcx-1] if pcx-1 < len(var_pct) else 0
        vy = var_pct[pcy-1] if pcy-1 < len(var_pct) else 0
        ax.set_xlabel(f"PC{pcx} ({vx:.1f}%)", fontsize=7)
        ax.set_ylabel(f"PC{pcy} ({vy:.1f}%)", fontsize=7)
        for sp in ["top","right"]: ax.spines[sp].set_visible(False)

    _scatter(axs[0], 1, 2); axs[0].set_title(f"{title_prefix}  PC1 vs PC2", fontsize=7, fontweight="bold")
    _scatter(axs[1], 1, 3); axs[1].set_title(f"{title_prefix}  PC1 vs PC3", fontsize=7, fontweight="bold")
    _scatter(axs[2], 2, 3); axs[2].set_title(f"{title_prefix}  PC2 vs PC3", fontsize=7, fontweight="bold")

    # Scree
    n_show = min(15, len(var_pct))
    ax_scree.bar(range(1,n_show+1), var_pct[:n_show],
                 color="#2256a8", width=0.7, edgecolor="white", linewidth=0.3)
    cum = np.cumsum(var_pct[:n_show])
    ax2 = ax_scree.twinx()
    ax2.plot(range(1,n_show+1), cum, "o-", color="#d62728",
             markersize=2.5, linewidth=0.8)
    ax2.axhline(80, color="#d62728", ls="--", lw=0.5, alpha=0.6)
    ax2.set_ylabel("Cumulative (%)", fontsize=6, color="#d62728")
    ax2.tick_params(axis="y", colors="#d62728", labelsize=5.5)
    ax2.set_ylim(0,105)
    ax_scree.set_xlabel("PC", fontsize=6.5)
    ax_scree.set_ylabel("Variance (%)", fontsize=6.5)
    ax_scree.tick_params(labelsize=5.5)
    for sp in ["top","right"]: ax_scree.spines[sp].set_visible(False)

    # Legend (grouped by SuperPop if available)
    ax_leg.axis("off")
    if show_legend:
        handles = []
        from collections import defaultdict
        if "SuperPop" in evc.columns and colour_by == "Population":
            sp_groups = defaultdict(list)
            for pop in pops:
                rows = evc[evc["Population"]==pop]
                sp = rows["SuperPop"].iloc[0] if len(rows) and not rows["SuperPop"].isna().all() else "?"
                sp_groups[sp].append(pop)
            for sp in sorted(sp_groups.keys()):
                handles.append(mpatches.Patch(
                    color=SUPERPOP_COL.get(sp,"#999"),
                    label=f"▶ {sp}", linewidth=0))
                for pop in sorted(sp_groups[sp]):
                    handles.append(mpatches.Patch(
                        color=POP_COL.get(pop,"#999"),
                        label=f"  {pop}"))
        else:
            handles = [mpatches.Patch(color=col_map.get(p,"#999"), label=p) for p in pops]
        ax_leg.legend(handles=handles, loc="center", fontsize=5,
                      frameon=False, ncol=1 if len(handles)>15 else 1,
                      title="Population", title_fontsize=6,
                      handlelength=1.2, handleheight=1.0)


def figure_pca(combined_eigenvec, combined_eigenval, population_file,
               cohort_a_covar, cohort_b_covar,
               out_dir, author, institute):
    """
    Produces THREE publication figures:
      F3a — Combined PCA: all samples, fine-grained 1000G population colours
      F3b — Cohort A PCA: EUR+EAS samples from covariates.tsv PC columns
      F3c — Cohort B PCA: AFR+AMR+SAS samples from covariates.tsv PC columns
    """
    print("[F3] PCA figures (combined + per-cohort)...")

    # ── Load combined eigenvec + population annotation ─────────────
    evc_all = _load_eigenvec(combined_eigenvec)
    evc_all, pop_col, spop_col = _annotate_populations(evc_all, population_file)

    # Load eigenvalues
    eval_path = str(combined_eigenval)
    eva = pd.read_csv(eval_path, header=None)
    vals = eva[0].values.astype(float)
    var_pct_all = 100 * vals / vals.sum()

    # ── F3a: Combined PCA ──────────────────────────────────────────
    # Always colour by super-population (5 groups: AFR AMR EAS EUR SAS)
    # This is clearer than 26 fine-grained codes for any audience.
    if spop_col and spop_col in evc_all.columns:
        colour_by_all = spop_col
        pops_all = sorted(evc_all[spop_col].dropna().unique())
        col_map_all = {p: SUPERPOP_COL.get(p, "#999") for p in pops_all}
    elif pop_col and pop_col in evc_all.columns:
        # Derive super-pop from fine-grained code via POP_TO_SPOP map
        evc_all["SuperPop_derived"] = evc_all[pop_col].map(POP_TO_SPOP).fillna("Unknown")
        colour_by_all = "SuperPop_derived"
        pops_all = [p for p in ["AFR","AMR","EAS","EUR","SAS","Unknown"]
                    if p in evc_all["SuperPop_derived"].values]
        col_map_all = {p: SUPERPOP_COL.get(p, "#999") for p in pops_all}
    else:
        # Derive cohort from SuperPop if available (EUR+EAS → A, AFR+AMR+SAS → B)
        if "SuperPop" in evc_all.columns:
            spop_to_coh = {"EUR":"Cohort_A","EAS":"Cohort_A",
                           "AFR":"Cohort_B","AMR":"Cohort_B","SAS":"Cohort_B"}
            evc_all["Cohort"] = evc_all["SuperPop"].map(spop_to_coh).fillna("Unknown")
        elif "Population" in evc_all.columns:
            evc_all["Cohort"] = evc_all["Population"].map(POP_TO_SPOP).map(
                {"EUR":"Cohort_A","EAS":"Cohort_A",
                 "AFR":"Cohort_B","AMR":"Cohort_B","SAS":"Cohort_B"}).fillna("Unknown")
        else:
            evc_all["Cohort"] = "Unknown"
        pops_all = ["Cohort_A","Cohort_B"]
        col_map_all = COHORT_COL
        colour_by_all = "Cohort"

    fig_a = plt.figure(figsize=(22/2.54, 16/2.54))
    gs_a  = gridspec.GridSpec(2, 3, figure=fig_a, hspace=0.44, wspace=0.30)

    def scatter_ax(ax, evc, colour_by, col_map, pops, pcx, pcy, var_pct):
        px, py = f"PC{pcx}", f"PC{pcy}"
        if px not in evc.columns or py not in evc.columns: return
        for pop in pops:
            sub = evc[evc[colour_by]==pop] if colour_by in evc.columns else evc
            if sub.empty: continue
            ax.scatter(sub[px], sub[py], c=col_map.get(pop,"#999"),
                       s=3, alpha=0.7, linewidths=0, label=pop, rasterized=True)
        vx = var_pct[pcx-1] if pcx-1<len(var_pct) else 0
        vy = var_pct[pcy-1] if pcy-1<len(var_pct) else 0
        ax.set_xlabel(f"PC{pcx} ({vx:.1f}%)", fontsize=7)
        ax.set_ylabel(f"PC{pcy} ({vy:.1f}%)", fontsize=7)
        for sp in ["top","right"]: ax.spines[sp].set_visible(False)

    def draw_legend(ax, evc, colour_by, col_map, pops):
        """Simple 5-colour super-population legend (AFR/AMR/EAS/EUR/SAS)."""
        ax.axis("off")
        handles = [
            mpatches.Patch(color=col_map.get(p,"#999"), label=p)
            for p in pops if p in col_map
        ]
        ax.legend(handles=handles, loc="center",
                  fontsize=8, frameon=True,
                  framealpha=0.9, edgecolor="#d0cfc6",
                  title="Super-population", title_fontsize=8.5,
                  handlelength=1.4, labelspacing=0.6)

    def draw_scree(ax, var_pct, n=15):
        n_show = min(n, len(var_pct))
        ax.bar(range(1,n_show+1), var_pct[:n_show],
               color="#2256a8", width=0.7, edgecolor="white", linewidth=0.3)
        cum = np.cumsum(var_pct[:n_show])
        ax2 = ax.twinx()
        ax2.plot(range(1,n_show+1), cum, "o-", color="#d62728", markersize=2.5, lw=0.8)
        ax2.axhline(80, color="#d62728", ls="--", lw=0.5, alpha=0.6)
        ax2.set_ylabel("Cumulative (%)", fontsize=6, color="#d62728")
        ax2.tick_params(axis="y", colors="#d62728", labelsize=5.5)
        ax2.set_ylim(0,105)
        ax.set_xlabel("PC", fontsize=6.5); ax.set_ylabel("Variance (%)", fontsize=6.5)
        ax.tick_params(labelsize=5.5)
        for sp in ["top","right"]: ax.spines[sp].set_visible(False)

    # F3a scatter panels
    ax_12 = fig_a.add_subplot(gs_a[0,0]); scatter_ax(ax_12, evc_all, colour_by_all, col_map_all, pops_all, 1, 2, var_pct_all); ax_12.set_title("PC1 vs PC2", fontsize=8, fontweight="bold")
    ax_13 = fig_a.add_subplot(gs_a[0,1]); scatter_ax(ax_13, evc_all, colour_by_all, col_map_all, pops_all, 1, 3, var_pct_all); ax_13.set_title("PC1 vs PC3", fontsize=8, fontweight="bold")
    ax_23 = fig_a.add_subplot(gs_a[0,2]); scatter_ax(ax_23, evc_all, colour_by_all, col_map_all, pops_all, 2, 3, var_pct_all); ax_23.set_title("PC2 vs PC3", fontsize=8, fontweight="bold")
    ax_sc = fig_a.add_subplot(gs_a[1,:2]); draw_scree(ax_sc, var_pct_all)
    ax_lg = fig_a.add_subplot(gs_a[1,2]);  draw_legend(ax_lg, evc_all, colour_by_all, col_map_all, pops_all)

    fig_a.text(0.99, 0.01, f"{author} · {institute}", ha="right", va="bottom",
               fontsize=4, color="0.6", style="italic")
    fig_a.suptitle("Combined PCA — All Samples (EUR+EAS+AFR+AMR+SAS)",
                   fontsize=9, fontweight="bold", y=1.01)
    save(fig_a, out_dir, "F3a_pca_combined")
    plt.close(fig_a)
    print("  [F3a] Combined PCA done")

    # ── F3b: Cohort A (EUR+EAS) from covariates ────────────────────
    for covar_path, cohort_name, cohort_pops, superpops, fig_tag in [
        (cohort_a_covar, "Cohort A (EUR+EAS)",
         ["CEU","TSI","FIN","GBR","IBS","CHB","JPT","CHS","CDX","KHV"],
         ["EUR","EAS"], "F3b_pca_cohort_A"),
        (cohort_b_covar, "Cohort B (AFR+AMR+SAS)",
         ["YRI","LWK","GWD","MSL","ESN","ASW","ACB",
          "MXL","PUR","CLM","PEL",
          "GIH","PJL","BEB","STU","ITU"],
         ["AFR","AMR","SAS"], "F3c_pca_cohort_B"),
    ]:
        if not covar_path or not Path(covar_path).exists():
            print(f"  [WARN] {covar_path} not found — skipping {fig_tag}")
            continue

        # Load PCs from covariates
        cov_pcs = _load_covar_pcs(covar_path)
        if cov_pcs is None:
            print(f"  [WARN] No PC columns in {covar_path}")
            continue

        # Merge with population annotation
        if spop_col and spop_col in evc_all.columns:
            # Merge super-pop directly
            cov_pcs = cov_pcs.merge(
                evc_all[["IID", spop_col]], on="IID", how="left")
            cov_pcs = cov_pcs.rename(columns={spop_col: "SuperPop"})
            present_pops = [p for p in ["AFR","AMR","EAS","EUR","SAS"]
                            if p in cov_pcs["SuperPop"].dropna().values]
            col_map_c = {p: SUPERPOP_COL.get(p, "#999") for p in present_pops}
            colour_by_c = "SuperPop"
        elif pop_col and pop_col in evc_all.columns:
            # Derive super-pop from fine-grained code
            cov_pcs = cov_pcs.merge(
                evc_all[["IID", pop_col]], on="IID", how="left")
            cov_pcs["SuperPop"] = cov_pcs[pop_col].map(POP_TO_SPOP).fillna("Unknown")
            present_pops = [p for p in ["AFR","AMR","EAS","EUR","SAS"]
                            if p in cov_pcs["SuperPop"].dropna().values]
            col_map_c = {p: SUPERPOP_COL.get(p, "#999") for p in present_pops}
            colour_by_c = "SuperPop"
        else:
            present_pops = [cohort_name]
            col_map_c = {cohort_name: COHORT_COL.get(
                "Cohort_A" if "A" in fig_tag else "Cohort_B","#2256a8")}
            colour_by_c = "Cohort"
            cov_pcs["Cohort"] = cohort_name

        if not present_pops:
            print(f"  [WARN] No matching population codes found for {fig_tag}")
            continue

        # Variance from per-cohort covariate PCs (approximate from spread)
        pc_cols_c = [c for c in cov_pcs.columns if c.startswith("PC")]
        var_pct_c = np.array([cov_pcs[c].var() for c in pc_cols_c])
        if var_pct_c.sum() > 0:
            var_pct_c = 100 * var_pct_c / var_pct_c.sum()
        else:
            var_pct_c = np.zeros(len(pc_cols_c))

        fig_c = plt.figure(figsize=(22/2.54, 16/2.54))
        gs_c  = gridspec.GridSpec(2, 3, figure=fig_c, hspace=0.44, wspace=0.30)

        ax_12c = fig_c.add_subplot(gs_c[0,0]); scatter_ax(ax_12c, cov_pcs, colour_by_c, col_map_c, present_pops, 1, 2, var_pct_c); ax_12c.set_title("PC1 vs PC2", fontsize=8, fontweight="bold")
        ax_13c = fig_c.add_subplot(gs_c[0,1]); scatter_ax(ax_13c, cov_pcs, colour_by_c, col_map_c, present_pops, 1, 3, var_pct_c); ax_13c.set_title("PC1 vs PC3", fontsize=8, fontweight="bold")
        ax_23c = fig_c.add_subplot(gs_c[0,2]); scatter_ax(ax_23c, cov_pcs, colour_by_c, col_map_c, present_pops, 2, 3, var_pct_c); ax_23c.set_title("PC2 vs PC3", fontsize=8, fontweight="bold")
        ax_scc = fig_c.add_subplot(gs_c[1,:2]); draw_scree(ax_scc, var_pct_c)
        ax_lgc = fig_c.add_subplot(gs_c[1,2]);  draw_legend(ax_lgc, cov_pcs, colour_by_c, col_map_c, present_pops)

        fig_c.text(0.99, 0.01, f"{author} · {institute}", ha="right", va="bottom",
                   fontsize=4, color="0.6", style="italic")
        fig_c.suptitle(f"PCA — {cohort_name}", fontsize=9, fontweight="bold", y=1.01)
        save(fig_c, out_dir, fig_tag)
        plt.close(fig_c)
        print(f"  [{fig_tag}] {cohort_name} PCA done")

    print("[F3] All PCA figures complete")


# F4 — λGC heatmap + GW-hit summary barplot
# ══════════════════════════════════════════════════════════════════
def figure_lgc_summary(summary_tsvs, out_dir, author, institute):
    print("[F4] λGC heatmap + summary barplot...")

    rows = []
    for f in summary_tsvs:
        if not Path(f).exists(): continue
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows.append(row)
    if not rows:
        print("  [F4] No summary data — skipping")
        return

    df = pd.DataFrame(rows)
    for c in ["lambda_gc","n_gw_sig","n_variants","min_p"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    fig, axes = plt.subplots(1, 2, figsize=(18/2.54, 8/2.54),
                              gridspec_kw={"width_ratios":[1.2,1]})

    # ── λGC heatmap ──
    ax = axes[0]
    traits  = df["trait"].unique().tolist()
    cohorts = df["cohort"].unique().tolist()
    matrix  = np.zeros((len(cohorts), len(traits)))
    for i, coh in enumerate(cohorts):
        for j, tr in enumerate(traits):
            v = df[(df["cohort"]==coh) & (df["trait"]==tr)]["lambda_gc"]
            matrix[i,j] = v.values[0] if len(v) else np.nan

    # Colour: green near 1.0, red above 1.1
    cmap = LinearSegmentedColormap.from_list("lgc",
        ["#1a9641","#a6d96a","#ffffbf","#fdae61","#d73027"], N=256)
    vmin, vmax = 0.9, max(1.3, np.nanmax(matrix))
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax,
                   aspect="auto")

    # Annotate cells
    for i in range(len(cohorts)):
        for j in range(len(traits)):
            v = matrix[i,j]
            if not np.isnan(v):
                fc = "white" if v > 1.15 else "black"
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=6.5, color=fc, fontweight="bold")

    ax.set_xticks(range(len(traits)))
    ax.set_xticklabels([TRAIT_LABEL.get(t,t) for t in traits],
                       rotation=35, ha="right", fontsize=6)
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels([c.replace("_"," ") for c in cohorts], fontsize=7)
    ax.set_title("Genomic Inflation (λGC)", fontsize=7.5, fontweight="bold", pad=8)

    # Threshold line annotation
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.ax.tick_params(labelsize=5.5)
    cbar.ax.axhline(1.0,  color="green", lw=1, ls="--", alpha=0.7)
    cbar.ax.axhline(1.1,  color="red",   lw=1, ls="--", alpha=0.7)
    cbar.set_label("λGC", fontsize=6)

    # ── GW-hit barplot ──
    ax2 = axes[1]
    trait_order = df.groupby("trait")["n_gw_sig"].sum().sort_values(ascending=True).index.tolist()
    x_pos = np.arange(len(trait_order))
    bar_w = 0.35
    for ci, (coh, col) in enumerate(COHORT_COL.items()):
        vals = [df[(df["cohort"]==coh)&(df["trait"]==t)]["n_gw_sig"].sum() for t in trait_order]
        ax2.barh(x_pos + ci*bar_w, vals, height=bar_w,
                 color=col, label=coh.replace("_"," "),
                 edgecolor="white", linewidth=0.3)

    ax2.set_yticks(x_pos + bar_w/2)
    ax2.set_yticklabels([TRAIT_LABEL.get(t,t) for t in trait_order], fontsize=6.5)
    ax2.set_xlabel("GW-significant variants (p < 5×10⁻⁸)", fontsize=6.5)
    ax2.set_title("GW-significant Hits per Trait", fontsize=7.5, fontweight="bold")
    ax2.legend(fontsize=6, frameon=False)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    for sp in ["top","right"]: ax2.spines[sp].set_visible(False)

    fig.text(0.99, 0.01, f"{author} · {institute}",
             ha="right", va="bottom", fontsize=4, color="0.6", style="italic")
    fig.suptitle("Genomic Inflation and GWAS Summary", fontsize=9, fontweight="bold", y=1.02)

    save(fig, out_dir, "F4_lgc_summary")
    plt.close(fig)
    print("  [F4] Done")

# ══════════════════════════════════════════════════════════════════
# F5 — Locus zoom for top hits
# ══════════════════════════════════════════════════════════════════
def compute_ld(ref_bfile, snp_list, lead_snp, plink_bin="plink2"):
    """Compute LD r² between lead SNP and all SNPs in list using plink2."""
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        snp_file = Path(tmp) / "snps.txt"
        snp_file.write_text("\n".join(snp_list))
        cmd = [
            plink_bin, "--bfile", ref_bfile,
            "--extract", str(snp_file),
            "--ld-snp", lead_snp,
            "--r2", "inter-chr",
            "--ld-window-r2", "0",
            "--out", str(Path(tmp) / "ld"),
            "--no-psam-pheno"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        ld_file = Path(tmp) / "ld.vcor2"
        if not ld_file.exists():
            # Try plink 1.9 syntax
            cmd19 = [
                plink_bin.replace("plink2","plink"),
                "--bfile", ref_bfile,
                "--extract", str(snp_file),
                "--ld-snp", lead_snp,
                "--r2", "--ld-window-r2", "0",
                "--ld-window", "1000",
                "--ld-window-kb", "1000",
                "--out", str(Path(tmp) / "ld")
            ]
            subprocess.run(cmd19, capture_output=True, text=True)
            ld_file = Path(tmp) / "ld.ld"

        if not ld_file.exists():
            return {}

        ld_df = pd.read_csv(ld_file, sep=r"\s+")
        snp_col = next((c for c in ld_df.columns if "SNP_B" in c or "ID_B" in c), None)
        r2_col  = next((c for c in ld_df.columns if "R2" in c.upper()), None)
        if snp_col and r2_col:
            return dict(zip(ld_df[snp_col], ld_df[r2_col]))
    return {}

def figure_locus_zoom(sumstat_dir, summary_tsvs, ref_bfile,
                      out_dir, author, institute, n_loci=6):
    print("[F5] Locus zoom...")

    # Find top hits across all cohort×trait combinations
    top_hits = []
    summary_rows = []
    for f in summary_tsvs:
        if not Path(f).exists(): continue
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                summary_rows.append(row)

    for row in summary_rows:
        if int(row.get("n_gw_sig","0")) > 0:
            cohort = row["cohort"]; trait = row["trait"]
            gz = Path(sumstat_dir)/cohort/"sumstats"/f"{cohort}_{trait}_sumstats.tsv.gz"
            if not gz.exists(): continue
            df = load_sumstat(gz, max_pval=5e-8)
            if df.empty: continue
            df = df.sort_values("P")
            top_hits.append({
                "cohort":cohort,"trait":trait,"df":df,
                "lead_snp":df.iloc[0]["SNP"],
                "lead_chr":int(df.iloc[0]["CHR"]),
                "lead_bp":int(df.iloc[0]["BP"]),
                "min_p":float(df.iloc[0]["P"])
            })

    if not top_hits:
        print("  [F5] No GW-significant hits — skipping locus zoom")
        return

    # Sort by p-value, take top N unique loci
    top_hits.sort(key=lambda x: x["min_p"])
    loci = []
    for h in top_hits:
        # Avoid duplicate loci (same chr+bp±2Mb)
        if all(abs(h["lead_chr"]-l["lead_chr"])>0 or abs(h["lead_bp"]-l["lead_bp"])>2e6 for l in loci):
            loci.append(h)
        if len(loci) >= n_loci: break

    if not loci:
        print("  [F5] No unique loci for zoom")
        return

    ncols = min(3, len(loci))
    nrows = (len(loci) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols*6/2.54, nrows*5/2.54),
                              squeeze=False)

    # LD colour map (blue→green→yellow→orange→red)
    ld_cmap = LinearSegmentedColormap.from_list("ld",
        ["#313695","#4575b4","#74add1","#ffffbf","#fdae61","#f46d43","#d73027"], N=256)

    plink_bin = "plink2"
    try:
        subprocess.run([plink_bin,"--version"], capture_output=True)
    except:
        plink_bin = "plink"

    for li, locus in enumerate(loci):
        row_i, col_i = li // ncols, li % ncols
        ax = axes[row_i][col_i]

        chrom = locus["lead_chr"]
        bp    = locus["lead_bp"]
        window = 500_000
        sub = locus["df"].copy()
        sub = sub[(sub["CHR"]==chrom) &
                  (sub["BP"] >= bp-window) &
                  (sub["BP"] <= bp+window)].copy()

        if sub.empty: ax.set_visible(False); continue

        sub["neg_log10p"] = -np.log10(sub["P"].clip(1e-300))
        lead_snp = locus["lead_snp"]

        # Compute LD if bfile available
        ld_r2 = {}
        if ref_bfile and Path(str(ref_bfile)+".bim").exists():
            try:
                ld_r2 = compute_ld(ref_bfile, sub["SNP"].tolist(), lead_snp, plink_bin)
            except Exception as e:
                print(f"    [WARN] LD computation failed: {e}")

        if ld_r2:
            sub["r2"] = sub["SNP"].map(ld_r2).fillna(0.0)
            sub["r2"] = sub["r2"].clip(0, 1)
            sc = ax.scatter(sub["BP"]/1e6, sub["neg_log10p"],
                            c=sub["r2"], cmap=ld_cmap, vmin=0, vmax=1,
                            s=8, linewidths=0.2, edgecolors="0.3",
                            rasterized=True)
            # Colorbar
            cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
            cb.set_label("r²", fontsize=5); cb.ax.tick_params(labelsize=5)
        else:
            ax.scatter(sub["BP"]/1e6, sub["neg_log10p"],
                       c="#2166ac", s=6, linewidths=0,
                       alpha=0.8, rasterized=True)

        # Lead SNP
        lead_row = sub[sub["SNP"]==lead_snp]
        if not lead_row.empty:
            ax.scatter(lead_row["BP"].values/1e6, lead_row["neg_log10p"].values,
                       c="red", s=20, zorder=5, marker="D",
                       edgecolors="white", linewidths=0.5)
            ax.annotate(lead_snp,
                xy=(lead_row["BP"].values[0]/1e6, lead_row["neg_log10p"].values[0]),
                xytext=(0, 5), textcoords="offset points",
                fontsize=4.5, ha="center", color="#c0392b",
                path_effects=[pe.withStroke(linewidth=0.8, foreground="white")])

        ax.axhline(-np.log10(5e-8), color=GW_RED, ls="--", lw=0.6, alpha=0.7)

        title_trait = TRAIT_LABEL.get(locus["trait"], locus["trait"])
        ax.set_title(f"{title_trait} · chr{chrom}:{bp//1e6:.1f} Mb\n{locus['cohort'].replace('_',' ')}",
                     fontsize=6, fontweight="bold", pad=3)
        ax.set_xlabel(f"Chromosome {chrom} position (Mb)", fontsize=6)
        ax.set_ylabel("−log₁₀(P)", fontsize=6)
        for sp in ["top","right"]: ax.spines[sp].set_visible(False)

    # Hide empty axes
    for li in range(len(loci), nrows*ncols):
        axes[li//ncols][li%ncols].set_visible(False)

    fig.text(0.99, 0.01, f"{author} · {institute}",
             ha="right", va="bottom", fontsize=4, color="0.6", style="italic")
    fig.suptitle("Regional Association Plots (Locus Zoom)",
                 fontsize=9, fontweight="bold", y=1.02)

    save(fig, out_dir, "F5_locus_zoom")
    plt.close(fig)
    print(f"  [F5] Done. {len(loci)} loci plotted")

# ══════════════════════════════════════════════════════════════════
# F6 — Forest plot for GW-significant hits
# ══════════════════════════════════════════════════════════════════
def figure_forest(sumstat_dir, summary_tsvs, out_dir, author, institute,
                  max_hits=80):
    print("[F6] Forest plot...")

    rows = []
    for f in summary_tsvs:
        if not Path(f).exists(): continue
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows.append(row)
    if not rows:
        print("  [F6] No summary — skipping"); return

    all_hits = []
    for row in rows:
        if int(row.get("n_gw_sig","0")) == 0: continue
        cohort = row["cohort"]; trait = row["trait"]
        gz = Path(sumstat_dir)/cohort/"sumstats"/f"{cohort}_{trait}_sumstats.tsv.gz"
        if not gz.exists(): continue
        df = load_sumstat(gz, max_pval=5e-8)
        if df.empty: continue
        df["cohort"] = cohort; df["trait"] = trait
        all_hits.append(df)

    if not all_hits:
        print("  [F6] No GW hits — skipping"); return

    hits = pd.concat(all_hits, ignore_index=True)
    hits = hits.sort_values(["trait","BETA"], ascending=[True,False])

    # Cap at max_hits for readability
    if len(hits) > max_hits:
        # Take top hits by p-value per trait
        hits = hits.sort_values("P").groupby("trait").head(max_hits//7).reset_index(drop=True)
        hits = hits.sort_values(["trait","BETA"], ascending=[True,False])

    n = len(hits)
    if n == 0:
        print("  [F6] No hits after filtering — skipping"); return

    fig_h = max(6, n * 0.22 + 1.5)
    fig, ax = plt.subplots(figsize=(14/2.54, min(fig_h, 28)/2.54))

    y_pos = np.arange(n)
    betas = hits["BETA"].values.astype(float)
    ses   = hits["SE"].values.astype(float)
    ci_lo = betas - 1.96 * ses
    ci_hi = betas + 1.96 * ses
    pvals = hits["P"].values.astype(float)

    # Colour by trait
    trait_list = hits["trait"].tolist()
    colours = [TRAIT_COL.get(t,"#333333") for t in trait_list]

    # Error bars
    for yi, (b, lo, hi, col) in enumerate(zip(betas, ci_lo, ci_hi, colours)):
        ax.plot([lo, hi], [yi, yi], color=col, lw=0.9, solid_capstyle="round")
        ax.scatter([b], [yi], color=col, s=10, zorder=3, linewidths=0)

    ax.axvline(0, color="#666", lw=0.6, ls="--")

    # Y tick labels
    ylabels = [f"{row['SNP'][:14]}  {TRAIT_LABEL.get(row['trait'],row['trait'])[:12]}"
               f"  {row['cohort'].replace('Cohort_','C')}"
               for _, row in hits.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ylabels, fontsize=4.5)
    ax.set_xlabel("Effect size (β) ± 95% CI", fontsize=7)
    ax.set_title("GW-Significant Hits — Effect Sizes", fontsize=8, fontweight="bold")
    ax.set_ylim(-0.8, n - 0.2)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)

    # Trait legend
    handles = [mpatches.Patch(color=TRAIT_COL.get(t,"#333"),
                               label=TRAIT_LABEL.get(t,t))
               for t in sorted(set(trait_list))]
    ax.legend(handles=handles, loc="lower right", fontsize=5.5,
              frameon=True, framealpha=0.9, edgecolor="0.8")

    # P-value annotations on right
    ax_r = ax.twinx()
    ax_r.set_ylim(-0.8, n-0.2)
    ax_r.set_yticks(y_pos)
    ax_r.set_yticklabels([f"{p:.1e}" if p < 0.001 else f"{p:.4f}"
                          for p in pvals], fontsize=4)
    ax_r.set_ylabel("P-value", fontsize=6)
    ax_r.spines["top"].set_visible(False)

    fig.text(0.99, 0.01, f"{author} · {institute}",
             ha="right", va="bottom", fontsize=4, color="0.6", style="italic")

    save(fig, out_dir, "F6_forest_plot")
    plt.close(fig)
    print(f"  [F6] Done. {n} hits plotted")

# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(description="Generate publication-quality GWAS figures")
    ap.add_argument("--sumstat_dir",    required=True)
    ap.add_argument("--eigenvec",         required=True,  help="Combined eigenvec (all samples)")
    ap.add_argument("--eigenval",         required=True,  help="Combined eigenval")
    ap.add_argument("--cohort_a_covar",   default="",     help="Cohort A covariates.tsv (contains PC columns)")
    ap.add_argument("--cohort_b_covar",   default="",     help="Cohort B covariates.tsv (contains PC columns)")
    ap.add_argument("--summary_tsvs",     required=True, nargs="+")
    ap.add_argument("--out_dir",          required=True)
    ap.add_argument("--ref_bfile",        default="")
    ap.add_argument("--population_file",  default="")
    ap.add_argument("--cohorts",          default="Cohort_A,Cohort_B")
    ap.add_argument("--quant_traits",     default="ldl_cholesterol,bmi,crp_log,height_cm")
    ap.add_argument("--binary_traits",    default="cad,t2d,hypertension")
    ap.add_argument("--author",           default="Nadeem Khan")
    ap.add_argument("--institute",        default="INRS-CAFSB, Laval, QC, Canada")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    cohorts = [c.strip() for c in args.cohorts.split(",")]
    traits  = ([t.strip() for t in args.quant_traits.split(",")]
               + [t.strip() for t in args.binary_traits.split(",")])

    print(f"\n{'='*60}")
    print(f"  Generating publication figures")
    print(f"  Author    : {args.author}")
    print(f"  Institute : {args.institute}")
    print(f"  Output    : {out}")
    print(f"  Cohorts   : {cohorts}")
    print(f"  Traits    : {traits}")
    print(f"{'='*60}\n")
    t0 = datetime.now()

    figure_miami(args.sumstat_dir, traits, cohorts,
                 out, args.author, args.institute)

    figure_qq(args.sumstat_dir, traits, cohorts,
              out, args.author, args.institute)

    figure_pca(
        combined_eigenvec  = args.eigenvec,
        combined_eigenval  = args.eigenval,
        population_file    = args.population_file or None,
        cohort_a_covar     = args.cohort_a_covar or None,
        cohort_b_covar     = args.cohort_b_covar or None,
        out_dir            = out,
        author             = args.author,
        institute          = args.institute,
    )

    figure_lgc_summary(args.summary_tsvs,
                       out, args.author, args.institute)

    figure_locus_zoom(args.sumstat_dir, args.summary_tsvs,
                      args.ref_bfile or None,
                      out, args.author, args.institute)

    figure_forest(args.sumstat_dir, args.summary_tsvs,
                  out, args.author, args.institute)

    elapsed = (datetime.now()-t0).seconds
    print(f"\n{'='*60}")
    print(f"  All figures complete in {elapsed}s")
    print(f"  Output: {out}/F1_* … F6_*")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
