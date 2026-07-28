"""
Fluorescence Intensity Quantification - ImageJ IntDen-equivalent
Channel-agnostic: handles GFP, RFP, DAPI, CY5, mCherry, or any custom tag.
Groups files by well name and replicate; no fixed channel pairing required.
"""

import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage.filters import threshold_otsu

# ── known channel display styles (extend freely) ─────────────────────────────
_CHANNEL_STYLE = {
    "GFP":     {"cmap": "Greens",  "color": "#2ca02c"},
    "RFP":     {"cmap": "Reds",    "color": "#d62728"},
    "DAPI":    {"cmap": "Blues",   "color": "#1f77b4"},
    "CY5":     {"cmap": "Purples", "color": "#9467bd"},
    "MCHERRY": {"cmap": "Reds",    "color": "#e377c2"},
    "YFP":     {"cmap": "YlOrBr", "color": "#bcbd22"},
    "CFP":     {"cmap": "GnBu",   "color": "#17becf"},
}
_DEFAULT_STYLE = {"cmap": "gray", "color": "#7f7f7f"}

# RGB channel index fallback for multi-channel TIFFs
_CH_INDEX = {"GFP": 1, "RFP": 0, "DAPI": 2}

_INHIBITORS = ["WORTMAN", "UO126", "SB90", "JNK"]

# regex that matches a channel tag anywhere in the stem (case-insensitive)
_CHANNEL_RE = re.compile(
    r"_(" + "|".join([
        "GFP", "RFP", "DAPI", "CY5", "MCHERRY", "YFP", "CFP",
        "TRITC", "FITC", "A488", "A555", "A647", "TXR",
    ]) + r")(?=_\d+$|$)",
    re.IGNORECASE,
)


# ── channel detection ─────────────────────────────────────────────────────────

def detect_channel(stem: str) -> str | None:
    """Return the channel tag found in a file stem, or None."""
    m = _CHANNEL_RE.search(stem)
    return m.group(1).upper() if m else None


def well_base(stem: str) -> str:
    """Strip channel tag (and trailing replicate) to get the well base name."""
    s = _CHANNEL_RE.sub("", stem)          # remove channel tag
    s = re.sub(r"_\d+$", "", s)            # remove replicate suffix
    return s.strip("_")


# ── file grouping ─────────────────────────────────────────────────────────────

def group_files(img_dir: Path) -> dict[str, dict[str, list[Path]]]:
    """
    Returns  { well_base: { CHANNEL: [path, ...], ... }, ... }
    Each well can have multiple channels and multiple replicates per channel.
    Files whose channel cannot be detected are skipped with a warning.
    """
    groups: dict[str, dict[str, list[Path]]] = defaultdict(lambda: defaultdict(list))
    skipped = []
    for p in sorted(img_dir.glob("*.tif")):
        ch = detect_channel(p.stem)
        if ch is None:
            skipped.append(p.name)
            continue
        wb = well_base(p.stem)
        groups[wb][ch].append(p)

    if skipped:
        print(f"Skipped (no channel tag detected): {', '.join(skipped)}")

    return {k: dict(v) for k, v in sorted(groups.items())}


# ── image loading ─────────────────────────────────────────────────────────────

def load_channel(path: Path, channel_key: str) -> np.ndarray:
    img = tifffile.imread(str(path))
    if img.ndim == 2:
        return img.astype(float)
    if img.ndim == 3 and img.shape[2] in (3, 4):   # H×W×C
        idx = _CH_INDEX.get(channel_key, 0)
        idx = min(idx, img.shape[2] - 1)
        return img[:, :, idx].astype(float)
    if img.ndim == 3:                                # C×H×W stack
        idx = _CH_INDEX.get(channel_key, 0)
        idx = min(idx, img.shape[0] - 1)
        return img[idx].astype(float)
    raise ValueError(f"Unexpected image shape {img.shape} in {path}")


# ── metrics ───────────────────────────────────────────────────────────────────

def quantify(channel: np.ndarray) -> dict:
    # 1. Otsu threshold separates foreground (cells) from background
    thresh      = threshold_otsu(channel)
    fg_mask     = channel > thresh                        # True  = foreground pixel
    bg_mask     = ~fg_mask                                # True  = background pixel

    # 2. Estimate background level from pixels BELOW the threshold
    bg_mean     = float(channel[bg_mask].mean()) if bg_mask.any() else 0.0

    # 3. Corrected Total Cell Fluorescence (CTCF) — ImageJ standard:
    #    IntDen(foreground) − background_mean × foreground_area
    fg_pixels   = channel[fg_mask]
    fg_area     = int(fg_mask.sum())
    intden_raw  = float(fg_pixels.sum())
    intden      = intden_raw - bg_mean * fg_area          # ← background subtracted here

    return dict(
        area   = fg_area,
        mean   = float(fg_pixels.mean()) if fg_area else 0.0,
        intden = intden,
        bg_mean= round(bg_mean, 4),
        thresh = round(float(thresh), 4),
    )


# ── per-replicate mask plot ───────────────────────────────────────────────────

def save_mask_plot(well: str, ch_rep_arrays: dict[str, list[np.ndarray]], out_dir: Path) -> None:
    """
    One column per replicate image (GFP_1, GFP_2, RFP_1, RFP_2 …).
    Row 0 = raw image, Row 1 = Otsu foreground mask.
    """
    # build ordered list of (channel, rep_index, array)
    entries = [
        (ch, rep_idx, arr)
        for ch, arrs in sorted(ch_rep_arrays.items())
        for rep_idx, arr in enumerate(arrs)
    ]
    n_cols = len(entries)
    fig, axes = plt.subplots(2, n_cols, figsize=(n_cols * 3.5, 7))
    fig.suptitle(f"{well} — per-replicate masks", fontsize=12)

    if n_cols == 1:
        axes = axes[:, np.newaxis]

    for col, (ch, rep_idx, arr) in enumerate(entries):
        style  = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)
        thresh = threshold_otsu(arr)
        mask   = arr > thresh

        # row 0 – raw image
        axes[0, col].imshow(arr, cmap=style["cmap"], vmin=0, vmax=arr.max() or 1)
        axes[0, col].set_title(f"{ch} rep{rep_idx + 1}\nraw", fontsize=8)
        axes[0, col].axis("off")

        # row 1 – foreground mask
        axes[1, col].imshow(mask, cmap="gray", vmin=0, vmax=1)
        axes[1, col].set_title(f"{ch} rep{rep_idx + 1}\nmask (thresh={thresh:.0f})", fontsize=8)
        axes[1, col].axis("off")

    plt.tight_layout()
    safe = re.sub(r"[^\w]", "_", well)
    fig.savefig(out_dir / f"{safe}_masks.png", dpi=120)
    plt.close(fig)


# ── per-well plot (all channels + histograms + pairwise scatter) ──────────────

def save_well_plot(well: str, channels: dict[str, np.ndarray], out_dir: Path) -> None:
    ch_names = list(channels.keys())
    n = len(ch_names)

    # rows: images | masks (bg subtraction) | histograms | pairwise scatter
    pairs = list(combinations(ch_names, 2))
    n_scatter = len(pairs)
    n_rows = 3 + (1 if n_scatter else 0)   # +1 row for foreground masks
    n_cols = max(n, n_scatter, 1)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(max(12, n_cols * 4), n_rows * 3.5))
    fig.suptitle(well, fontsize=12)

    # ensure axes is always 2-D
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    # row 0 – channel images
    for i, ch in enumerate(ch_names):
        arr = channels[ch]
        style = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)
        axes[0, i].imshow(arr, cmap=style["cmap"], vmin=0, vmax=arr.max() or 1)
        axes[0, i].set_title(f"{ch} channel"); axes[0, i].axis("off")
    for i in range(n, n_cols):
        axes[0, i].axis("off")

    # row 1 – foreground masks (white = kept, black = subtracted as background)
    for i, ch in enumerate(ch_names):
        arr   = channels[ch]
        thresh = threshold_otsu(arr)
        mask  = arr > thresh
        axes[1, i].imshow(mask, cmap="gray", vmin=0, vmax=1)
        axes[1, i].set_title(f"{ch} fg mask (thresh={thresh:.0f})\nwhite=kept  black=subtracted")
        axes[1, i].axis("off")
    for i in range(n, n_cols):
        axes[1, i].axis("off")

    # row 2 – histograms
    for i, ch in enumerate(ch_names):
        arr   = channels[ch]
        thresh = threshold_otsu(arr)
        style = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)
        axes[2, i].hist(arr.ravel(), bins=64, color=style["color"], alpha=0.75)
        axes[2, i].axvline(thresh, color="red", linewidth=1.2, linestyle="--",
                           label=f"Otsu={thresh:.0f}")
        axes[2, i].legend(fontsize=7)
        axes[2, i].set_title(f"{ch} distribution")
        axes[2, i].set_xlabel("Intensity"); axes[2, i].set_ylabel("Count")
    for i in range(n, n_cols):
        axes[2, i].axis("off")

    # row 3 – pairwise scatter
    if n_scatter:
        for i, (a, b) in enumerate(pairs):
            av, bv = channels[a].ravel(), channels[b].ravel()
            idx = np.random.choice(len(av), min(5000, len(av)), replace=False)
            axes[3, i].scatter(av[idx], bv[idx], s=1, alpha=0.3, color="purple")
            axes[3, i].set_title(f"{a} vs {b}")
            axes[3, i].set_xlabel(a); axes[3, i].set_ylabel(b)
        for i in range(n_scatter, n_cols):
            axes[3, i].axis("off")

    plt.tight_layout()
    safe = re.sub(r"[^\w]", "_", well)
    fig.savefig(out_dir / f"{safe}_analysis.png", dpi=120)
    plt.close(fig)


# ── summary plots (per-channel IntDen bars) ───────────────────────────────────

def save_summary_plots(df: pd.DataFrame, out_dir: Path, channels: list[str]) -> None:
    intden_cols = [f"{ch}_intden" for ch in channels if f"{ch}_intden" in df.columns]
    if not intden_cols:
        return

    wells = df["well"]
    x = np.arange(len(wells))
    w = 0.8 / len(intden_cols)

    fig, ax = plt.subplots(figsize=(max(8, len(wells) * 0.7), 5))
    for i, col in enumerate(intden_cols):
        ch = col.replace("_intden", "")
        color = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)["color"]
        ax.bar(x + i * w - (len(intden_cols) - 1) * w / 2,
               df[col], w, label=ch, color=color, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(wells, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("IntDen"); ax.set_title("Integrated Density per Well")
    ax.legend(); plt.tight_layout()
    fig.savefig(out_dir / "summary_intden.png", dpi=120); plt.close(fig)

    # one ratio plot per detected pair
    ratio_cols = [c for c in df.columns if c.endswith("_ratio")]
    for rc in ratio_cols:
        fig, ax = plt.subplots(figsize=(max(8, len(wells) * 0.7), 5))
        ax.bar(x, df[rc], color="teal", alpha=0.8)
        ax.axhline(1, color="k", linestyle="--", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(wells, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(rc.replace("_", " ")); ax.set_title(rc.replace("_", " ") + " per Well")
        plt.tight_layout()
        fig.savefig(out_dir / f"summary_{rc}.png", dpi=120); plt.close(fig)


# ── condition parsing ─────────────────────────────────────────────────────────

def parse_condition(well: str) -> tuple[str, str]:
    stem = re.sub(r"_\d+$", "", well)
    inhibitor = "Control"
    for inh in _INHIBITORS:
        if re.search(rf"_{inh}$", stem, re.IGNORECASE):
            inhibitor = inh
            stem = re.sub(rf"_{inh}", "", stem, flags=re.IGNORECASE)
            break
    return stem, inhibitor


# ── condition-level plots ─────────────────────────────────────────────────────

def _bar_with_replicates(ax, cond_df, col, color, ylabel, title):
    grp   = cond_df.groupby("condition")[col]
    means = grp.mean(); sds = grp.std(ddof=1).fillna(0)
    conds = means.index.tolist(); x = np.arange(len(conds))
    ax.bar(x, means, color=color, alpha=0.75, zorder=2)
    ax.errorbar(x, means, yerr=sds, fmt="none", color="black",
                capsize=4, linewidth=1.2, zorder=3)
    for i, cond in enumerate(conds):
        vals = cond_df.loc[cond_df["condition"] == cond, col].dropna()
        ax.scatter(np.full(len(vals), i), vals, color="black", s=25, zorder=4, alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(conds, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(ylabel); ax.set_title(title)


def save_condition_plots(df: pd.DataFrame, out_dir: Path, channels: list[str]) -> None:
    parsed = df["well"].apply(parse_condition)
    df = df.copy()
    df["treatment_group"] = [p[0] for p in parsed]
    df["inhibitor"]       = [p[1] for p in parsed]
    df["condition"]       = df["treatment_group"] + " / " + df["inhibitor"]

    intden_cols = [f"{ch}_intden" for ch in channels if f"{ch}_intden" in df.columns]
    ratio_cols  = [c for c in df.columns if c.endswith("_ratio")]
    metric_cols = intden_cols + ratio_cols

    if not metric_cols:
        return

    # per-channel IntDen by condition
    for col in intden_cols:
        ch = col.replace("_intden", "")
        color = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)["color"]
        fig, ax = plt.subplots(figsize=(max(8, df["condition"].nunique() * 0.9), 5))
        _bar_with_replicates(ax, df, col, color, "IntDen", f"{ch} IntDen by Condition")
        plt.tight_layout()
        fig.savefig(out_dir / f"condition_{col}.png", dpi=120); plt.close(fig)

    # per-ratio by condition
    for col in ratio_cols:
        fig, ax = plt.subplots(figsize=(max(8, df["condition"].nunique() * 0.9), 5))
        _bar_with_replicates(ax, df, col, "#17becf",
                             col.replace("_", " "), col.replace("_", " ") + " by Condition")
        ax.axhline(1, color="k", linestyle="--", linewidth=0.8)
        plt.tight_layout()
        fig.savefig(out_dir / f"condition_{col}.png", dpi=120); plt.close(fig)

    # inhibitor comparison
    for tg, grp_df in df.groupby("treatment_group"):
        if grp_df["inhibitor"].nunique() < 2:
            continue
        plot_cols = (intden_cols + ratio_cols)[:4]   # cap at 4 subplots
        fig, axes = plt.subplots(1, len(plot_cols), figsize=(4 * len(plot_cols), 5))
        if len(plot_cols) == 1:
            axes = [axes]
        fig.suptitle(f"Inhibitor comparison - {tg}", fontsize=12)
        for ax, col in zip(axes, plot_cols):
            ch = col.replace("_intden", "").replace("_ratio", "")
            color = _CHANNEL_STYLE.get(ch, _DEFAULT_STYLE)["color"]
            inh_df = grp_df.copy(); inh_df["condition"] = inh_df["inhibitor"]
            _bar_with_replicates(ax, inh_df, col, color,
                                 col.replace("_", " "), col.replace("_", " "))
        plt.tight_layout()
        safe = re.sub(r"[^\w]", "_", tg)
        fig.savefig(out_dir / f"inhibitor_comparison_{safe}.png", dpi=120); plt.close(fig)

    # heatmap (z-scored)
    heat = (df.groupby("condition")[metric_cols].mean()
              .apply(lambda c: (c - c.mean()) / c.std(ddof=1) if c.std(ddof=1) else c))
    fig, ax = plt.subplots(figsize=(len(metric_cols) * 1.4, max(4, len(heat) * 0.55)))
    im = ax.imshow(heat.values, aspect="auto", cmap="RdYlGn")
    ax.set_xticks(range(len(metric_cols)))
    ax.set_xticklabels([m.replace("_", "\n") for m in metric_cols], fontsize=8)
    ax.set_yticks(range(len(heat))); ax.set_yticklabels(heat.index, fontsize=8)
    plt.colorbar(im, ax=ax, label="z-score")
    ax.set_title("Condition × Metric Heatmap (z-scored)")
    plt.tight_layout()
    fig.savefig(out_dir / "condition_heatmap.png", dpi=120); plt.close(fig)

    # replicate consistency (first intden column)
    if intden_cols:
        first_col = intden_cols[0]
        df["rep"] = df["well"].str.extract(r"_(\d+)$")[0].fillna("1")
        pivot = df.pivot_table(index="condition", columns="rep",
                               values=first_col, aggfunc="first")
        if pivot.shape[1] >= 2:
            cols = pivot.columns.tolist()
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(pivot[cols[0]], pivot[cols[1]], color="green", alpha=0.8, s=50)
            for cond, row in pivot.iterrows():
                ax.annotate(cond, (row[cols[0]], row[cols[1]]),
                            fontsize=6, ha="left", va="bottom")
            lim = max(pivot.max().max() * 1.05, 1)
            ax.plot([0, lim], [0, lim], "k--", linewidth=0.8)
            ax.set_xlabel(f"Rep {cols[0]} {first_col}")
            ax.set_ylabel(f"Rep {cols[1]} {first_col}")
            ax.set_title(f"Replicate Consistency ({first_col})")
            plt.tight_layout()
            fig.savefig(out_dir / "replicate_consistency.png", dpi=120); plt.close(fig)

    # condition summary CSV
    cond_summary = (df.groupby("condition")[metric_cols]
                      .agg(["mean", "std", "count"]).round(4))
    cond_summary.columns = ["_".join(c) for c in cond_summary.columns]
    cond_summary.to_csv(out_dir / "condition_summary.csv")
    print(f"Condition summary saved -> {out_dir / 'condition_summary.csv'}")


# ── main ──────────────────────────────────────────────────────────────────────

def main(img_dir: str = "images") -> None:
    img_dir = Path(img_dir)
    out_dir = img_dir.parent / "results"
    out_dir.mkdir(exist_ok=True)

    groups = group_files(img_dir)
    if not groups:
        sys.exit("No files with recognised channel tags found. Check file naming.")

    # collect all channel names seen across the whole experiment
    all_channels: list[str] = sorted({ch for chs in groups.values() for ch in chs})
    print(f"Detected channels: {', '.join(all_channels)}")
    print(f"Found {len(groups)} wells -> results in {out_dir}\n")

    records = []
    for well, ch_files in groups.items():
        print(f"  Processing: {well}  channels: {', '.join(sorted(ch_files))}")

        # load all replicate arrays per channel
        ch_rep_arrays: dict[str, list[np.ndarray]] = {}
        ch_paths:      dict[str, str]              = {}
        for ch, paths in ch_files.items():
            ch_rep_arrays[ch] = [load_channel(p, ch) for p in paths]
            ch_paths[ch]      = "|".join(p.name for p in paths)

        # averaged arrays used only for the well plot
        ch_arrays_avg: dict[str, np.ndarray] = {
            ch: (np.mean(arrs, axis=0) if len(arrs) > 1 else arrs[0])
            for ch, arrs in ch_rep_arrays.items()
        }

        # one row per replicate so std across replicates is meaningful
        n_reps = max(len(v) for v in ch_rep_arrays.values())
        for rep_idx in range(n_reps):
            row: dict = {"well": well, "replicate": rep_idx + 1}
            for ch in all_channels:
                if ch in ch_rep_arrays:
                    arrs = ch_rep_arrays[ch]
                    arr  = arrs[rep_idx] if rep_idx < len(arrs) else arrs[-1]
                    m    = quantify(arr)
                    row[f"{ch}_file"]    = ch_paths[ch]
                    row[f"{ch}_mean"]    = round(m["mean"],   4)
                    row[f"{ch}_intden"]  = round(m["intden"], 2)
                    row[f"{ch}_bg_mean"] = m["bg_mean"]
                    row[f"{ch}_thresh"]  = m["thresh"]
                    if "image_area" not in row:
                        row["image_area"] = m["area"]
                else:
                    row[f"{ch}_file"]   = ""
                    row[f"{ch}_mean"]   = np.nan
                    row[f"{ch}_intden"] = np.nan

            # pairwise ratios
            present = [ch for ch in all_channels if ch in ch_rep_arrays]
            for a, b in combinations(present, 2):
                denom = row[f"{b}_intden"]
                row[f"{a}_{b}_ratio"] = (
                    round(row[f"{a}_intden"] / denom, 6) if denom else np.nan
                )

            records.append(row)

        save_well_plot(well, ch_arrays_avg, out_dir)
        save_mask_plot(well, ch_rep_arrays, out_dir)

    df = pd.DataFrame(records)
    csv_path = out_dir / "fluorescence_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved -> {csv_path}")

    save_summary_plots(df, out_dir, all_channels)
    save_condition_plots(df, out_dir, all_channels)
    print(f"Summary + condition plots saved -> {out_dir}")
    print("\n" + df.to_string(index=False))


if __name__ == "__main__":
    img_dir = sys.argv[1] if len(sys.argv) > 1 else "images"
    main(img_dir)
