"""Draw the benchmark figure from benchmark_summary.csv: one row per section, histoalign vs 10x median error.
usage: python plot.py   (writes benchmark_dumbbell.png and .svg next to this file)"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator

HERE = Path(__file__).resolve().parent
NAMES = {
    "breast_rep1": "Breast IDC, rep 1", "breast_rep2": "Breast IDC, rep 2", "liver_normal": "Liver, normal",
    "liver_cancer": "Liver, cancer", "heart": "Heart", "tonsil_flh": "Tonsil, FLH", "tonsil_rfh": "Tonsil, RFH",
    "kidney_normal": "Kidney, normal", "kidney_prcc": "Kidney, papillary RCC", "pancreas_cancer": "Pancreas, PDAC",
    "pancreas_cancer_mm": "Pancreas, PDAC (strips)", "lung_cancer_mm": "Lung, adeno (strips)",
    "colon_normal": "Colon, normal", "colon_cancer": "Colon, adeno", "crc_io": "Colon, adeno (IO)",
    "lung_nsclc_io": "Lung, NSCLC (IO)", "brain_gbm": "Brain, GBM (IO)", "ovary_io": "Ovary, cancer (IO)",
    "kidney_ccrcc": "Kidney, ccRCC (protein)", "ovary_5k": "Ovary, serous (5K)", "lymph_node_5k": "Lymph node (5K)",
    "skin_melanoma_5k": "Skin, melanoma (5K)", "cervix_5k": "Cervix, cancer (5K)",
    "prostate_5k": "Prostate, adeno (5K)", "breast_5k": "Breast, cancer (5K)", "lung_5k": "Lung, adeno (5K)",
    "mouse_brain_5k": "Mouse brain (5K)", "mouse_colon": "Mouse colon", "mouse_femur_formic": "Mouse femur, formic",
    "mouse_femur_edta": "Mouse femur, EDTA", "breast_idc_addon": "Breast IDC, add-on", "breast_idc": "Breast IDC/DCIS",
    "breast_ilc": "Breast ILC"}
BLUE, GREY, LINE, TEXT, TEXT2, GRID = "#2a78d6", "#b9b8b0", "#d9d8d2", "#0b0b0b", "#52514e", "#ecebe7"

df = pd.read_csv(HERE / "benchmark_summary.csv")
df = df[df.status == "ok"].sort_values("tenx_median_um", ascending=True).reset_index(drop=True)
ident = df.baseline.str.startswith("identity").to_numpy()
y = np.arange(len(df))
wins = int((df.ours_median_um < df.tenx_median_um).sum())

fig, ax = plt.subplots(figsize=(7.6, 9.4), dpi=200)
fig.subplots_adjust(left=0.27, right=0.96, top=0.875, bottom=0.095)
for t in (0.2, 0.3, 0.5, 1, 2, 3):
    ax.axvline(t, color=GRID, lw=0.8, zorder=0)
ax.axvline(3, color="#8a8984", lw=0.9, ls=(0, (4, 3)), zorder=0)
ax.hlines(y, df.ours_median_um, df.tenx_median_um, color=LINE, lw=2, zorder=1)
ax.scatter(df.tenx_median_um[~ident], y[~ident], s=42, color=GREY, zorder=2, lw=0)
ax.scatter(df.tenx_median_um[ident], y[ident], s=34, facecolor="white", edgecolor=GREY, lw=1.6, zorder=2)
ax.scatter(df.ours_median_um, y, s=42, color=BLUE, zorder=3, lw=0)
ax.set_xscale("log")
ax.set_xlim(0.15, 3.4)
ax.xaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 1, 2, 3]))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xticklabels(["0.2", "0.3", "0.5", "1", "2", "3"])
ax.set_yticks(y, df.label.map(NAMES).fillna(df.label))
ax.set_ylim(-0.7, len(df) - 0.3)
ax.tick_params(length=0, labelsize=8.5, colors=TEXT2)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xlabel("Median offset, µm (log scale; lower is better)", fontsize=8.5, color=TEXT2)
fig.text(0.96, 0.022, "Dashed line at 3 µm: squares no better than chance pairing are scored as 3 µm", ha="right",
         va="bottom", fontsize=7.5, color=TEXT2)
fig.text(0.02, 0.975, f"histoalign beats 10x's own alignment on {wins} of {len(df)} public Xenium sections",
         fontsize=11.5, weight="bold", color=TEXT, va="top")
fig.text(0.02, 0.948, "Median nucleus offset per 200 µm square, H&E vs Xenium; median across sections "
         f"{df.ours_median_um.median():.2f} µm vs {df.tenx_median_um.median():.2f} µm",
         fontsize=8.5, color=TEXT2, va="top")
handles = [Line2D([], [], marker="o", ls="", ms=6.5, color=BLUE, label="histoalign"),
           Line2D([], [], marker="o", ls="", ms=6.5, color=GREY, label="10x published alignment"),
           Line2D([], [], marker="o", ls="", ms=6, mfc="white", mec=GREY, mew=1.6,
                  label="10x baseline is identity (H&E shipped pre-resampled)")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.93), ncol=3, frameon=False, fontsize=8,
           handletextpad=0.3, columnspacing=1.2, labelcolor=TEXT)
for ext in ("png", "svg"):
    fig.savefig(HERE / f"benchmark_dumbbell.{ext}", facecolor="white")
print("wrote benchmark_dumbbell.png/.svg")
