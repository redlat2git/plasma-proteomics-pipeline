###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 22. Generate the LOCO forest plot
# Requires: Script 17 Source Data
# Produces: publication-ready Figure 4e panel
# Data policy: participant-level inputs and predictions remain local.
###############################################################################

from pathlib import Path
import sys

_REPO_HINT = Path(__file__).resolve()
for _candidate in [_REPO_HINT.parent, *_REPO_HINT.parents]:
    if (_candidate / ".redlat-root").exists() or (_candidate / ".here").exists():
        PROJECT_ROOT = _candidate
        break
else:
    raise FileNotFoundError("Repository root not found. Set REDLAT_PROJECT_ROOT.")
sys.path.insert(0, str(PROJECT_ROOT / "python"))
from redlat_ml.config import load_config, require_files
CONFIG = load_config(__file__)

# ============================================================
# 1. IMPORTS AND USER SETTINGS
# ============================================================

from pathlib import Path
from typing import Optional
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon

# Rendering settings are intentionally identical to Script 17.
FIGURE_WIDTH_CM = 12.0
FIGURE_HEIGHT_CM = 9.0
X_LIMITS = (0.50, 1.00)
PANEL_LETTER = "e"
TEXT_COLOR = "#000000"
REFERENCE_COLOR = "#888888"
POOLED_GUIDE_COLOR = "#B0B0B0"
POOLED_COLOR = "#0072B2"
PREDICTION_COLOR = "#444444"
SHOW_PREDICTION_INTERVAL = True
X_TICKS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
COUNTRY_COLORS = {
    "Argentina": "#795548",
    "Chile": "#8B0000",
    "Colombia": "#56B4E9",
    "Mexico": "#8E7CC3",
    "Peru": "#4D4D4D",
}

PROJECT_DIR_OVERRIDE: Optional[str] = None
PROJECT_DIR = CONFIG.project_root

SOURCE_DIR = (
    CONFIG.publication_root / "loco_meta_analysis" / "source_data"
)
OUTPUT_DIR = (
    CONFIG.publication_root / "figures" / "figure_4" / "loco_forest"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COUNTRY_FILE = SOURCE_DIR / "LOCO_country_AUC_estimates.csv"
SUMMARY_FILE = SOURCE_DIR / "LOCO_AUC_meta_analysis_summary.csv"

for path, label in [
    (COUNTRY_FILE, "Country AUC table"),
    (SUMMARY_FILE, "Meta-analysis summary"),
]:
    if not path.exists():
        raise FileNotFoundError(
            f"{label} not found:\n{path}\n\n"
            "Run scripts/ML/17_ML_meta_analyze_loco.py first."
        )

country = pd.read_csv(COUNTRY_FILE)
summary = pd.read_csv(SUMMARY_FILE)

required_country = {
    "Country", "N_CN", "N_AD", "AUC",
    "CI95_low", "CI95_high", "Weight_random_percent",
}
missing = required_country - set(country.columns)
if missing:
    raise ValueError(f"Missing country columns: {sorted(missing)}")

country["Country"] = pd.Categorical(
    country["Country"],
    categories=["Argentina", "Chile", "Colombia", "Mexico", "Peru"],
    ordered=True,
)
country = country.sort_values("Country").reset_index(drop=True)

random_summary = summary[
    summary["Model"].astype(str).str.startswith("Random")
].iloc[0]

print(country[
    ["Country", "N_CN", "N_AD", "AUC",
     "CI95_low", "CI95_high", "Weight_random_percent"]
].to_string(index=False))

# ============================================================
# 3. DRAW FINAL FOREST PLOT
# ============================================================

def diamond_coordinates(
    center: float,
    lower: float,
    upper: float,
    y: float,
    half_height: float = 0.14,
) -> np.ndarray:
    return np.array([
        [lower, y],
        [center, y + half_height],
        [upper, y],
        [center, y - half_height],
    ])


def auc_ci_text(row: pd.Series) -> str:
    return (
        f"{row['AUC']:.2f} "
        f"[{row['CI95_low']:.2f}, {row['CI95_high']:.2f}]"
    )


fig = plt.figure(
    figsize=(FIGURE_WIDTH_CM / 2.54, FIGURE_HEIGHT_CM / 2.54),
    dpi=300,
)

grid = GridSpec(
    1,
    3,
    figure=fig,
    width_ratios=[1.80, 3.60, 1.55],
    left=0.055,
    right=0.985,
    top=0.78,
    bottom=0.19,
    wspace=0.04,
)

ax_left = fig.add_subplot(grid[0, 0])
ax_forest = fig.add_subplot(grid[0, 1])
ax_right = fig.add_subplot(grid[0, 2])

n_country = len(country)
country_y = np.arange(n_country, 0, -1, dtype=float)
pooled_y = -0.05
prediction_y = -0.82
heterogeneity_y = -1.52

y_min = -1.85
y_max = n_country + 0.75

for axis in [ax_left, ax_forest, ax_right]:
    axis.set_ylim(y_min, y_max)

# ------------------------------------------------------------
# Left table: country and sample sizes
# ------------------------------------------------------------
ax_left.set_xlim(0, 1)
ax_left.axis("off")

header_y = n_country + 0.52

ax_left.text(
    0.02, header_y, "Held-out country",
    fontsize=6.2, fontweight="bold",
    color=TEXT_COLOR, ha="left", va="center",
)
ax_left.text(
    0.98, header_y, "CN / AD",
    fontsize=6.2, fontweight="bold",
    color=TEXT_COLOR, ha="right", va="center",
)

for i, row in country.iterrows():
    y = country_y[i]
    ax_left.text(
        0.02, y, str(row["Country"]),
        fontsize=6.2, color=TEXT_COLOR,
        ha="left", va="center",
    )
    ax_left.text(
        0.98, y,
        f"{int(row['N_CN'])} / {int(row['N_AD'])}",
        fontsize=5.9, color=TEXT_COLOR,
        ha="right", va="center",
    )

ax_left.text(
    0.02, pooled_y, "Random-effects estimate",
    fontsize=6.2, fontweight="bold",
    color=TEXT_COLOR, ha="left", va="center",
)

if SHOW_PREDICTION_INTERVAL:
    ax_left.text(
        0.02, prediction_y, "95% prediction interval",
        fontsize=5.8, color=TEXT_COLOR,
        ha="left", va="center",
    )

# ------------------------------------------------------------
# Forest axis
# ------------------------------------------------------------
ax_forest.set_xlim(*X_LIMITS)
ax_forest.set_yticks([])

ax_forest.xaxis.set_ticks_position("top")
ax_forest.xaxis.set_label_position("top")
ax_forest.set_xticks(X_TICKS)
ax_forest.set_xticklabels(
    [f"{value:.1f}" for value in X_TICKS],
    fontsize=6,
)
ax_forest.set_xlabel("ROC AUC", fontsize=7, labelpad=4)
ax_forest.tick_params(
    axis="x",
    direction="out",
    length=2.5,
    width=0.5,
    pad=2,
)

for side in ["left", "right", "bottom"]:
    ax_forest.spines[side].set_visible(False)
ax_forest.spines["top"].set_linewidth(0.6)

# Chance and pooled guides
ax_forest.axvline(
    0.50,
    linestyle=(0, (2.5, 2.5)),
    linewidth=0.8,
    color=REFERENCE_COLOR,
    zorder=1,
)
ax_forest.axvline(
    random_summary["AUC"],
    linestyle=(0, (1.2, 2.0)),
    linewidth=0.7,
    color=POOLED_GUIDE_COLOR,
    zorder=1,
)

weights = country["Weight_random_percent"].to_numpy(float)
marker_sizes = 4.4 + 3.7 * (weights / weights.max())

for i, row in country.iterrows():
    y = country_y[i]
    color = COUNTRY_COLORS[str(row["Country"])]

    ax_forest.errorbar(
        row["AUC"],
        y,
        xerr=np.array([
            [row["AUC"] - row["CI95_low"]],
            [row["CI95_high"] - row["AUC"]],
        ]),
        fmt="s",
        markersize=marker_sizes[i],
        color=color,
        ecolor=color,
        markerfacecolor=color,
        markeredgecolor=color,
        elinewidth=0.9,
        capsize=2.0,
        capthick=0.7,
        zorder=3,
    )

pooled_diamond = Polygon(
    diamond_coordinates(
        random_summary["AUC"],
        random_summary["CI95_low"],
        random_summary["CI95_high"],
        pooled_y,
    ),
    closed=True,
    facecolor=POOLED_COLOR,
    edgecolor=POOLED_COLOR,
    linewidth=0.6,
    zorder=4,
)
ax_forest.add_patch(pooled_diamond)

if SHOW_PREDICTION_INTERVAL:
    prediction_low = random_summary["Prediction_low"]
    prediction_high = random_summary["Prediction_high"]

    ax_forest.hlines(
        prediction_y,
        prediction_low,
        prediction_high,
        color=PREDICTION_COLOR,
        linewidth=1.0,
        zorder=2,
    )
    ax_forest.vlines(
        [prediction_low, prediction_high],
        prediction_y - 0.10,
        prediction_y + 0.10,
        color=PREDICTION_COLOR,
        linewidth=0.8,
        zorder=2,
    )
    ax_forest.plot(
        random_summary["AUC"],
        prediction_y,
        marker="o",
        markersize=3.2,
        color=PREDICTION_COLOR,
        zorder=3,
    )

# ------------------------------------------------------------
# Right table: numerical estimates
# ------------------------------------------------------------
ax_right.set_xlim(0, 1)
ax_right.axis("off")

ax_right.text(
    0.02, header_y, "AUC [95% CI]",
    fontsize=6.2, fontweight="bold",
    color=TEXT_COLOR, ha="left", va="center",
)

for i, row in country.iterrows():
    ax_right.text(
        0.02, country_y[i], auc_ci_text(row),
        fontsize=5.9, color=TEXT_COLOR,
        ha="left", va="center",
    )

ax_right.text(
    0.02, pooled_y,
    (
        f"{random_summary['AUC']:.2f} "
        f"[{random_summary['CI95_low']:.2f}, "
        f"{random_summary['CI95_high']:.2f}]"
    ),
    fontsize=6.0, fontweight="bold",
    color=TEXT_COLOR, ha="left", va="center",
)

if SHOW_PREDICTION_INTERVAL:
    ax_right.text(
        0.02, prediction_y,
        (
            f"[{random_summary['Prediction_low']:.2f}, "
            f"{random_summary['Prediction_high']:.2f}]"
        ),
        fontsize=5.8, color=TEXT_COLOR,
        ha="left", va="center",
    )

# ------------------------------------------------------------
# Title, panel letter and heterogeneity
# ------------------------------------------------------------
fig.text(
    0.018, 0.935, PANEL_LETTER,
    fontsize=10, fontweight="bold",
    color=TEXT_COLOR, ha="left", va="top",
)
fig.text(
    0.075, 0.935,
    "Meta-analytic synthesis of country-held-out ROC AUCs",
    fontsize=7.8, color=TEXT_COLOR,
    ha="left", va="top",
)

q_p = float(random_summary["Q_P"])
q_p_text = "<0.001" if q_p < 0.001 else f"{q_p:.3f}"

heterogeneity_line_1 = (
    f"Heterogeneity: Q = {random_summary['Q']:.2f}, "
    f"df = {int(random_summary['Q_df'])}, "
    f"P = {q_p_text}"
)
heterogeneity_line_2 = (
    f"I² = {random_summary['I2_percent']:.1f}%; "
    f"τ²$_{{logit-AUC}}$ = {random_summary['Tau2_logit']:.3f}"
)

# Place heterogeneity below the three axes as a two-line annotation.
fig.text(
    0.075, 0.075,
    heterogeneity_line_1,
    fontsize=5.9, color=TEXT_COLOR,
    ha="left", va="bottom",
)
fig.text(
    0.075, 0.045,
    heterogeneity_line_2,
    fontsize=5.9, color=TEXT_COLOR,
    ha="left", va="bottom",
)

pdf_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_Nature_style_v3.pdf"
svg_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_Nature_style_v3.svg"
png_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_Nature_style_v3.png"
tiff_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_Nature_style_v3.tiff"

fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
fig.savefig(png_path, bbox_inches="tight", dpi=600)
fig.savefig(
    tiff_path,
    bbox_inches="tight",
    dpi=600,
    pil_kwargs={"compression": "tiff_lzw"},
)
plt.show()

print("Saved:")
print(pdf_path)
print(svg_path)
print(png_path)
print(tiff_path)


