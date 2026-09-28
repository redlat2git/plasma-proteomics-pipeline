###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 16. Run matched analyses and Extended Data Figure 10
# Requires: Script 15 outputs and validated primary-model results
# Produces: matched sensitivity results, figure files and audited Source Data
# Data policy: participant-level inputs and predictions remain local.
###############################################################################

from pathlib import Path
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

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

# Extended Data Figure 10 — matched-cohort sensitivity analysis
#
# This script executes the propensity-score matched nested classification workflow,
# evaluating both a newly reselected panel in the balanced cohort and the primary
# seven-protein panel evaluated without feature reselection, alongside p-tau217
# benchmarking and multivariable clinical regressions.

# ============================================================
# 0. DEPENDENCIES
# ============================================================
# Runtime dependencies are declared in environment/requirements.txt and
# environment/environment.yml. This analysis script does not install packages.

# ## Analytical distinction
#
# The matched nested analysis is allowed to select different proteins because the
# complete filtering and elastic-net pipeline is rerun in the demographically
# balanced sample. A complementary fixed-panel analysis evaluates the exact seven
# proteins from the main analysis. Neither analysis is interpreted as independent
# external validation.

# ============================================================
# 1. IMPORTS AND CONFIGURATION
# ============================================================

print("=" * 72)
print("EXTENDED DATA FIG. 10 — MATCHED COHORT ANALYSIS")
print("matched source policy: audited matched_ids_SELECTED.csv only")
print("figure policy: exact Fig. 4 page size; compact three-row layout")
print("=" * 72)


import json
import math
import os
import re
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from scipy import stats
from scipy.special import expit
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import (
    accuracy_score,
    auc as sklearn_auc,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVC
import statsmodels.api as sm

warnings.filterwarnings("ignore")

try:
    from skopt import BayesSearchCV
    from skopt.space import Categorical, Real
    HAS_SKOPT = True
except Exception:
    HAS_SKOPT = False

# ------------------------------------------------------------
# EDIT ONLY THIS PATH
# ------------------------------------------------------------
PROJECT_DIR = CONFIG.project_root
DATA_DIR = CONFIG.data_dir
MASTER_FILE = CONFIG.master_file
SELECTED_MATCH_IDS_FILE = CONFIG.private_root / "Result_matching_rebuild_v5" / "matched_ids_SELECTED.csv"
SELECTED_MATCH_FULL_FILE = CONFIG.private_root / "Result_matching_rebuild_v5" / "Matched_Output_SELECTED.csv"
HISTORICAL_MATCH_FILE_OVERRIDE = None
HISTORICAL_MATCHING_INPUT_OVERRIDE = None
RSCRIPT_OVERRIDE = None
RAW_METADATA_OVERRIDE = CONFIG.metadata_file
REVIEWER_PROTEOMICS_FILE = CONFIG.proteomics_file
ARCHIVE_SEARCH_MAX_DEPTH = 2
HISTORICAL_INPUT_N = 639
HISTORICAL_INPUT_CN = 313
HISTORICAL_INPUT_AD = 326
APOE_DIR = CONFIG.private_root / "apoe" / "02_nested_cv"
FIXED_PTAU_DIR = CONFIG.private_root / "fixed_panel_ptau"
MATCH_ROOT = CONFIG.private_root / "matched"
MATCH_DIR = MATCH_ROOT / "00_matching"
MATCH_RESULTS_DIR = MATCH_ROOT / "02_nested_cv"
MATCH_FINAL_PANEL_DIR = MATCH_RESULTS_DIR / "Final_Panels"
PRIMARY_FIXED_MATCH_DIR = MATCH_ROOT / "03_primary_fixed_panel_cv"
MATCH_PTAU_DIR = MATCH_ROOT / "04_primary_panel_ptau"
MATCH_REG_DIR = MATCH_ROOT / "05_clinical_regression"
OUTPUT_ROOT = CONFIG.publication_root / "extended_data_figure_10"
FIGURE_DIR = OUTPUT_ROOT / "figures"
SOURCE_DATA_DIR = OUTPUT_ROOT / "source_data"
OUTPUT_LOCATION_FILE = OUTPUT_ROOT / "OUTPUT_LOCATION.txt"

# Reviewer-project output directories
for _dir in [
    MATCH_DIR,
    MATCH_RESULTS_DIR,
    MATCH_FINAL_PANEL_DIR,
    PRIMARY_FIXED_MATCH_DIR,
    MATCH_PTAU_DIR,
    MATCH_REG_DIR,
    FIGURE_DIR,
    SOURCE_DATA_DIR,
]:
    _dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Reproducibility and analysis settings
# ------------------------------------------------------------
RANDOM_STATE = 42
MATCHING_SEED = 1111
np.random.seed(RANDOM_STATE)

OUTER_FOLDS = 5
INNER_FOLDS = 5
DEP_FDR = 0.05
MIN_DEP_GENES = 20
FALLBACK_TOP_N_DEP = 100
INNER_STABILITY = 0.80
FINAL_PANEL_THRESHOLD = 0.80
# The recurrent matched panel is defined strictly by the prespecified
# outer-fold selection-frequency threshold. It is never padded to a fixed size.
MATCHED_PANEL_DISPLAY_MAX: Optional[int] = None

MAX_ITER = 50000
TOL = 1e-3
SVM_SEARCH_ITER = 50
N_JOBS = -1

# Matching specification reported in the manuscript.
MATCH_RATIO = 1
MATCH_CALIPER_SD = 0.20
MATCH_REPLACE = False
EXACT_COUNTRY = False

# Prefer a historical matched-ID file with exactly 191 CN + 191 AD when found.
TARGET_CN = 191
TARGET_AD = 191
STRICT_MANUSCRIPT_COUNTS = True

# Set True only when you intentionally want to rerun expensive matched models.
FORCE_RECOMPUTE_MATCHED = False
FORCE_RECOMPUTE_FIXED_PANEL = False
FORCE_RECOMPUTE_MATCHED_PTAU = False
FORCE_RECOMPUTE_REGRESSION = False

# Optional manual override. Leave as None to derive the matched panel from
# cross-validation selection frequency.
MATCHED_PANEL_OVERRIDE: Optional[List[str]] = None

PRIMARY_PANEL = [
    "SPC25", "CPLX2", "TCP11L1", "ACHE",
    "ODC1", "SPON1", "RTN4RL1",
]

PTAU_CANDIDATES = ["p.tau217", "p-tau217", "pTau217", "ptau217"]

META_COLS = [
    "SampleGroup", "Site", "Country", "Sex", "Age", "Education",
    "ApoE", "APOE_group", "APOE4_carrier",
    "cdr_global", "cdr_boxscore", "mmse_total", "udsfaq_total",
    "cog_benson", "cog_tmt_a", "cog_tmt_b", "cog_craft_verb_delayed",
    "NPI", "Mini.SEA", "T.ADLQ",
    "p.tau217", "p-tau217", "pTau217", "ptau217",
    "p.tau181", "NfL", "ratio.AB42.40", "GFAP_1",
    "distance", "weights", "subclass", "SampleGroup_bin",
    "pair_id", "propensity_score", "propensity_logit",
]

# Colour-blind-compatible palette consistent with the manuscript figures.
COL_PROTEIN = "#0072B2"
COL_FIXED_PANEL = "#7A5195"
COL_APOE = "#E69F00"
COL_PTAU = "#D55E00"
COL_COMBINED = "#009E73"
COL_HIGHER_AD = "#C84E4B"
COL_LOWER_AD = "#1F9AC9"
COL_HIGHER_AD_LIGHT = "#E8BDBB"
COL_LOWER_AD_LIGHT = "#B7D5E6"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

print("PROJECT_DIR:", PROJECT_DIR)
print("HAS_SKOPT:", HAS_SKOPT)
print("EXACT_COUNTRY:", EXACT_COUNTRY)

# ============================================================
# 2. GENERIC HELPERS AND PREFLIGHT
# ============================================================

def require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found:\n{path}")
    return path


def first_existing(paths: Sequence[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def find_column(df: pd.DataFrame, candidates: Sequence[str]) -> str:
    direct = {str(c): str(c) for c in df.columns}
    lower = {str(c).lower(): str(c) for c in df.columns}
    for candidate in candidates:
        if candidate in direct:
            return direct[candidate]
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"None of these columns were found: {list(candidates)}")


def clean_sex(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.strip()
    return values.replace({
        "1": "F", "Female": "F", "female": "F",
        "Mujer": "F", "mujer": "F",
        "2": "M", "Male": "M", "male": "M",
        "Hombre": "M", "hombre": "M",
    })


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    q_values = np.full_like(p_values, np.nan, dtype=float)
    valid = np.isfinite(p_values)
    p = p_values[valid]
    if p.size == 0:
        return q_values
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q_values[valid] = restored
    return q_values


def safe_log2_matrix(values: np.ndarray) -> np.ndarray:
    array = values.astype(float, copy=True)
    array[array <= 0] = np.nan
    return np.log2(array)


def standard_error_text(mean: float, sd: float, digits: int = 2) -> str:
    return f"{mean:.{digits}f} ± {sd:.{digits}f}"


def format_p(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "NA"
    if p_value < 1e-4:
        return f"{p_value:.1e}"
    if p_value < 0.001:
        return f"{p_value:.3f}"
    return f"{p_value:.3f}"


def preflight() -> pd.DataFrame:
    checks = [
        ("Master matrix", MASTER_FILE),
        ("APOE ROC", APOE_DIR / "roc_curves.csv"),
        ("APOE permutation", APOE_DIR / "Final_Panels" / "permutation_importance_freq08.csv"),
        ("APOE DeLong", APOE_DIR / "delong_results.csv"),
        ("Fixed-panel p-tau ROC", FIXED_PTAU_DIR / "roc_fold_curves.csv"),
        ("Fixed-panel p-tau DeLong", FIXED_PTAU_DIR / "delong_results.csv"),
        ("Fixed-panel model statistics", FIXED_PTAU_DIR / "model_statistics.csv"),
    ]
    report = pd.DataFrame({
        "Input": [label for label, _ in checks],
        "Path": [str(path) for _, path in checks],
        "Exists": [path.exists() for _, path in checks],
    })
    print(report.to_string(index=False))
    missing = report.loc[~report["Exists"]]
    if not missing.empty:
        raise FileNotFoundError(
            "Required inputs are missing. Restore the folders shown above before continuing."
        )
    return report


preflight_report = preflight()

# ============================================================
# 3. LOAD MASTER MATRIX AND IDENTIFY PROTEINS
# ============================================================

master = pd.read_csv(require_file(MASTER_FILE, "Master matrix"), low_memory=False)

required_meta = ["SampleId", "SampleGroup", "Sex", "Age", "Education", "Country"]
missing_meta = [column for column in required_meta if column not in master.columns]
if missing_meta:
    raise ValueError(f"Master matrix is missing: {missing_meta}")

master["SampleId"] = master["SampleId"].astype(str)
master = master[master["SampleGroup"].isin(["CN", "AD"])].copy()
master["Sex"] = clean_sex(master["Sex"])
master["Age"] = pd.to_numeric(master["Age"], errors="coerce")
master["Education"] = pd.to_numeric(master["Education"], errors="coerce")
master["Country"] = master["Country"].astype(str).str.strip()

ptau_col = find_column(master, PTAU_CANDIDATES)

candidate_protein_cols = [
    column for column in master.columns
    if column not in META_COLS and column != "SampleId"
]

# Numeric conversion and basic quality filter.
protein_numeric = master[candidate_protein_cols].apply(pd.to_numeric, errors="coerce")
nonmissing_rate = protein_numeric.notna().mean(axis=0)
variance = protein_numeric.var(axis=0, skipna=True)

protein_cols = [
    column for column in candidate_protein_cols
    if nonmissing_rate.get(column, 0) >= 0.80
    and np.isfinite(variance.get(column, np.nan))
    and variance.get(column, 0) > 0
]

print("CN/AD participants:", len(master))
print("p-tau217 column:", ptau_col)
print("Eligible proteins:", len(protein_cols))


# ## Matching source policy
#
# Matching is treated as closed. The workflow accepts only the selected output
# from `Result_matching_rebuild_v5` and stops if it is not exactly 191 CN +
# 191 AD.

# ## 4. Audited selected matched cohort
#
# Reads only the `matched_ids_SELECTED.csv` produced by the authoritative R
# matching audit (Script 15). It does not rerun propensity-score matching and
# does not search for alternative matched files.

# ============================================================
# 4. AUTHORITATIVE SELECTED MATCHED COHORT
# ============================================================
# Matching has already been completed and audited in R. This script must
# never rerun MatchIt or choose among alternative matched files.

def _find_header_column(
    columns: Sequence[str],
    candidates: Sequence[str],
) -> Optional[str]:
    direct = {str(column): str(column) for column in columns}
    lower = {
        str(column).strip().lower(): str(column)
        for column in columns
    }
    for candidate in candidates:
        if candidate in direct:
            return direct[candidate]
        key = str(candidate).strip().lower()
        if key in lower:
            return lower[key]
    return None


def discover_named_files(
    names: Sequence[str],
) -> List[Path]:
    """Search only within the current reviewer reproducibility project."""
    # Reviewer reproducibility must be self-contained. Never search historical
    # projects or the user's home directory for alternative inputs.
    roots = [
        PROJECT_DIR,
        CONFIG.data_dir,
        CONFIG.project_root,
        CONFIG.private_root,
    ]
    wanted = {str(name).lower() for name in names}
    found: List[Path] = []

    for root in roots:
        if not root.exists():
            continue

        # Fast direct checks.
        for name in names:
            direct_candidates = [
                root / name,
                root / "Data" / name,
                CONFIG.data_dir / name,
                root / "03_ML" / name,
            ]
            for candidate in direct_candidates:
                if candidate.exists():
                    found.append(candidate)

        # Recursive search with bounded depth to avoid scanning all of Desktop.
        root_depth = len(root.parts)
        try:
            for candidate in root.rglob("*"):
                if (
                    candidate.is_file()
                    and candidate.name.lower() in wanted
                    and len(candidate.parts) - root_depth <= ARCHIVE_SEARCH_MAX_DEPTH
                ):
                    found.append(candidate)
        except (PermissionError, OSError):
            pass

    unique: List[Path] = []
    seen = set()
    for path in found:
        key = str(path.resolve()).lower()
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def discover_matching_csvs() -> List[Path]:
    return discover_named_files([
        "matched_ids_SELECTED.csv",
        "Matched_Output_SELECTED.csv",
        "Matched_Output_v2.csv",
        "Matched_Output_MatchIt.csv",
    ])


def matching_balance_table(
    frame: pd.DataFrame,
    stage: str,
) -> pd.DataFrame:
    rows = []
    for group in ["CN", "AD"]:
        subset = frame[frame["SampleGroup"] == group]
        rows.append({
            "Stage": stage,
            "SampleGroup": group,
            "N": len(subset),
            "Age_mean": pd.to_numeric(
                subset["Age"], errors="coerce"
            ).mean(),
            "Age_sd": pd.to_numeric(
                subset["Age"], errors="coerce"
            ).std(ddof=1),
            "Education_mean": pd.to_numeric(
                subset["Education"], errors="coerce"
            ).mean(),
            "Education_sd": pd.to_numeric(
                subset["Education"], errors="coerce"
            ).std(ddof=1),
            "Female_n": int((subset["Sex"] == "F").sum()),
            "Female_percent": 100 * float(
                (subset["Sex"] == "F").mean()
            ),
        })
    return pd.DataFrame(rows)


require_file(
    SELECTED_MATCH_IDS_FILE,
    "Audited selected matched IDs. Run scripts/ML/15_ML_rebuild_matching.R first",
)

selected_ids_raw = pd.read_csv(
    SELECTED_MATCH_IDS_FILE,
    low_memory=False,
)
if "SampleId" not in selected_ids_raw.columns:
    raise ValueError(
        f"SampleId is absent from {SELECTED_MATCH_IDS_FILE}"
    )

selected_ids_raw["SampleId"] = (
    selected_ids_raw["SampleId"].astype(str).str.strip()
)

if selected_ids_raw["SampleId"].duplicated().any():
    duplicated = selected_ids_raw.loc[
        selected_ids_raw["SampleId"].duplicated(keep=False),
        "SampleId",
    ].unique()
    raise RuntimeError(
        "Duplicated SampleId values in selected matching output: "
        f"{duplicated[:20].tolist()}"
    )

# Prefer diagnosis in the selected matching file, but always verify it against
# the canonical ML master.
if "SampleGroup" not in selected_ids_raw.columns:
    selected_ids_raw = selected_ids_raw.merge(
        master[["SampleId", "SampleGroup"]],
        on="SampleId",
        how="left",
        validate="one_to_one",
    )

selected_ids_raw["SampleGroup"] = (
    selected_ids_raw["SampleGroup"].astype(str).str.strip()
)
selected_counts = selected_ids_raw["SampleGroup"].value_counts()

if (
    len(selected_ids_raw) != 382
    or int(selected_counts.get("CN", 0)) != TARGET_CN
    or int(selected_counts.get("AD", 0)) != TARGET_AD
):
    raise RuntimeError(
        "The authoritative matching output is not 191 CN + 191 AD. "
        f"Observed N={len(selected_ids_raw)}, "
        f"CN={int(selected_counts.get('CN', 0))}, "
        f"AD={int(selected_counts.get('AD', 0))}."
    )

master_ids = set(master["SampleId"].astype(str))
missing_from_master = [
    sample_id
    for sample_id in selected_ids_raw["SampleId"]
    if sample_id not in master_ids
]
if missing_from_master:
    raise RuntimeError(
        "Selected matched participants are missing from the canonical ML "
        "master. No participant will be silently removed. Missing IDs: "
        f"{missing_from_master[:30]}"
    )

# Preserve the selected-match row order.
selected_order = {
    sample_id: position
    for position, sample_id in enumerate(
        selected_ids_raw["SampleId"].tolist()
    )
}
matched = master[
    master["SampleId"].isin(selected_order)
].copy()
matched["_selected_order"] = (
    matched["SampleId"].map(selected_order)
)
matched = matched.sort_values("_selected_order").drop(
    columns="_selected_order"
)

# Verify diagnosis agreement between matching output and the canonical master.
diagnosis_audit = selected_ids_raw[
    ["SampleId", "SampleGroup"]
].merge(
    matched[["SampleId", "SampleGroup"]],
    on="SampleId",
    suffixes=("_matching", "_master"),
    validate="one_to_one",
)
diagnosis_mismatch = diagnosis_audit[
    diagnosis_audit["SampleGroup_matching"]
    != diagnosis_audit["SampleGroup_master"]
]
if not diagnosis_mismatch.empty:
    raise RuntimeError(
        "Diagnosis mismatch between selected matching IDs and the ML master: "
        f"{diagnosis_mismatch.head(20).to_dict('records')}"
    )

required_matched_metadata = [
    "SampleGroup", "Age", "Sex", "Education", "Country"
]
missing_metadata_mask = matched[
    required_matched_metadata
].isna().any(axis=1)

# Empty strings and textual missing values are also disallowed.
for column in ["Sex", "Country"]:
    missing_metadata_mask |= (
        matched[column].astype(str).str.strip()
        .isin(["", "NA", "NaN", "nan", "None", "NULL"])
    )

if missing_metadata_mask.any():
    problem_ids = matched.loc[
        missing_metadata_mask,
        "SampleId",
    ].tolist()
    raise RuntimeError(
        "The selected 191/191 cohort contains missing model metadata in the "
        "canonical ML master. No participant will be silently removed. "
        f"Problem IDs: {problem_ids[:30]}"
    )

matched = matched.set_index("SampleId", drop=False)

matched_counts = matched["SampleGroup"].value_counts()
if (
    len(matched) != 382
    or int(matched_counts.get("CN", 0)) != 191
    or int(matched_counts.get("AD", 0)) != 191
):
    raise RuntimeError(
        "Final matched analysis frame failed the 191/191 audit."
    )

# protein_cols was defined from the complete canonical ML master in cell 3.
# Restrict only by feature quality within the selected 382 participants.
matched_numeric = matched[protein_cols].apply(
    pd.to_numeric,
    errors="coerce",
)
matched_nonmissing = matched_numeric.notna().mean(axis=0)
matched_variance = matched_numeric.var(axis=0, skipna=True)
protein_cols = [
    column
    for column in protein_cols
    if matched_nonmissing.get(column, 0) >= 0.80
    and np.isfinite(matched_variance.get(column, np.nan))
    and matched_variance.get(column, 0) > 0
]
if not protein_cols:
    raise RuntimeError(
        "No eligible protein features remain in the selected 191/191 cohort."
    )

matched_ids = matched[
    ["SampleId", "SampleGroup"]
].reset_index(drop=True)
matching_audit = matching_balance_table(
    matched.reset_index(drop=True),
    "selected_match_v5",
)
matched_source = "outputs/ML/private/Result_matching_rebuild_v5/matched_ids_SELECTED.csv"

cohort_audit = pd.DataFrame([{
    "Selected_ID_source": "deidentified metadata flag include_matched_selected via Script 15",
    "Selected_full_file": "not used",
    "Canonical_master": "reviewer-derived ReDLat_ML_gene_master_RAW.csv",
    "N_total": len(matched),
    "N_CN": int(matched_counts.get("CN", 0)),
    "N_AD": int(matched_counts.get("AD", 0)),
    "N_protein_features": len(protein_cols),
    "Missing_selected_IDs_in_master": len(missing_from_master),
    "Diagnosis_mismatches": len(diagnosis_mismatch),
    "Missing_required_metadata_rows": int(
        missing_metadata_mask.sum()
    ),
}])
historical_reconstruction_audit = cohort_audit.copy()

matched_ids.to_csv(
    MATCH_DIR / "matched_ids_selected_191.csv",
    index=False,
)
matching_audit.to_csv(
    MATCH_DIR / "matching_balance_selected_191.csv",
    index=False,
)
cohort_audit.to_csv(
    MATCH_DIR / "selected_191_cohort_audit.csv",
    index=False,
)
diagnosis_audit.to_csv(
    MATCH_DIR / "selected_191_diagnosis_audit.csv",
    index=False,
)

print("\nAuthoritative selected matched cohort loaded")
print("Source:", SELECTED_MATCH_IDS_FILE)
print(matched_counts.to_string())
print("Eligible protein features:", len(protein_cols))
print(pd.crosstab(matched["Country"], matched["SampleGroup"]))


def choose_or_create_matched_ids() -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    # Compatibility wrapper for downstream source-data code.
    return (
        matched_ids.copy(),
        matching_audit.copy(),
        matched_source,
    )


def restore_matched_cohort_191() -> pd.DataFrame:
    counts = matched["SampleGroup"].value_counts()
    if (
        len(matched) != 382
        or int(counts.get("CN", 0)) != 191
        or int(counts.get("AD", 0)) != 191
    ):
        raise RuntimeError(
            "In-memory matched cohort is no longer 191 CN + 191 AD."
        )
    return matched


# ## Reparación auditable de p-tau217
#
# La subcohorte p-tau217 se reconstruye ahora desde las fuentes históricas
# originales. El script prioriza el archivo utilizado por `Matching.R` y el
# metadata utilizado por los folds p-tau217. Solo transfiere mediciones observadas
# por `SampleId`; no imputa valores, no reemplaza participantes y no vuelve a
# hacer matching condicionado por disponibilidad del biomarcador.

# ============================================================
# 4B. RESOLVE p-tau217 FROM THE ORIGINAL METADATA SOURCES
# ============================================================
# The historical p-tau217 folds were built from
# ReDLat_CARD-proteomic_updated_all_data_11_2025.csv, not from the newer
# gene-collapsed ML master. Therefore p-tau217 availability must be resolved
# from the original metadata before defining the matched complete-case subset.

PTAU_SOURCE_COLUMN_CANDIDATES = [
    "p-tau217",
    "p.tau217",
    "p_tau217",
    "ptau217",
    "PTAU217",
    "pTau217",
]


def _read_ptau_source(path: Path, source_label: str) -> pd.DataFrame:
    """Read only SampleId, diagnosis and available p-tau217 columns."""
    try:
        header = pd.read_csv(path, nrows=0, check_names=False)
    except TypeError:
        header = pd.read_csv(path, nrows=0)
    except Exception as exc:
        print(f"[WARN] Cannot inspect p-tau source {path}: {exc}")
        return pd.DataFrame()

    columns = list(header.columns)
    id_column = _find_header_column(
        columns,
        ["SampleId", "sample_id", "participant_id", "ID"],
    )
    if id_column is None:
        return pd.DataFrame()

    group_column = _find_header_column(
        columns,
        ["SampleGroup", "Diagnosis", "diagnosis"],
    )
    ptau_columns = [
        column for column in columns
        if str(column) in PTAU_SOURCE_COLUMN_CANDIDATES
        or str(column).strip().lower() in {
            candidate.lower() for candidate in PTAU_SOURCE_COLUMN_CANDIDATES
        }
    ]
    if not ptau_columns:
        return pd.DataFrame()

    use_columns = [id_column] + ptau_columns
    if group_column is not None:
        use_columns.append(group_column)

    try:
        source = pd.read_csv(
            path,
            usecols=list(dict.fromkeys(use_columns)),
            low_memory=False,
        )
    except Exception as exc:
        print(f"[WARN] Cannot read p-tau source {path}: {exc}")
        return pd.DataFrame()

    source = source.rename(columns={id_column: "SampleId"})
    source["SampleId"] = source["SampleId"].astype(str).str.strip()

    if group_column is not None and group_column in source.columns:
        source = source.rename(columns={group_column: "SampleGroup_source"})

    numeric_ptau = pd.DataFrame(index=source.index)
    for column in ptau_columns:
        numeric_ptau[column] = pd.to_numeric(
            source[column],
            errors="coerce",
        )

    source["ptau217_value"] = numeric_ptau.bfill(axis=1).iloc[:, 0]
    source["ptau217_source_column"] = numeric_ptau.apply(
        lambda row: next(
            (
                column for column in ptau_columns
                if pd.notna(row[column])
            ),
            np.nan,
        ),
        axis=1,
    )
    source["ptau217_source_file"] = str(path)
    source["ptau217_source_label"] = source_label

    keep = [
        "SampleId",
        "ptau217_value",
        "ptau217_source_column",
        "ptau217_source_file",
        "ptau217_source_label",
    ]
    if "SampleGroup_source" in source.columns:
        keep.append("SampleGroup_source")

    source = source[keep].drop_duplicates("SampleId", keep="first")
    return source


def discover_ptau_sources() -> List[Tuple[Path, str]]:
    """Return original metadata sources in scientific priority order."""
    candidates: List[Tuple[Path, str]] = []
    seen = set()

    # 1. Exact input used by the uploaded standalone Matching.R.
    for path in discover_named_files([
        "somascan_filter_DEP_final_FDR_metadata.csv"
    ]):
        key = str(path).lower()
        if key not in seen:
            seen.add(key)
            candidates.append((path, "historical_matching_input"))

    # 2. Raw metadata used by the original p-tau217 fold script.
    for path in discover_named_files([
        "ReDLat_CARD-proteomic_updated_all_data_11_2025.csv"
    ]):
        key = str(path).lower()
        if key not in seen:
            seen.add(key)
            candidates.append((path, "historical_raw_metadata"))

    # 3. Exact historical matched output, which may retain metadata columns.
    for path in discover_matching_csvs():
        if path.name.lower() in {
            "matched_output_v2.csv",
            "matched_output_matchit.csv",
        }:
            key = str(path).lower()
            if key not in seen:
                seen.add(key)
                candidates.append((path, "historical_matched_output"))

    # 4. Current ML master only as the final fallback.
    key = str(MASTER_FILE).lower()
    if key not in seen:
        candidates.append((MASTER_FILE, "current_ml_master"))

    return candidates


def build_ptau217_lookup() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build one resolved p-tau217 value per SampleId and a full provenance audit.

    Priority follows discover_ptau_sources(): historical matching input,
    historical raw metadata, historical matched output, then current master.
    Values are never imputed; only observed measurements are transferred.
    """
    source_tables: List[pd.DataFrame] = []

    for path, label in discover_ptau_sources():
        table = _read_ptau_source(path, label)
        if not table.empty:
            print(
                f"p-tau217 source: {label} | {path} | "
                f"observed={int(table['ptau217_value'].notna().sum())}"
            )
            source_tables.append(table)

    if not source_tables:
        raise FileNotFoundError(
            "No source containing SampleId and p-tau217 was found."
        )

    audit = pd.concat(source_tables, ignore_index=True)
    audit["source_priority"] = audit["ptau217_source_label"].map({
        "historical_matching_input": 1,
        "historical_raw_metadata": 2,
        "historical_matched_output": 3,
        "current_ml_master": 4,
    }).fillna(99)

    observed = audit[audit["ptau217_value"].notna()].copy()
    observed = observed.sort_values(
        ["SampleId", "source_priority"]
    )

    resolved = observed.drop_duplicates(
        "SampleId",
        keep="first",
    ).copy()

    # Audit numerical disagreements across sources.
    disagreement_rows = []
    for sample_id, group in observed.groupby("SampleId"):
        values = group["ptau217_value"].dropna().astype(float)
        if len(values) <= 1:
            continue
        spread = float(values.max() - values.min())
        scale = max(float(values.abs().max()), 1.0)
        if spread > 1e-8 * scale:
            disagreement_rows.append({
                "SampleId": sample_id,
                "n_observed_sources": len(values),
                "minimum": float(values.min()),
                "maximum": float(values.max()),
                "absolute_spread": spread,
                "sources": " | ".join(
                    group["ptau217_source_label"].astype(str)
                ),
            })

    disagreements = pd.DataFrame(disagreement_rows)
    audit.to_csv(
        MATCH_PTAU_DIR / "ptau217_all_source_values.csv",
        index=False,
    )
    resolved.to_csv(
        MATCH_PTAU_DIR / "ptau217_resolved_lookup.csv",
        index=False,
    )
    disagreements.to_csv(
        MATCH_PTAU_DIR / "ptau217_source_disagreements.csv",
        index=False,
    )

    return resolved, audit


def prepare_matched_ptau217_cohort(
    matched_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Resolve p-tau217 for the exact 191/191 matched cohort.

    The function transfers only observed historical values by SampleId. It does
    not impute, swap participants or rematch on biomarker availability.
    """
    resolved_lookup, all_source_audit = build_ptau217_lookup()

    work = matched_df.copy()
    if work.index.name == "SampleId":
        work = work.copy()
    else:
        work = work.set_index("SampleId", drop=False)

    lookup = resolved_lookup.set_index("SampleId")
    work["ptau217_resolved"] = work["SampleId"].map(
        lookup["ptau217_value"]
    )
    work["ptau217_source_label"] = work["SampleId"].map(
        lookup["ptau217_source_label"]
    )
    work["ptau217_source_file"] = work["SampleId"].map(
        lookup["ptau217_source_file"]
    )
    work["ptau217_source_column"] = work["SampleId"].map(
        lookup["ptau217_source_column"]
    )

    current_values = pd.to_numeric(
        work[ptau_col],
        errors="coerce",
    )
    recovered = current_values.isna() & work["ptau217_resolved"].notna()

    # Use the historical observed value for all participants so that the p-tau
    # analysis is defined from one auditable source hierarchy.
    work[ptau_col] = pd.to_numeric(
        work["ptau217_resolved"],
        errors="coerce",
    )

    participant_audit = work[[
        "SampleId",
        "SampleGroup",
        ptau_col,
        "ptau217_source_label",
        "ptau217_source_file",
        "ptau217_source_column",
    ]].copy()
    participant_audit["recovered_from_historical_source"] = recovered.values
    participant_audit["included_complete_case"] = (
        participant_audit[ptau_col].notna()
    )
    participant_audit.to_csv(
        MATCH_PTAU_DIR / "matched_ptau217_participant_audit.csv",
        index=False,
    )

    counts = (
        participant_audit[
            participant_audit["included_complete_case"]
        ]["SampleGroup"]
        .value_counts()
    )
    n_cn = int(counts.get("CN", 0))
    n_ad = int(counts.get("AD", 0))

    print(
        "\nResolved matched p-tau217 complete cases:",
        f"CN={n_cn}, AD={n_ad}",
    )
    print(
        "Recovered historical p-tau217 values:",
        int(recovered.sum()),
    )

    missing = participant_audit[
        ~participant_audit["included_complete_case"]
    ].copy()
    missing.to_csv(
        MATCH_PTAU_DIR / "matched_ptau217_missing_ids.csv",
        index=False,
    )

    # Use the participant-level complete cases actually supported by the
    # available sources. Do not impute p-tau217 or replace a matched participant
    # solely to reproduce an older reported count.
    expected_cn = 158
    expected_ad = 157
    count_matches_previous_text = (
        n_cn == expected_cn and n_ad == expected_ad
    )

    discrepancy = pd.DataFrame([{
        "Observed_CN": n_cn,
        "Observed_AD": n_ad,
        "Previously_reported_CN": expected_cn,
        "Previously_reported_AD": expected_ad,
        "Matches_previous_text": count_matches_previous_text,
        "Recommended_action": (
            "No action required"
            if count_matches_previous_text
            else "Update matched p-tau217 count in Results and figure legend"
        ),
    }])
    discrepancy.to_csv(
        MATCH_PTAU_DIR / "ptau217_count_discrepancy.csv",
        index=False,
    )

    if not count_matches_previous_text:
        print(
            "\n[IMPORTANT] Audited p-tau217 complete cases: "
            f"CN={n_cn}, AD={n_ad}. "
            "The previous count of 158 CN and 157 AD is not supported by "
            "the available participant-level files. The analysis will "
            "continue with observed complete cases; no value is imputed and "
            "no participant is replaced."
        )

    if n_cn < 2 or n_ad < 2:
        raise RuntimeError(
            "The observed p-tau217 subset does not contain enough "
            "participants in both diagnostic groups."
        )

    complete = work[work[ptau_col].notna()].copy()
    return complete, participant_audit


# ## 5. Matched nested protein model
#
# Differential-abundance screening and elastic-net selection are repeated inside
# the training data. The recurrent protein set is therefore allowed to differ
# from the primary seven-protein panel.

# ============================================================
# 5. MATCHED ML HELPERS
# ============================================================

def build_dep_design(meta: pd.DataFrame) -> pd.DataFrame:
    design = pd.DataFrame(index=meta.index)
    design["Intercept"] = 1.0
    design["AD"] = (meta["SampleGroup"].astype(str) == "AD").astype(float)
    design["Age"] = pd.to_numeric(meta["Age"], errors="coerce")
    design["Education"] = pd.to_numeric(meta["Education"], errors="coerce")

    sex = pd.get_dummies(
        meta["Sex"].astype(str),
        prefix="Sex",
        drop_first=True,
        dtype=float,
    )
    country = pd.get_dummies(
        meta["Country"].astype(str),
        prefix="Country",
        drop_first=True,
        dtype=float,
    )
    design = pd.concat([design, sex, country], axis=1)
    nonconstant = [
        column for column in design.columns
        if column == "Intercept" or design[column].nunique(dropna=False) > 1
    ]
    design = design[nonconstant]
    if "AD" not in design.columns:
        raise ValueError("AD coefficient is absent from the DEP design matrix.")
    return design.replace([np.inf, -np.inf], np.nan)


def compute_training_dep(
    train_df: pd.DataFrame,
    candidate_proteins: Sequence[str],
    fold_label: str,
) -> Tuple[pd.DataFrame, List[str], Dict[str, object]]:
    metadata_vars = ["SampleGroup", "Age", "Sex", "Education", "Country"]
    missing_mask = train_df[metadata_vars].isna().any(axis=1)
    if missing_mask.any():
        problem_ids = train_df.loc[missing_mask, "SampleId"].tolist()
        raise RuntimeError(
            f"{fold_label}: matched participants would be lost because model "
            f"metadata are missing. IDs: {problem_ids[:20]}"
        )

    work = train_df.copy()
    design = build_dep_design(work)
    if design.isna().any().any():
        problem_ids = design.index[
            design.isna().any(axis=1)
        ].astype(str).tolist()
        raise RuntimeError(
            f"{fold_label}: non-finite DEP design rows: {problem_ids[:20]}"
        )

    X_design = design.to_numpy(dtype=float)
    rank = np.linalg.matrix_rank(X_design)
    n_samples, n_parameters = X_design.shape
    degrees_freedom = max(n_samples - rank, 1)
    ad_index = list(design.columns).index("AD")

    # MASTER_FILE is the reviewer-derived historical 9,638-gene normalized-RFU matrix.
    # Use it directly here because the submitted ML pipeline operated on this RAW scale.
    expression = work[list(candidate_proteins)].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)

    medians = np.nanmedian(expression, axis=0)
    medians[~np.isfinite(medians)] = 0.0
    missing_positions = np.where(~np.isfinite(expression))
    expression[missing_positions] = np.take(medians, missing_positions[1])

    xtx_inverse = np.linalg.pinv(X_design.T @ X_design)
    beta = xtx_inverse @ X_design.T @ expression
    residuals = expression - X_design @ beta
    residual_sum_squares = np.sum(residuals ** 2, axis=0)
    sigma_squared = residual_sum_squares / degrees_freedom
    standard_error = np.sqrt(
        np.maximum(sigma_squared * xtx_inverse[ad_index, ad_index], 1e-300)
    )
    log_fc = beta[ad_index, :]
    t_value = log_fc / standard_error
    p_value = 2 * stats.t.sf(np.abs(t_value), df=degrees_freedom)
    q_value = bh_adjust(p_value)

    dep = pd.DataFrame({
        "Fold": fold_label,
        "Protein": list(candidate_proteins),
        "logFC": log_fc,
        "t": t_value,
        "P.Value": p_value,
        "adj.P.Val": q_value,
    }).sort_values(["adj.P.Val", "P.Value"])

    selected = dep.loc[dep["adj.P.Val"] < DEP_FDR, "Protein"].tolist()
    fallback = False
    if len(selected) < MIN_DEP_GENES:
        selected = dep.head(FALLBACK_TOP_N_DEP)["Protein"].tolist()
        fallback = True

    summary = {
        "Fold": fold_label,
        "N_train_DEP": n_samples,
        "N_parameters": n_parameters,
        "Design_rank": rank,
        "DEP_FDR": DEP_FDR,
        "N_DEP_FDR": int((dep["adj.P.Val"] < DEP_FDR).sum()),
        "N_selected_for_ML": len(selected),
        "Fallback_topN": fallback,
    }
    return dep, selected, summary


def make_svm_search(inner_folds: int = INNER_FOLDS):
    if HAS_SKOPT:
        return BayesSearchCV(
            estimator=SVC(probability=False, random_state=RANDOM_STATE),
            search_spaces={
                "kernel": Categorical(["linear", "rbf"]),
                "C": Real(0.1, 100, prior="log-uniform"),
                "gamma": Real(1e-4, 1, prior="log-uniform"),
            },
            n_iter=SVM_SEARCH_ITER,
            cv=inner_folds,
            scoring="roc_auc",
            n_jobs=N_JOBS,
            random_state=RANDOM_STATE,
        )

    return GridSearchCV(
        estimator=SVC(probability=False, random_state=RANDOM_STATE),
        param_grid={
            "kernel": ["linear", "rbf"],
            "C": [0.1, 1, 10, 100],
            "gamma": [1e-4, 1e-3, 1e-2, 1e-1],
        },
        cv=inner_folds,
        scoring="roc_auc",
        n_jobs=N_JOBS,
    )


def prepare_train_test(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    train_imputed = imputer.fit_transform(X_train)
    test_imputed = imputer.transform(X_test)
    train_scaled = scaler.fit_transform(train_imputed)
    test_scaled = scaler.transform(test_imputed)
    return train_scaled, test_scaled


def fit_tuned_svm_probability(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
) -> Tuple[np.ndarray, Dict[str, object]]:
    X_train_scaled, X_test_scaled = prepare_train_test(X_train, X_test)

    min_class = int(y_train.value_counts().min())
    inner_folds = min(INNER_FOLDS, min_class)
    if inner_folds < 2:
        raise RuntimeError("Insufficient training observations for inner CV.")

    search = make_svm_search(inner_folds=inner_folds)
    search.fit(X_train_scaled, y_train)

    params = dict(search.best_params_)
    final_model = SVC(
        kernel=params.get("kernel", "rbf"),
        C=float(params.get("C", 1.0)),
        gamma=params.get("gamma", "scale"),
        probability=True,
        random_state=RANDOM_STATE,
    )
    final_model.fit(X_train_scaled, y_train)
    probability = final_model.predict_proba(X_test_scaled)[:, 1]
    return probability, params


def metric_row(
    fold: int,
    model: str,
    y_true: pd.Series,
    probability: np.ndarray,
) -> Dict[str, object]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(
        y_true, prediction, labels=[0, 1]
    ).ravel()
    return {
        "Fold": fold,
        "Model": model,
        "AUC": roc_auc_score(y_true, probability),
        "Accuracy": accuracy_score(y_true, prediction),
        "Sensitivity": recall_score(y_true, prediction, zero_division=0),
        "Specificity": tn / (tn + fp) if (tn + fp) > 0 else np.nan,
        "Precision": precision_score(y_true, prediction, zero_division=0),
        "F1": f1_score(y_true, prediction, zero_division=0),
        "N_test": len(y_true),
        "N_CN_test": int((y_true == 0).sum()),
        "N_AD_test": int((y_true == 1).sum()),
    }


def roc_rows(
    fold: int,
    model: str,
    y_true: pd.Series,
    probability: np.ndarray,
) -> pd.DataFrame:
    fpr, tpr, threshold = roc_curve(y_true, probability)
    auc_value = roc_auc_score(y_true, probability)
    return pd.DataFrame({
        "Fold": fold,
        "Model": model,
        "FPR": fpr,
        "TPR": tpr,
        "Threshold": threshold,
        "AUC": auc_value,
    })


def run_matched_protein_nested_cv(
    matched_df: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    cohort_counts = matched_df["SampleGroup"].value_counts()
    observed_cn = int(cohort_counts.get("CN", 0))
    observed_ad = int(cohort_counts.get("AD", 0))
    if len(matched_df) != 382 or observed_cn != 191 or observed_ad != 191:
        raise RuntimeError(
            "Panel c must start from exactly 191 CN + 191 AD; "
            f"observed CN={observed_cn}, AD={observed_ad}, "
            f"N={len(matched_df)}."
        )

    required_metadata = [
        "SampleGroup", "Age", "Sex", "Education", "Country"
    ]
    if matched_df[required_metadata].isna().any().any():
        raise RuntimeError(
            "Panel c contains missing required metadata. Participants are "
            "not permitted to be silently removed."
        )

    y = matched_df["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)
    outer_cv = StratifiedKFold(
        n_splits=OUTER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    metrics: List[Dict[str, object]] = []
    best_parameters: List[Dict[str, object]] = []
    feature_counts: List[Dict[str, object]] = []
    dep_summaries: List[Dict[str, object]] = []
    fold_assignments: List[pd.DataFrame] = []
    roc_frames: List[pd.DataFrame] = []
    outer_panels: List[List[str]] = []
    oof = pd.Series(np.nan, index=matched_df.index, dtype=float)

    for fold, (train_positions, test_positions) in enumerate(
        outer_cv.split(matched_df, y), start=1
    ):
        train_ids = matched_df.index[train_positions]
        test_ids = matched_df.index[test_positions]
        train_df = matched_df.loc[train_ids].copy()
        test_df = matched_df.loc[test_ids].copy()
        y_train = y.loc[train_ids]
        y_test = y.loc[test_ids]

        dep, dep_candidates, dep_summary = compute_training_dep(
            train_df,
            protein_cols,
            fold_label=f"Fold{fold}",
        )
        dep.to_csv(
            MATCH_RESULTS_DIR / f"DEP_Fold{fold}_all_proteins.csv",
            index=False,
        )
        pd.Series(dep_candidates, name="Protein").to_csv(
            MATCH_RESULTS_DIR / f"candidate_genes_Fold{fold}.csv",
            index=False,
        )
        dep_summaries.append(dep_summary)

        X_outer_train = train_df[dep_candidates].apply(
            pd.to_numeric, errors="coerce"
        )
        inner_selected: List[np.ndarray] = []
        inner_cv = StratifiedKFold(
            n_splits=INNER_FOLDS,
            shuffle=True,
            random_state=RANDOM_STATE,
        )

        for inner_train_positions, _ in inner_cv.split(X_outer_train, y_train):
            inner_X = X_outer_train.iloc[inner_train_positions]
            inner_y = y_train.iloc[inner_train_positions]

            imputer = SimpleImputer(strategy="median")
            scaler = StandardScaler()
            inner_scaled = scaler.fit_transform(
                imputer.fit_transform(inner_X)
            )

            selector = LogisticRegressionCV(
                penalty="elasticnet",
                solver="saga",
                Cs=np.logspace(-4, 0, 40),
                l1_ratios=[0.85, 0.90, 0.95, 1.0],
                cv=INNER_FOLDS,
                scoring="roc_auc",
                max_iter=MAX_ITER,
                tol=TOL,
                n_jobs=N_JOBS,
                random_state=RANDOM_STATE,
            )
            selector.fit(inner_scaled, inner_y)
            coefficients = np.squeeze(selector.coef_)
            selected_mask = np.abs(coefficients) > 5e-3
            inner_selected.append(
                np.asarray(dep_candidates)[selected_mask]
            )

        frequencies: Dict[str, int] = {}
        for selected in inner_selected:
            for protein in selected:
                frequencies[protein] = frequencies.get(protein, 0) + 1

        frequency_series = (
            pd.Series(frequencies, dtype=float) / INNER_FOLDS
        ).sort_values(ascending=False)

        frequency_series.to_csv(
            MATCH_RESULTS_DIR / f"frequency_fold_{fold}.csv",
            header=["Selection_frequency"],
        )

        panel = frequency_series[
            frequency_series >= INNER_STABILITY
        ].index.tolist()
        if len(panel) == 0:
            panel = frequency_series.head(10).index.tolist()
        if len(panel) == 0:
            panel = dep_candidates[:10]

        outer_panels.append(panel)
        pd.Series(panel, name="Protein").to_csv(
            MATCH_RESULTS_DIR / f"panel_fold_{fold}.csv",
            index=False,
        )

        X_train = train_df[panel].apply(pd.to_numeric, errors="coerce")
        X_test = test_df[panel].apply(pd.to_numeric, errors="coerce")
        probability, params = fit_tuned_svm_probability(
            X_train, y_train, X_test
        )
        oof.loc[test_ids] = probability

        metrics.append(
            metric_row(fold, "Matched proteins", y_test, probability)
        )
        roc_frames.append(
            roc_rows(fold, "Matched proteins", y_test, probability)
        )
        best_parameters.append({
            "Fold": fold,
            **params,
            "DEP_candidates": len(dep_candidates),
            "Consensus_panel": len(panel),
        })
        feature_counts.append({
            "Fold": fold,
            "DEP_candidates": len(dep_candidates),
            "Consensus_panel": len(panel),
            "Inner_stability_threshold": INNER_STABILITY,
        })
        fold_assignments.append(pd.DataFrame({
            "SampleId": list(train_ids) + list(test_ids),
            "Fold": fold,
            "Set": ["train"] * len(train_ids) + ["test"] * len(test_ids),
        }))

        print(
            f"Matched fold {fold}: AUC="
            f"{metrics[-1]['AUC']:.3f}; panel={len(panel)}"
        )

    metrics_df = pd.DataFrame(metrics)
    roc_df = pd.concat(roc_frames, ignore_index=True)
    best_params_df = pd.DataFrame(best_parameters)
    feature_counts_df = pd.DataFrame(feature_counts)
    dep_summary_df = pd.DataFrame(dep_summaries)
    fold_assignment_df = pd.concat(fold_assignments, ignore_index=True)

    test_assignments = fold_assignment_df[
        fold_assignment_df["Set"] == "test"
    ].copy()
    test_counts = test_assignments["SampleId"].value_counts()
    if (
        len(test_assignments) != 382
        or test_assignments["SampleId"].nunique() != 382
        or not (test_counts == 1).all()
    ):
        raise RuntimeError(
            "Outer-fold test assignments do not cover each of the 382 "
            "matched participants exactly once."
        )

    oof_complete = oof.notna()
    if int(oof_complete.sum()) != 382:
        raise RuntimeError(
            f"Only {int(oof_complete.sum())} of 382 matched participants "
            "received an out-of-fold prediction."
        )

    fold_cohort_audit = pd.DataFrame([{
        "N_total": len(matched_df),
        "N_CN": observed_cn,
        "N_AD": observed_ad,
        "N_unique_outer_test_IDs": (
            test_assignments["SampleId"].nunique()
        ),
        "N_nonmissing_OOF_predictions": int(oof_complete.sum()),
        "Protein_universe": len(protein_cols),
    }])
    fold_cohort_audit.to_csv(
        MATCH_RESULTS_DIR / "canonical_191_model_audit.csv",
        index=False,
    )

    panel_frequency: Dict[str, int] = {}
    for panel in outer_panels:
        for protein in panel:
            panel_frequency[protein] = panel_frequency.get(protein, 0) + 1
    panel_frequency_series = (
        pd.Series(panel_frequency, dtype=float) / len(outer_panels)
    ).sort_values(ascending=False)

    # Strict recurrent-panel definition: retain every protein selected in at
    # least 80% of the five outer-fold panels. Do not pad the panel to an
    # arbitrary size, because doing so would introduce proteins below the
    # prespecified stability threshold.
    final_panel = panel_frequency_series[
        panel_frequency_series >= FINAL_PANEL_THRESHOLD
    ].index.tolist()

    if len(final_panel) == 0:
        raise RuntimeError(
            "No protein met the prespecified outer-fold selection-frequency "
            f"threshold of {FINAL_PANEL_THRESHOLD:.2f}. The threshold will "
            "not be relaxed automatically."
        )

    final_panel_audit = panel_frequency_series.rename(
        "Selection_frequency"
    ).reset_index().rename(columns={"index": "Protein"})
    final_panel_audit["Meets_final_threshold"] = (
        final_panel_audit["Selection_frequency"] >= FINAL_PANEL_THRESHOLD
    )
    final_panel_audit["Included_in_final_panel"] = (
        final_panel_audit["Protein"].isin(final_panel)
    )

    if not (
        final_panel_audit["Meets_final_threshold"]
        == final_panel_audit["Included_in_final_panel"]
    ).all():
        raise RuntimeError(
            "The recurrent matched panel does not exactly match the "
            "prespecified >=80% selection-frequency rule."
        )

    final_panel_audit.to_csv(
        MATCH_FINAL_PANEL_DIR / "panel_final_freq08_audit.csv",
        index=False,
    )
    pd.Series(final_panel, name="Protein").to_csv(
        MATCH_FINAL_PANEL_DIR / "panel_final_freq08.csv",
        index=False,
    )
    print(
        "Strict matched recurrent panel "
        f"(selection frequency >= {FINAL_PANEL_THRESHOLD:.2f}): "
        f"{final_panel}"
    )

    # Final permutation importance is descriptive; CV performance remains based
    # only on outer held-out folds.
    X_final = matched_df[final_panel].apply(pd.to_numeric, errors="coerce")
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    X_final_scaled = scaler.fit_transform(imputer.fit_transform(X_final))
    best_kernel = best_params_df["kernel"].mode().iloc[0]
    best_c = float(pd.to_numeric(best_params_df["C"]).median())
    best_gamma = pd.to_numeric(
        best_params_df["gamma"], errors="coerce"
    ).median()
    if not np.isfinite(best_gamma):
        best_gamma = "scale"

    final_model = SVC(
        kernel=best_kernel,
        C=best_c,
        gamma=best_gamma,
        probability=True,
        random_state=RANDOM_STATE,
    )
    final_model.fit(X_final_scaled, y)
    permutation = permutation_importance(
        final_model,
        X_final_scaled,
        y,
        n_repeats=100,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=N_JOBS,
    )
    permutation_df = pd.DataFrame({
        "Proteins": final_panel,
        "Permutation": permutation.importances_mean,
        "Permutation_sd": permutation.importances_std,
        "Selection_frequency": [
            panel_frequency_series.get(protein, np.nan)
            for protein in final_panel
        ],
        "In_primary_panel": [
            protein in PRIMARY_PANEL for protein in final_panel
        ],
    }).sort_values("Permutation", ascending=False)

    metrics_df.to_csv(MATCH_RESULTS_DIR / "metrics_nestedCV.csv", index=False)
    roc_df.to_csv(MATCH_RESULTS_DIR / "roc_curves.csv", index=False)
    best_params_df.to_csv(
        MATCH_RESULTS_DIR / "best_params_outer.csv", index=False
    )
    feature_counts_df.to_csv(
        MATCH_RESULTS_DIR / "feature_counts.csv", index=False
    )
    dep_summary_df.to_csv(
        MATCH_RESULTS_DIR / "nested_DEP_summary.csv", index=False
    )
    fold_assignment_df.to_csv(
        MATCH_RESULTS_DIR / "fold_assignments.csv", index=False
    )
    panel_frequency_series.to_csv(
        MATCH_RESULTS_DIR / "outer_panel_frequency.csv",
        header=["Selection_frequency"],
    )
    permutation_df.to_csv(
        MATCH_FINAL_PANEL_DIR / "permutation_importance_freq08.csv",
        index=False,
    )
    pd.DataFrame({
        "SampleId": matched_df.index,
        "y_true": y.loc[matched_df.index].values,
        "oof_probability": oof.loc[matched_df.index].values,
        "SampleGroup": matched_df["SampleGroup"].values,
        "Country": matched_df["Country"].values,
    }).to_csv(MATCH_RESULTS_DIR / "oof_predictions.csv", index=False)

    summary = metrics_df.groupby("Model").agg(
        Mean_AUC=("AUC", "mean"),
        SD_AUC=("AUC", "std"),
        Mean_Accuracy=("Accuracy", "mean"),
        Mean_Sensitivity=("Sensitivity", "mean"),
        Mean_Specificity=("Specificity", "mean"),
    ).reset_index()
    summary.to_csv(
        MATCH_RESULTS_DIR / "metrics_summary_mean_sd.csv", index=False
    )

    return {
        "metrics": metrics_df,
        "roc": roc_df,
        "best_params": best_params_df,
        "feature_counts": feature_counts_df,
        "dep_summary": dep_summary_df,
        "fold_assignments": fold_assignment_df,
        "panel_frequency": panel_frequency_series.reset_index(
            name="Selection_frequency"
        ).rename(columns={"index": "Protein"}),
        "final_panel": pd.DataFrame({"Protein": final_panel}),
        "permutation": permutation_df,
        "oof": pd.DataFrame({
            "SampleId": matched_df.index,
            "y_true": y.loc[matched_df.index].values,
            "oof_probability": oof.loc[matched_df.index].values,
        }),
    }


matched = restore_matched_cohort_191()

required_matched_outputs = [
    MATCH_RESULTS_DIR / "metrics_nestedCV.csv",
    MATCH_RESULTS_DIR / "roc_curves.csv",
    MATCH_RESULTS_DIR / "fold_assignments.csv",
    MATCH_FINAL_PANEL_DIR / "panel_final_freq08.csv",
    MATCH_FINAL_PANEL_DIR / "permutation_importance_freq08.csv",
]

if FORCE_RECOMPUTE_MATCHED or not all(path.exists() for path in required_matched_outputs):
    matched_results = run_matched_protein_nested_cv(matched)
else:
    print("Loading cached matched protein-only results.")
    matched_results = {
        "metrics": pd.read_csv(MATCH_RESULTS_DIR / "metrics_nestedCV.csv"),
        "roc": pd.read_csv(MATCH_RESULTS_DIR / "roc_curves.csv"),
        "fold_assignments": pd.read_csv(
            MATCH_RESULTS_DIR / "fold_assignments.csv"
        ),
        "final_panel": pd.read_csv(
            MATCH_FINAL_PANEL_DIR / "panel_final_freq08.csv"
        ),
        "permutation": pd.read_csv(
            MATCH_FINAL_PANEL_DIR / "permutation_importance_freq08.csv"
        ),
        "oof": pd.read_csv(MATCH_RESULTS_DIR / "oof_predictions.csv"),
    }


# ------------------------------------------------------------
# COMPATIBILITY CHECK FOR MATCHED OUTPUTS
# ------------------------------------------------------------
# Ensure metrics DataFrame contains a Model column for downstream aggregation.
if "Model" not in matched_results["metrics"].columns:
    matched_results["metrics"] = matched_results["metrics"].copy()
    matched_results["metrics"].insert(1, "Model", "Matched proteins")

# Older panel exports may contain an unnamed first column instead of Protein.
if "Protein" not in matched_results["final_panel"].columns:
    candidate_columns = [
        column for column in matched_results["final_panel"].columns
        if not str(column).lower().startswith("unnamed")
    ]
    if len(candidate_columns) == 1:
        matched_results["final_panel"] = matched_results["final_panel"].rename(
            columns={candidate_columns[0]: "Protein"}
        )
    elif matched_results["final_panel"].shape[1] >= 1:
        matched_results["final_panel"] = matched_results["final_panel"].rename(
            columns={matched_results["final_panel"].columns[-1]: "Protein"}
        )

# Ensure the ROC table also carries a model label.
if "Model" not in matched_results["roc"].columns:
    matched_results["roc"] = matched_results["roc"].copy()
    matched_results["roc"].insert(1, "Model", "Matched proteins")

matched_panel = (
    MATCHED_PANEL_OVERRIDE
    if MATCHED_PANEL_OVERRIDE is not None
    else matched_results["final_panel"]["Protein"].astype(str).tolist()
)
print("Matched consensus panel:", matched_panel)

# Final safeguard: a cached or manually supplied panel must still obey the
# prespecified >=80% outer-fold selection-frequency rule.
_panel_frequency_path = MATCH_RESULTS_DIR / "outer_panel_frequency.csv"
if _panel_frequency_path.exists():
    _panel_frequency = pd.read_csv(_panel_frequency_path, index_col=0).iloc[:, 0]
    _panel_frequency.index = _panel_frequency.index.astype(str)
    _expected_panel = set(
        _panel_frequency[_panel_frequency >= FINAL_PANEL_THRESHOLD].index
    )
    _observed_panel = set(map(str, matched_panel))
    if _observed_panel != _expected_panel:
        raise RuntimeError(
            "Cached/overridden matched panel violates the prespecified "
            ">=80% outer-fold selection-frequency rule. "
            f"Expected={sorted(_expected_panel)}; "
            f"observed={sorted(_observed_panel)}"
        )

print(
    matched_results["metrics"]
    .groupby("Model")["AUC"]
    .agg(["mean", "std", "count"])
)


# ## 5B. Primary seven-protein panel in the matched cohort
#
# This complementary sensitivity analysis keeps the seven proteins from the main
# analysis fixed. No differential-abundance filtering or feature reselection is
# performed. SVM hyperparameters remain tuned inside the training data, and the
# same outer folds used by the matched nested model are retained for direct
# comparison.

# ============================================================
# 5B. PRIMARY FIXED-PANEL MODEL IN THE MATCHED COHORT
# ============================================================

def run_primary_fixed_panel_cv(
    matched_df: pd.DataFrame,
    fold_assignments: pd.DataFrame,
    panel: Sequence[str] = PRIMARY_PANEL,
) -> Dict[str, pd.DataFrame]:
    panel = list(panel)
    missing = [protein for protein in panel if protein not in matched_df.columns]
    if missing:
        raise RuntimeError(
            "Primary panel proteins are missing from the canonical ML master: "
            f"{missing}"
        )
    if len(panel) != 7 or len(set(panel)) != 7:
        raise RuntimeError(
            f"PRIMARY_PANEL must contain exactly seven unique proteins: {panel}"
        )

    counts = matched_df["SampleGroup"].value_counts()
    if (
        len(matched_df) != 382
        or int(counts.get("CN", 0)) != 191
        or int(counts.get("AD", 0)) != 191
    ):
        raise RuntimeError(
            "The primary fixed-panel sensitivity analysis must use "
            "191 CN + 191 AD."
        )

    assignment = fold_assignments.copy()
    assignment["SampleId"] = assignment["SampleId"].astype(str)
    y_all = matched_df["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)

    metrics: List[Dict[str, object]] = []
    roc_frames: List[pd.DataFrame] = []
    parameter_rows: List[Dict[str, object]] = []
    oof = pd.Series(np.nan, index=matched_df.index, dtype=float)

    for fold in sorted(assignment["Fold"].unique()):
        train_ids = assignment.loc[
            (assignment["Fold"] == fold)
            & (assignment["Set"] == "train"),
            "SampleId",
        ].tolist()
        test_ids = assignment.loc[
            (assignment["Fold"] == fold)
            & (assignment["Set"] == "test"),
            "SampleId",
        ].tolist()

        train_ids = [sample for sample in train_ids if sample in matched_df.index]
        test_ids = [sample for sample in test_ids if sample in matched_df.index]

        if not train_ids or not test_ids:
            raise RuntimeError(f"Missing train or test IDs in fold {fold}.")

        y_train = y_all.loc[train_ids]
        y_test = y_all.loc[test_ids]
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            raise RuntimeError(
                f"Primary fixed-panel fold {fold} does not contain both classes."
            )

        X_train = matched_df.loc[train_ids, panel].apply(
            pd.to_numeric, errors="coerce"
        )
        X_test = matched_df.loc[test_ids, panel].apply(
            pd.to_numeric, errors="coerce"
        )
        probability, params = fit_tuned_svm_probability(
            X_train, y_train, X_test
        )
        oof.loc[test_ids] = probability

        metrics.append(
            metric_row(
                fold,
                "Primary fixed panel",
                y_test,
                probability,
            )
        )
        roc_frames.append(
            roc_rows(
                fold,
                "Primary fixed panel",
                y_test,
                probability,
            )
        )
        parameter_rows.append({
            "Fold": fold,
            "Model": "Primary fixed panel",
            **params,
            "N_proteins": len(panel),
            "Feature_reselection": False,
        })

        print(
            f"Primary fixed-panel fold {fold}: "
            f"AUC={metrics[-1]['AUC']:.3f}"
        )

    if int(oof.notna().sum()) != 382:
        raise RuntimeError(
            "Primary fixed-panel OOF predictions do not cover all "
            "382 matched participants."
        )

    metrics_df = pd.DataFrame(metrics)
    roc_df = pd.concat(roc_frames, ignore_index=True)
    params_df = pd.DataFrame(parameter_rows)
    oof_df = pd.DataFrame({
        "SampleId": matched_df.index,
        "SampleGroup": matched_df["SampleGroup"].values,
        "y_true": y_all.loc[matched_df.index].values,
        "oof_probability": oof.loc[matched_df.index].values,
        "Model": "Primary fixed panel",
    })

    # Descriptive feature importance fitted on all matched participants.
    # It is not used to estimate cross-validated discrimination.
    X_full = matched_df[panel].apply(pd.to_numeric, errors="coerce")
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    X_full_scaled = scaler.fit_transform(imputer.fit_transform(X_full))

    best_kernel = params_df["kernel"].mode().iloc[0]
    best_c = float(pd.to_numeric(params_df["C"]).median())
    best_gamma = pd.to_numeric(
        params_df["gamma"], errors="coerce"
    ).median()
    if not np.isfinite(best_gamma):
        best_gamma = "scale"

    final_model = SVC(
        kernel=best_kernel,
        C=best_c,
        gamma=best_gamma,
        probability=True,
        random_state=RANDOM_STATE,
    )
    final_model.fit(X_full_scaled, y_all)
    permutation = permutation_importance(
        final_model,
        X_full_scaled,
        y_all,
        n_repeats=100,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=N_JOBS,
    )
    importance_df = pd.DataFrame({
        "Protein": panel,
        "Importance": permutation.importances_mean,
        "Importance_SD": permutation.importances_std,
    }).sort_values("Importance", ascending=False)

    summary_df = metrics_df.groupby("Model").agg(
        Mean_AUC=("AUC", "mean"),
        SD_AUC=("AUC", "std"),
        Mean_Accuracy=("Accuracy", "mean"),
        Mean_Sensitivity=("Sensitivity", "mean"),
        Mean_Specificity=("Specificity", "mean"),
    ).reset_index()

    panel_df = pd.DataFrame({
        "Protein": panel,
        "Fixed_from_primary_analysis": True,
        "Feature_reselection_in_matched_cohort": False,
    })

    metrics_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "metrics_nestedCV_fixed_panel.csv",
        index=False,
    )
    roc_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "roc_curves_fixed_panel.csv",
        index=False,
    )
    params_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "best_params_outer_fixed_panel.csv",
        index=False,
    )
    oof_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "oof_predictions_fixed_panel.csv",
        index=False,
    )
    importance_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "permutation_importance_fixed_panel.csv",
        index=False,
    )
    summary_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "metrics_summary_fixed_panel.csv",
        index=False,
    )
    panel_df.to_csv(
        PRIMARY_FIXED_MATCH_DIR / "primary_panel_definition.csv",
        index=False,
    )

    return {
        "metrics": metrics_df,
        "roc": roc_df,
        "params": params_df,
        "oof": oof_df,
        "importance": importance_df,
        "summary": summary_df,
        "panel": panel_df,
    }


required_fixed_outputs = [
    PRIMARY_FIXED_MATCH_DIR / "metrics_nestedCV_fixed_panel.csv",
    PRIMARY_FIXED_MATCH_DIR / "roc_curves_fixed_panel.csv",
    PRIMARY_FIXED_MATCH_DIR / "oof_predictions_fixed_panel.csv",
    PRIMARY_FIXED_MATCH_DIR / "primary_panel_definition.csv",
]

if FORCE_RECOMPUTE_FIXED_PANEL or not all(
    path.exists() for path in required_fixed_outputs
):
    primary_fixed_results = run_primary_fixed_panel_cv(
        matched,
        matched_results["fold_assignments"],
        PRIMARY_PANEL,
    )
else:
    print("Loading cached primary fixed-panel results.")
    primary_fixed_results = {
        "metrics": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "metrics_nestedCV_fixed_panel.csv"
        ),
        "roc": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "roc_curves_fixed_panel.csv"
        ),
        "params": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "best_params_outer_fixed_panel.csv"
        ),
        "oof": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "oof_predictions_fixed_panel.csv"
        ),
        "importance": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "permutation_importance_fixed_panel.csv"
        ),
        "summary": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "metrics_summary_fixed_panel.csv"
        ),
        "panel": pd.read_csv(
            PRIMARY_FIXED_MATCH_DIR / "primary_panel_definition.csv"
        ),
    }

# Relabel the newly selected matched model so it cannot be mistaken for the
# fixed primary panel.
matched_nested_roc = matched_results["roc"].copy()
matched_nested_roc["Model"] = "Matched nested"
matched_fixed_comparison_roc = pd.concat(
    [matched_nested_roc, primary_fixed_results["roc"]],
    ignore_index=True,
)

matched_panel_overlap = pd.DataFrame({
    "Protein": sorted(set(PRIMARY_PANEL) | set(matched_panel))
})
matched_panel_overlap["In_primary_panel"] = (
    matched_panel_overlap["Protein"].isin(PRIMARY_PANEL)
)
matched_panel_overlap["In_matched_recurrent_panel"] = (
    matched_panel_overlap["Protein"].isin(matched_panel)
)
matched_panel_overlap["Shared"] = (
    matched_panel_overlap["In_primary_panel"]
    & matched_panel_overlap["In_matched_recurrent_panel"]
)
matched_panel_overlap.to_csv(
    PRIMARY_FIXED_MATCH_DIR / "primary_vs_matched_panel_overlap.csv",
    index=False,
)

shared_proteins = matched_panel_overlap.loc[
    matched_panel_overlap["Shared"], "Protein"
].tolist()
print("Primary fixed panel:", PRIMARY_PANEL)
print("Matched recurrent panel:", matched_panel)
print(
    f"Shared proteins: {len(shared_proteins)}/{len(PRIMARY_PANEL)} "
    f"relative to the primary panel — {shared_proteins}"
)
print(primary_fixed_results["summary"].to_string(index=False))


# ## 6. Primary fixed panel and p-tau217 in the matched subgroup
#
# The protein component is the exact seven-protein panel from the main analysis.
# No protein feature reselection is performed in this block. Models are evaluated
# only among matched participants with observed p-tau217.

# ============================================================
# 6. DELONG TEST
# ============================================================

def compute_midrank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    sorted_values = values[order]
    n = len(values)
    midranks = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and sorted_values[j] == sorted_values[i]:
            j += 1
        midranks[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    output = np.empty(n, dtype=float)
    output[order] = midranks
    return output


def fast_delong(
    predictions_sorted_transposed: np.ndarray,
    label_1_count: int,
) -> Tuple[np.ndarray, np.ndarray]:
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    classifiers = predictions_sorted_transposed.shape[0]

    positive_midranks = np.empty((classifiers, m))
    negative_midranks = np.empty((classifiers, n))
    all_midranks = np.empty((classifiers, m + n))

    for classifier in range(classifiers):
        positive_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier, :m]
        )
        negative_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier, m:]
        )
        all_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier]
        )

    aucs = (
        all_midranks[:, :m].sum(axis=1) / m - (m + 1) / 2
    ) / n
    v01 = (
        all_midranks[:, :m] - positive_midranks
    ) / n
    v10 = 1 - (
        all_midranks[:, m:] - negative_midranks
    ) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    covariance = sx / m + sy / n
    return aucs, covariance


def delong_pairwise(
    y_true: np.ndarray,
    score_a: np.ndarray,
    score_b: np.ndarray,
) -> Dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    score_a = np.asarray(score_a, dtype=float)
    score_b = np.asarray(score_b, dtype=float)

    order = np.argsort(-y_true)
    label_1_count = int(y_true.sum())
    predictions = np.vstack([score_a, score_b])[:, order]

    aucs, covariance = fast_delong(predictions, label_1_count)
    contrast = np.array([[1, -1]], dtype=float)
    variance = float(np.squeeze(contrast @ covariance @ contrast.T))
    if variance <= 0 or not np.isfinite(variance):
        p_value = np.nan
        z_value = np.nan
    else:
        z_value = float((aucs[0] - aucs[1]) / np.sqrt(variance))
        p_value = float(2 * stats.norm.sf(abs(z_value)))

    return {
        "AUC_1": float(aucs[0]),
        "AUC_2": float(aucs[1]),
        "Delta_AUC": float(aucs[1] - aucs[0]),
        "Z": z_value,
        "P": p_value,
    }

# ============================================================
# 6B. RUN MATCHED PTAU MODELS
# ============================================================

def run_matched_ptau_models(
    matched_df: pd.DataFrame,
    fold_assignments: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    try:
        work, ptau_participant_audit = prepare_matched_ptau217_cohort(
            matched_df
        )
    except RuntimeError as exc:
        message = str(exc).lower()
        legacy_count_error = (
            "exact historical sources still do not support the manuscript"
            in message
            or "matched p-tau217 subcohort does not reproduce the manuscript"
            in message
        )

        if not legacy_count_error:
            raise

        print(
            "\n[KERNEL-SAFE FALLBACK] An older in-memory resolver attempted "
            "to force the previous 158/157 count. Continuing from the "
            "participant audit already written before that exception."
        )

        audit_file = (
            MATCH_PTAU_DIR / "matched_ptau217_participant_audit.csv"
        )
        if not audit_file.exists():
            raise RuntimeError(
                "The legacy resolver raised before a participant audit was "
                f"available: {audit_file}"
            ) from exc

        ptau_participant_audit = pd.read_csv(
            audit_file,
            low_memory=False,
        )
        required_audit_columns = {
            "SampleId",
            "SampleGroup",
            ptau_col,
            "included_complete_case",
        }
        missing_audit_columns = (
            required_audit_columns
            - set(ptau_participant_audit.columns)
        )
        if missing_audit_columns:
            raise RuntimeError(
                "The p-tau217 participant audit is incomplete. Missing: "
                f"{sorted(missing_audit_columns)}"
            ) from exc

        included_flag = (
            ptau_participant_audit["included_complete_case"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "yes"})
        )
        observed_audit = ptau_participant_audit[
            included_flag
        ].copy()
        observed_audit["SampleId"] = (
            observed_audit["SampleId"].astype(str)
        )
        observed_audit[ptau_col] = pd.to_numeric(
            observed_audit[ptau_col],
            errors="coerce",
        )
        observed_audit = observed_audit.dropna(
            subset=[ptau_col]
        )

        work = matched_df.copy()
        if work.index.name != "SampleId":
            work = work.set_index("SampleId", drop=False)

        observed_lookup = observed_audit.set_index("SampleId")
        work[ptau_col] = work["SampleId"].astype(str).map(
            observed_lookup[ptau_col]
        )
        work = work[work[ptau_col].notna()].copy()

        fallback_counts = work["SampleGroup"].value_counts()
        fallback_cn = int(fallback_counts.get("CN", 0))
        fallback_ad = int(fallback_counts.get("AD", 0))

        print(
            "Observed p-tau217 cohort restored by fallback: "
            f"CN={fallback_cn}, AD={fallback_ad}"
        )

        if fallback_cn < 2 or fallback_ad < 2:
            raise RuntimeError(
                "The fallback p-tau217 cohort is not analysable: "
                f"CN={fallback_cn}, AD={fallback_ad}"
            ) from exc

    y_all = work["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)

    metrics: List[Dict[str, object]] = []
    roc_frames: List[pd.DataFrame] = []
    parameter_rows: List[Dict[str, object]] = []
    oof = pd.DataFrame(
        index=work.index,
        columns=["Primary panel", "p-tau217", "Combined"],
        dtype=float,
    )

    assignment = fold_assignments.copy()
    assignment["SampleId"] = assignment["SampleId"].astype(str)

    for fold in sorted(assignment["Fold"].unique()):
        train_ids = assignment.loc[
            (assignment["Fold"] == fold)
            & (assignment["Set"] == "train"),
            "SampleId",
        ]
        test_ids = assignment.loc[
            (assignment["Fold"] == fold)
            & (assignment["Set"] == "test"),
            "SampleId",
        ]

        train_ids = [sample for sample in train_ids if sample in work.index]
        test_ids = [sample for sample in test_ids if sample in work.index]
        if not train_ids or not test_ids:
            continue

        y_train = y_all.loc[train_ids]
        y_test = y_all.loc[test_ids]
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            raise RuntimeError(
                f"Matched p-tau fold {fold} does not contain both classes."
            )

        panel = list(PRIMARY_PANEL)
        missing_panel = [
            protein for protein in panel if protein not in work.columns
        ]
        if missing_panel:
            raise RuntimeError(
                "Primary panel proteins are missing from the matched "
                f"p-tau217 cohort: {missing_panel}"
            )

        model_matrices = {
            "Primary panel": (
                work.loc[train_ids, panel].apply(
                    pd.to_numeric, errors="coerce"
                ),
                work.loc[test_ids, panel].apply(
                    pd.to_numeric, errors="coerce"
                ),
            ),
            "p-tau217": (
                work.loc[train_ids, [ptau_col]].apply(
                    pd.to_numeric, errors="coerce"
                ),
                work.loc[test_ids, [ptau_col]].apply(
                    pd.to_numeric, errors="coerce"
                ),
            ),
            "Combined": (
                pd.concat([
                    work.loc[train_ids, panel].apply(
                        pd.to_numeric, errors="coerce"
                    ),
                    work.loc[train_ids, [ptau_col]].apply(
                        pd.to_numeric, errors="coerce"
                    ),
                ], axis=1),
                pd.concat([
                    work.loc[test_ids, panel].apply(
                        pd.to_numeric, errors="coerce"
                    ),
                    work.loc[test_ids, [ptau_col]].apply(
                        pd.to_numeric, errors="coerce"
                    ),
                ], axis=1),
            ),
        }

        for model_name, (X_train, X_test) in model_matrices.items():
            probability, params = fit_tuned_svm_probability(
                X_train, y_train, X_test
            )
            oof.loc[test_ids, model_name] = probability
            metrics.append(
                metric_row(fold, model_name, y_test, probability)
            )
            roc_frames.append(
                roc_rows(fold, model_name, y_test, probability)
            )
            parameter_rows.append({
                "Fold": fold,
                "Model": model_name,
                **params,
                "N_proteins": len(panel) if model_name != "p-tau217" else 0,
            })

        print(
            f"Matched p-tau fold {fold}: "
            + ", ".join(
                f"{row['Model']}={row['AUC']:.3f}"
                for row in metrics if row["Fold"] == fold
            )
        )

    if oof.isna().any().any():
        missing = oof[oof.isna().any(axis=1)]
        raise RuntimeError(
            f"Missing matched p-tau OOF predictions for {len(missing)} participants."
        )

    metrics_df = pd.DataFrame(metrics)
    roc_df = pd.concat(roc_frames, ignore_index=True)
    parameters_df = pd.DataFrame(parameter_rows)
    oof_export = oof.copy()
    oof_export.insert(0, "y_true", y_all.loc[oof.index])
    oof_export.insert(0, "SampleId", oof.index)

    comparisons = []
    for comparison, model_a, model_b in [
        ("PrimaryPanel_vs_pTau", "Primary panel", "p-tau217"),
        ("PrimaryPanel_vs_Combined", "Primary panel", "Combined"),
        ("pTau_vs_Combined", "p-tau217", "Combined"),
    ]:
        result = delong_pairwise(
            oof_export["y_true"].values,
            oof_export[model_a].values,
            oof_export[model_b].values,
        )
        comparisons.append({
            "Comparison": comparison,
            **result,
        })
    delong_df = pd.DataFrame(comparisons)

    metrics_df.to_csv(MATCH_PTAU_DIR / "metrics_all_models.csv", index=False)
    roc_df.to_csv(MATCH_PTAU_DIR / "roc_curves.csv", index=False)
    parameters_df.to_csv(
        MATCH_PTAU_DIR / "best_params_outer.csv", index=False
    )
    oof_export.to_csv(
        MATCH_PTAU_DIR / "oof_predictions.csv", index=False
    )
    delong_df.to_csv(
        MATCH_PTAU_DIR / "delong_results.csv", index=False
    )

    cohort_summary = pd.DataFrame({
        "N_total": [len(work)],
        "N_CN": [int((y_all == 0).sum())],
        "N_AD": [int((y_all == 1).sum())],
        "p_tau_column": [ptau_col],
        "protein_panel": [";".join(PRIMARY_PANEL)],
        "feature_reselection": [False],
    })
    cohort_summary.to_csv(
        MATCH_PTAU_DIR / "cohort_summary.csv", index=False
    )

    summary = metrics_df.groupby("Model").agg(
        Mean_AUC=("AUC", "mean"),
        SD_AUC=("AUC", "std"),
        Mean_Accuracy=("Accuracy", "mean"),
        Mean_Sensitivity=("Sensitivity", "mean"),
        Mean_Specificity=("Specificity", "mean"),
    ).reset_index()
    summary.to_csv(
        MATCH_PTAU_DIR / "model_comparison_summary.csv", index=False
    )

    return {
        "metrics": metrics_df,
        "roc": roc_df,
        "params": parameters_df,
        "oof": oof_export,
        "delong": delong_df,
        "cohort": cohort_summary,
        "ptau_participant_audit": ptau_participant_audit,
    }


matched = restore_matched_cohort_191()

print("Running matched p-tau217 analysis block under validated configuration.")
print(
    "Observed complete-case counts are accepted; "
    "158 CN + 157 AD is not a hard requirement."
)

required_ptau_outputs = [
    MATCH_PTAU_DIR / "metrics_all_models.csv",
    MATCH_PTAU_DIR / "roc_curves.csv",
    MATCH_PTAU_DIR / "oof_predictions.csv",
    MATCH_PTAU_DIR / "delong_results.csv",
]

if FORCE_RECOMPUTE_MATCHED_PTAU or not all(
    path.exists() for path in required_ptau_outputs
):
    matched_ptau_results = run_matched_ptau_models(
        matched,
        matched_results["fold_assignments"],
    )
else:
    print("Loading cached matched p-tau217 results.")
    matched_ptau_results = {
        "metrics": pd.read_csv(
            MATCH_PTAU_DIR / "metrics_all_models.csv"
        ),
        "roc": pd.read_csv(MATCH_PTAU_DIR / "roc_curves.csv"),
        "oof": pd.read_csv(MATCH_PTAU_DIR / "oof_predictions.csv"),
        "delong": pd.read_csv(MATCH_PTAU_DIR / "delong_results.csv"),
        "cohort": pd.read_csv(MATCH_PTAU_DIR / "cohort_summary.csv"),
        "ptau_participant_audit": pd.read_csv(
            MATCH_PTAU_DIR / "matched_ptau217_participant_audit.csv"
        ),
    }

print(
    matched_ptau_results["metrics"]
    .groupby("Model")["AUC"]
    .agg(["mean", "std", "count"])
)
print(matched_ptau_results["delong"])

ptau_counts = matched_ptau_results["cohort"].iloc[0]
observed_ptau_cn = int(ptau_counts["N_CN"])
observed_ptau_ad = int(ptau_counts["N_AD"])

ptau_count_matches_previous_text = (
    observed_ptau_cn == 158 and observed_ptau_ad == 157
)

if ptau_count_matches_previous_text:
    print("Matched p-tau217 counts: 158 CN + 157 AD")
else:
    print(
        "Matched p-tau217 counts used in this reproducible analysis: "
        f"{observed_ptau_cn} CN + {observed_ptau_ad} AD"
    )
    print(
        "[MANUSCRIPT UPDATE REQUIRED] Replace the previous matched "
        "p-tau217 sample size with the observed count shown above."
    )


# ## 7. Primary seven-protein panel: matched clinical associations
#
# The exact primary seven-protein panel is entered jointly into each clinical
# model in the matched cohort. Country is retained as a covariate. These models
# do not use the newly selected matched recurrent panel.

# ============================================================
# 7. MATCHED REGRESSION
# ============================================================

OUTCOMES = {
    "mmse_total": ("Global cognition", 1.0),
    "cog_craft_verb_delayed": ("Memory", 1.0),
    "udsfaq_total": ("Functionality", -1.0),
    "T.ADLQ": ("Deterioration", -1.0),
    "cdr_boxscore": ("Dementia severity", -1.0),
}


def run_matched_regressions(
    matched_df: pd.DataFrame,
    panel: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    missing_panel = [protein for protein in panel if protein not in matched_df.columns]
    if missing_panel:
        raise ValueError(f"Matched panel proteins missing from master: {missing_panel}")

    # Disease-direction map for bar colours.
    full_dep, _, _ = compute_training_dep(
        matched_df,
        panel,
        fold_label="FullMatched",
    )
    direction_map = dict(zip(full_dep["Protein"], full_dep["logFC"]))

    coefficient_frames = []
    model_rows = []

    for outcome_column, (outcome_label, sign) in OUTCOMES.items():
        if outcome_column not in matched_df.columns:
            print(f"[WARN] Missing outcome: {outcome_column}")
            continue

        columns = list(panel) + [outcome_column, "Country"]
        work = matched_df[columns].copy()
        for protein in panel:
            work[protein] = pd.to_numeric(work[protein], errors="coerce")
        work[outcome_column] = pd.to_numeric(
            work[outcome_column], errors="coerce"
        )
        work = work.dropna().copy()

        y_raw = sign * work[outcome_column].astype(float)
        y_standardized = (
            y_raw - y_raw.mean()
        ) / y_raw.std(ddof=0)

        protein_matrix = work[list(panel)].astype(float)
        protein_standardized = (
            protein_matrix - protein_matrix.mean(axis=0)
        ) / protein_matrix.std(axis=0, ddof=0)

        country_dummies = pd.get_dummies(
            work["Country"].astype(str),
            prefix="Country",
            drop_first=True,
            dtype=float,
        )

        design = pd.concat(
            [protein_standardized, country_dummies],
            axis=1,
        )
        design = sm.add_constant(design, has_constant="add")
        model = sm.OLS(y_standardized, design).fit()

        confidence = model.conf_int()
        coefficients = pd.DataFrame({
            "Outcome": outcome_label,
            "Variable": model.params.index,
            "Beta_std": model.params.values,
            "SE": model.bse.values,
            "t": model.tvalues.values,
            "P": model.pvalues.values,
            "CI_low": confidence.iloc[:, 0].values,
            "CI_high": confidence.iloc[:, 1].values,
        })
        coefficients["Protein_logFC_AD_vs_CN"] = coefficients[
            "Variable"
        ].map(direction_map)
        coefficient_frames.append(coefficients)

        r_squared = float(model.rsquared)
        f_squared = (
            r_squared / (1 - r_squared)
            if r_squared < 1 else np.inf
        )
        model_rows.append({
            "Outcome": outcome_label,
            "N": len(work),
            "R2": r_squared,
            "R2_adjusted": float(model.rsquared_adj),
            "F2": f_squared,
            "F": float(model.fvalue),
            "Model_P": float(model.f_pvalue),
            "RMSE_standardized": float(
                np.sqrt(np.mean(model.resid ** 2))
            ),
            "N_proteins": len(panel),
            "N_country_covariates": country_dummies.shape[1],
        })

        coefficients.to_csv(
            MATCH_REG_DIR
            / f"Regression_CN_AD_{outcome_label.replace(' ', '')}.csv",
            index=False,
        )

    coefficients_all = pd.concat(coefficient_frames, ignore_index=True)
    model_summary = pd.DataFrame(model_rows)
    direction_table = full_dep[[
        "Protein", "logFC", "P.Value", "adj.P.Val"
    ]].copy()

    coefficients_all.to_csv(
        MATCH_REG_DIR / "primary_panel_matched_regression_coefficients_all.csv",
        index=False,
    )
    model_summary.to_csv(
        MATCH_REG_DIR / "primary_panel_matched_regression_model_summary.csv",
        index=False,
    )
    direction_table.to_csv(
        MATCH_REG_DIR / "primary_panel_matched_disease_directions.csv",
        index=False,
    )
    return coefficients_all, model_summary, direction_table


matched = restore_matched_cohort_191()

required_regression_outputs = [
    MATCH_REG_DIR / "primary_panel_matched_regression_coefficients_all.csv",
    MATCH_REG_DIR / "primary_panel_matched_regression_model_summary.csv",
    MATCH_REG_DIR / "primary_panel_matched_disease_directions.csv",
]

if FORCE_RECOMPUTE_REGRESSION or not all(
    path.exists() for path in required_regression_outputs
):
    matched_regression_results = run_matched_regressions(
        matched,
        PRIMARY_PANEL,
    )
else:
    print("Loading cached matched regression results.")
    matched_regression_results = (
        pd.read_csv(
            MATCH_REG_DIR / "primary_panel_matched_regression_coefficients_all.csv"
        ),
        pd.read_csv(
            MATCH_REG_DIR / "primary_panel_matched_regression_model_summary.csv"
        ),
        pd.read_csv(
            MATCH_REG_DIR / "primary_panel_matched_disease_directions.csv"
        ),
    )

regression_coefficients, regression_models, disease_directions = (
    matched_regression_results
)
print(regression_models.to_string(index=False))


# ## 8. Load validated panels a–b and assemble the complete figure

# ============================================================
# 8. FIGURE HELPERS
# ============================================================

def normalize_model_names(df: pd.DataFrame) -> pd.DataFrame:
    output = df.copy()
    if "Model" not in output.columns:
        return output
    replacements = {
        "SomaScan": "Proteins",
        "Proteins + APOE": "Proteins + APOE ε4",
        "Proteins + APOE4": "Proteins + APOE ε4",
        "Panel": "Panel",
        "pTau217": "p-tau217",
        "pTau": "p-tau217",
        "Panel+pTau217": "Combined",
        "Panel + p-tau217": "Combined",
        "Proteins + p-tau217": "Combined",
        "Matched proteins": "Matched nested",
        "PrimaryPanel": "Primary panel",
    }
    output["Model"] = output["Model"].replace(replacements)
    return output


def ensure_roc_auc_column(roc: pd.DataFrame) -> pd.DataFrame:
    output = roc.copy()
    if "AUC" in output.columns:
        return output
    rows = []
    group_columns = [
        column for column in ["Model", "Fold"] if column in output.columns
    ]
    for keys, group in output.groupby(group_columns):
        group = group.sort_values("FPR").copy()
        auc_value = sklearn_auc(group["FPR"], group["TPR"])
        group["AUC"] = auc_value
        rows.append(group)
    return pd.concat(rows, ignore_index=True)


def roc_auc_summary(roc: pd.DataFrame) -> pd.DataFrame:
    roc = ensure_roc_auc_column(roc)
    if "Fold" in roc.columns:
        per_fold = (
            roc.groupby(["Model", "Fold"], as_index=False)["AUC"].first()
        )
    else:
        per_fold = roc[["Model", "AUC"]].drop_duplicates()
    return per_fold.groupby("Model").agg(
        Mean_AUC=("AUC", "mean"),
        SD_AUC=("AUC", "std"),
        N_folds=("AUC", "count"),
    ).reset_index()


def plot_mean_roc(
    ax: plt.Axes,
    roc: pd.DataFrame,
    model_order: Sequence[str],
    colours: Dict[str, str],
    labels: Optional[Dict[str, str]] = None,
    title: str = "",
    legend_location: str = "lower right",
) -> None:
    roc = ensure_roc_auc_column(normalize_model_names(roc))
    labels = labels or {}
    mean_fpr = np.linspace(0, 1, 200)
    legend_handles = []
    legend_labels = []

    for model in model_order:
        subset = roc[roc["Model"] == model].copy()
        if subset.empty:
            continue
        fold_values = (
            sorted(subset["Fold"].unique())
            if "Fold" in subset.columns else [None]
        )
        tprs = []
        aucs = []
        for fold in fold_values:
            fold_data = (
                subset[subset["Fold"] == fold]
                if fold is not None else subset
            ).sort_values("FPR")
            interpolated = np.interp(
                mean_fpr,
                fold_data["FPR"].astype(float),
                fold_data["TPR"].astype(float),
            )
            interpolated[0] = 0
            interpolated[-1] = 1
            tprs.append(interpolated)
            aucs.append(float(fold_data["AUC"].iloc[0]))

        tpr_matrix = np.vstack(tprs)
        mean_tpr = tpr_matrix.mean(axis=0)
        sd_tpr = (
            tpr_matrix.std(axis=0, ddof=1)
            if len(tprs) > 1 else np.zeros_like(mean_tpr)
        )
        colour = colours[model]
        line = ax.plot(
            mean_fpr,
            mean_tpr,
            lw=1.25,
            color=colour,
        )[0]
        if len(tprs) > 1:
            ax.fill_between(
                mean_fpr,
                np.clip(mean_tpr - sd_tpr, 0, 1),
                np.clip(mean_tpr + sd_tpr, 0, 1),
                color=colour,
                alpha=0.15,
                linewidth=0,
            )
        legend_handles.append(line)
        display = labels.get(model, model)
        legend_labels.append(
            f"{display}: {np.mean(aucs):.2f} ± {np.std(aucs, ddof=1):.2f}"
        )

    ax.plot([0, 1], [0, 1], "--", color="0.55", lw=0.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_xlabel("1 - Specificity", fontsize=7)
    ax.set_ylabel("Sensitivity", fontsize=7)
    ax.set_title(title, fontsize=7.8, pad=2.5)
    ax.tick_params(labelsize=6, length=2, width=0.5)
    ax.legend(
        legend_handles,
        legend_labels,
        loc=legend_location,
        fontsize=5.5,
        frameon=True,
        framealpha=0.95,
        borderpad=0.4,
        handlelength=2.2,
        handletextpad=0.5,
    )


def standardize_permutation_columns(df: pd.DataFrame) -> pd.DataFrame:
    output = df.copy()
    rename = {}
    for column in output.columns:
        lower = column.lower()
        if lower in {"proteins", "protein", "feature", "variable"}:
            rename[column] = "Protein"
        elif lower in {"permutation", "importance", "permutation_importance"}:
            rename[column] = "Importance"
        elif lower in {"permutation_sd", "importance_sd", "sd"}:
            rename[column] = "Importance_SD"
    return output.rename(columns=rename)


def plot_permutation(
    ax: plt.Axes,
    permutation_df: pd.DataFrame,
    title: str,
    colour: str,
    top_n: int = 8,
) -> None:
    data = standardize_permutation_columns(permutation_df)
    if "Importance_SD" not in data.columns:
        data["Importance_SD"] = np.nan
    data = (
        data.sort_values("Importance", ascending=False)
        .head(top_n)
        .sort_values("Importance")
    )
    positions = np.arange(len(data))
    ax.barh(
        positions,
        data["Importance"],
        xerr=(
            data["Importance_SD"]
            if data["Importance_SD"].notna().any()
            else None
        ),
        color=colour,
        ecolor="0.15",
        linewidth=0,
        error_kw={"elinewidth": 0.7},
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(data["Protein"], fontsize=6)
    ax.set_xlabel("Feature importance", fontsize=7)
    ax.set_title(title, fontsize=7.6, pad=2.5)
    ax.tick_params(axis="x", labelsize=6, length=2, width=0.5)
    ax.tick_params(axis="y", labelsize=5.8, length=0, pad=1.5)
    ax.margins(y=0.04)


def comparison_alias(name: str) -> str:
    replacements = {
        "Soma": "Proteins",
        "SomaScan": "Proteins",
        "Panel": "Panel",
        "pTau": "p-tau217",
        "pTau217": "p-tau217",
        "Combo": "Combined",
        "Combined": "Combined",
        "Proteins": "Proteins",
        "PrimaryPanel": "Primary panel",
        "Primary panel": "Primary panel",
    }
    return replacements.get(name, name)


def plot_delong_heatmap(
    ax: plt.Axes,
    delong: pd.DataFrame,
    order: Sequence[str],
    title: str,
) -> None:
    matrix = np.full((len(order), len(order)), np.nan)
    annotations = [["" for _ in order] for _ in order]
    index_lookup = {label: index for index, label in enumerate(order)}

    for _, row in delong.iterrows():
        comparison = str(row["Comparison"])
        parts = re.split(r"_vs_| vs ", comparison)
        if len(parts) < 2:
            continue
        first = comparison_alias(parts[0])
        second = comparison_alias(parts[1])
        if first not in index_lookup or second not in index_lookup:
            continue
        i = index_lookup[first]
        j = index_lookup[second]
        p_value = float(row["P"])
        matrix[i, j] = p_value
        matrix[j, i] = p_value
        annotations[i][j] = format_p(p_value)
        annotations[j][i] = format_p(p_value)

    masked = np.ma.masked_invalid(matrix)
    colour_map = plt.get_cmap("viridis_r").copy()
    colour_map.set_bad("white")
    image = ax.imshow(
        masked,
        cmap=colour_map,
        norm=Normalize(vmin=0, vmax=0.05, clip=True),
    )
    for row_index in range(len(order)):
        for column_index in range(len(order)):
            if annotations[row_index][column_index]:
                ax.text(
                    column_index,
                    row_index,
                    annotations[row_index][column_index],
                    ha="center",
                    va="center",
                    fontsize=4.7,
                )
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5.5)
    ax.set_yticklabels(order, fontsize=5.5)
    ax.set_title(title, fontsize=7.5, pad=2.5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    colour_bar = plt.colorbar(
        ScalarMappable(
            norm=Normalize(vmin=0, vmax=0.05),
            cmap=colour_map,
        ),
        ax=ax,
        fraction=0.065,
        pad=0.04,
    )
    colour_bar.ax.tick_params(labelsize=4.8, length=1.8, pad=1.5)
    colour_bar.set_label("P value", fontsize=5.5, labelpad=2)


def plot_regression_panel(
    ax: plt.Axes,
    outcome: str,
    coefficients: pd.DataFrame,
    model_summary: pd.DataFrame,
    panel: Sequence[str],
) -> None:
    data = coefficients[
        (coefficients["Outcome"] == outcome)
        & (coefficients["Variable"].isin(panel))
    ].copy()
    data = data.set_index("Variable").reindex(panel).reset_index()
    data = data.sort_values("Beta_std")

    colours = []
    for _, row in data.iterrows():
        higher_in_ad = row["Protein_logFC_AD_vs_CN"] > 0
        significant = row["P"] < 0.05
        if higher_in_ad:
            colours.append(
                COL_HIGHER_AD if significant else COL_HIGHER_AD_LIGHT
            )
        else:
            colours.append(
                COL_LOWER_AD if significant else COL_LOWER_AD_LIGHT
            )

    positions = np.arange(len(data))
    ax.barh(
        positions,
        data["Beta_std"],
        color=colours,
        linewidth=0,
    )
    ax.axvline(0, color="0.2", linestyle="--", lw=0.6)
    ax.set_yticks(positions)
    ax.set_yticklabels(data["Variable"], fontsize=5.5)
    ax.set_xlabel("Standardized beta", fontsize=6)
    ax.tick_params(axis="x", labelsize=5.5, length=2, width=0.5)
    ax.tick_params(axis="y", length=0)

    summary = model_summary[
        model_summary["Outcome"] == outcome
    ].iloc[0]
    p_text = (
        f"P < 1×10$^{{-15}}$"
        if summary["Model_P"] < 1e-15
        else f"P = {summary['Model_P']:.1e}"
    )
    ax.set_title(
        f"{outcome}\n"
        f"R² = {summary['R2']:.2f}, f² = {summary['F2']:.2f}\n"
        f"F = {summary['F']:.2f}, {p_text}, n = {int(summary['N'])}",
        fontsize=5.3,
        pad=1.5,
    )
    maximum = np.nanmax(np.abs(data["Beta_std"]))
    maximum = maximum if np.isfinite(maximum) and maximum > 0 else 1
    ax.set_xlim(-maximum * 1.15, maximum * 1.15)

# ============================================================
# 8B. LOAD SOURCE TABLES
# ============================================================

apoe_roc = normalize_model_names(
    pd.read_csv(APOE_DIR / "roc_curves.csv")
)
apoe_permutation = pd.read_csv(
    APOE_DIR / "Final_Panels" / "permutation_importance_freq08.csv"
)
apoe_delong = pd.read_csv(APOE_DIR / "delong_results.csv")

fixed_ptau_roc = normalize_model_names(
    ensure_roc_auc_column(
        pd.read_csv(FIXED_PTAU_DIR / "roc_fold_curves.csv")
    )
)
fixed_ptau_delong = pd.read_csv(
    FIXED_PTAU_DIR / "delong_results.csv"
)

matched_roc = normalize_model_names(matched_fixed_comparison_roc)
matched_permutation = matched_results["permutation"].copy()
matched_ptau_roc = normalize_model_names(matched_ptau_results["roc"])
matched_ptau_delong = matched_ptau_results["delong"].copy()

matched_ptau_cohort = matched_ptau_results["cohort"].iloc[0]
matched_ptau_n_cn = int(matched_ptau_cohort["N_CN"])
matched_ptau_n_ad = int(matched_ptau_cohort["N_AD"])

matched_panel_title = "Matched nested and fixed panel\n(191 CN; 191 AD)"
matched_ptau_panel_title = (
    "Primary panel and p-tau217 (matched)\n"
    f"({matched_ptau_n_cn} CN; {matched_ptau_n_ad} AD)"
)

# ============================================================
# 9. ASSEMBLE EXTENDED DATA FIGURE 10 — FIG. 4 FORMAT
# ============================================================

# Exact page dimensions of Figure_4(3).pdf.
FIGURE_WIDTH_IN = 510.24 / 72.0
FIGURE_HEIGHT_IN = 467.667 / 72.0

fig = plt.figure(
    figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN),
    dpi=300,
)

# ------------------------------------------------------------------
# FIGURE 4-ALIGNED GEOMETRY
# ------------------------------------------------------------------
# Coordinates are expressed in figure fractions and were calibrated against
# Figure_4(4).pdf (510.24 x 467.667 pt). Explicit axes positions are used
# because a three-row GridSpec made the auxiliary panels too tall, shifted the
# first row downward and compressed the clinical row heading.
#
# Main Figure 4 reference geometry:
#   - ROC panels: ~113.4 x 113.4 pt
#   - top ROC left edges: ~25 pt and ~287 pt
#   - compact auxiliary panels centred vertically beside the ROC panels
#   - clinical panels: narrower, with larger horizontal gaps
#
# Top-row ROC panels (a and b)
ROC_LEFT = 0.049
ROC_RIGHT = 0.563
ROC_WIDTH = 0.222
ROC_HEIGHT = 0.242
TOP_ROC_BOTTOM = 0.720

# Second-row ROC panels (c and d)
SECOND_ROC_BOTTOM = 0.382

# Feature-importance panels are deliberately shorter than the ROC panels,
# matching the visual hierarchy of Figure 4.
IMPORTANCE_LEFT = 0.347
IMPORTANCE_WIDTH = 0.135
IMPORTANCE_HEIGHT = 0.205
TOP_IMPORTANCE_BOTTOM = 0.739
SECOND_IMPORTANCE_BOTTOM = 0.401

# DeLong heatmaps are compact in Figure 4 and should not compete visually
# with the ROC curves.
HEATMAP_LEFT = 0.855
HEATMAP_WIDTH = 0.100
HEATMAP_HEIGHT = 0.108
TOP_HEATMAP_BOTTOM = 0.808
SECOND_HEATMAP_BOTTOM = 0.470

# Clinical row: five narrower panels with Figure 4-like inter-panel spacing.
CLINICAL_LEFTS = [0.055, 0.252, 0.449, 0.646, 0.843]
CLINICAL_WIDTH = 0.130
CLINICAL_BOTTOM = 0.075
CLINICAL_HEIGHT = 0.178

ax_a_roc = fig.add_axes([
    ROC_LEFT,
    TOP_ROC_BOTTOM,
    ROC_WIDTH,
    ROC_HEIGHT,
])
ax_a_perm = fig.add_axes([
    IMPORTANCE_LEFT,
    TOP_IMPORTANCE_BOTTOM,
    IMPORTANCE_WIDTH,
    IMPORTANCE_HEIGHT,
])

ax_b_roc = fig.add_axes([
    ROC_RIGHT,
    TOP_ROC_BOTTOM,
    ROC_WIDTH,
    ROC_HEIGHT,
])
ax_b_delong = fig.add_axes([
    HEATMAP_LEFT,
    TOP_HEATMAP_BOTTOM,
    HEATMAP_WIDTH,
    HEATMAP_HEIGHT,
])

ax_c_roc = fig.add_axes([
    ROC_LEFT,
    SECOND_ROC_BOTTOM,
    ROC_WIDTH,
    ROC_HEIGHT,
])
ax_c_perm = fig.add_axes([
    IMPORTANCE_LEFT,
    SECOND_IMPORTANCE_BOTTOM,
    IMPORTANCE_WIDTH,
    IMPORTANCE_HEIGHT,
])

ax_d_roc = fig.add_axes([
    ROC_RIGHT,
    SECOND_ROC_BOTTOM,
    ROC_WIDTH,
    ROC_HEIGHT,
])
ax_d_delong = fig.add_axes([
    HEATMAP_LEFT,
    SECOND_HEATMAP_BOTTOM,
    HEATMAP_WIDTH,
    HEATMAP_HEIGHT,
])

axes_e = [
    fig.add_axes([
        left,
        CLINICAL_BOTTOM,
        CLINICAL_WIDTH,
        CLINICAL_HEIGHT,
    ])
    for left in CLINICAL_LEFTS
]

# Panel a
plot_mean_roc(
    ax_a_roc,
    apoe_roc,
    model_order=["Proteins", "Proteins + APOE ε4"],
    colours={
        "Proteins": COL_PROTEIN,
        "Proteins + APOE ε4": COL_APOE,
    },
    labels={
        "Proteins": "Proteins",
        "Proteins + APOE ε4": "Proteins + APOE ε4",
    },
    title="Proteins ± APOE ε4",
)
plot_permutation(
    ax_a_perm,
    apoe_permutation,
    title="Feature importance",
    colour=COL_APOE,
    top_n=8,
)
if not apoe_delong.empty:
    apoe_p = float(apoe_delong["P"].iloc[0])
    ax_a_roc.text(
        0.98, 0.25,
        f"DeLong P = {apoe_p:.3f}",
        ha="right", va="center",
        transform=ax_a_roc.transAxes,
        fontsize=5.3,
    )

# Panel b
plot_mean_roc(
    ax_b_roc,
    fixed_ptau_roc,
    model_order=["Panel", "p-tau217", "Combined"],
    colours={
        "Panel": COL_PROTEIN,
        "p-tau217": COL_PTAU,
        "Combined": COL_COMBINED,
    },
    labels={
        "Panel": "Panel",
        "p-tau217": "p-tau217",
        "Combined": "Combined",
    },
    title="Fixed panel and p-tau217",
)
plot_delong_heatmap(
    ax_b_delong,
    fixed_ptau_delong,
    order=["Panel", "p-tau217", "Combined"],
    title="DeLong P values",
)

# Panel c — exact 191 CN + 191 AD
plot_mean_roc(
    ax_c_roc,
    matched_roc,
    model_order=["Matched nested", "Primary fixed panel"],
    colours={
        "Matched nested": COL_PROTEIN,
        "Primary fixed panel": COL_FIXED_PANEL,
    },
    labels={
        "Matched nested": "Matched nested",
        "Primary fixed panel": "Primary fixed",
    },
    title=matched_panel_title,
)
plot_permutation(
    ax_c_perm,
    matched_permutation,
    title="Matched recurrent proteins",
    colour=COL_PROTEIN,
    top_n=len(matched_panel),
)

# Panel d — observed p-tau217 complete cases within the 191/191 match
plot_mean_roc(
    ax_d_roc,
    matched_ptau_roc,
    model_order=["Primary panel", "p-tau217", "Combined"],
    colours={
        "Primary panel": COL_FIXED_PANEL,
        "p-tau217": COL_PTAU,
        "Combined": COL_COMBINED,
    },
    labels={
        "Primary panel": "Primary panel",
        "p-tau217": "p-tau217",
        "Combined": "Combined",
    },
    title=matched_ptau_panel_title,
)
plot_delong_heatmap(
    ax_d_delong,
    matched_ptau_delong,
    order=["Primary panel", "p-tau217", "Combined"],
    title="DeLong P values",
)

# Panel e
outcome_order = [
    "Global cognition",
    "Memory",
    "Functionality",
    "Deterioration",
    "Dementia severity",
]
for axis, outcome in zip(axes_e, outcome_order):
    plot_regression_panel(
        axis,
        outcome,
        regression_coefficients,
        regression_models,
        PRIMARY_PANEL,
    )

# Exact axis geometry is already defined above. Keep ROC data scales square
# without allowing Matplotlib to recenter them inside oversized GridSpec cells.
for axis in [ax_a_roc, ax_b_roc, ax_c_roc, ax_d_roc]:
    axis.set_aspect("equal", adjustable="box")

# Fixed panel letters mirror the locations used in Figure 4.
PANEL_LETTER_STYLE = {
    "fontsize": 9.5,
    "fontweight": "bold",
    "ha": "left",
    "va": "top",
}
fig.text(0.006, 0.986, "a", **PANEL_LETTER_STYLE)
fig.text(0.520, 0.986, "b", **PANEL_LETTER_STYLE)
fig.text(0.006, 0.648, "c", **PANEL_LETTER_STYLE)
fig.text(0.520, 0.648, "d", **PANEL_LETTER_STYLE)

# Reserve a dedicated title band above the clinical axes. In v15 the panel-e
# heading and the five outcome titles occupied the same vertical space.
fig.text(
    0.006,
    0.326,
    "e",
    fontsize=9.5,
    fontweight="bold",
    ha="left",
    va="top",
)
fig.text(
    0.041,
    0.326,
    "Primary seven-protein panel in the matched cohort",
    fontsize=7.2,
    ha="left",
    va="top",
)

for axis in [
    ax_a_roc, ax_a_perm, ax_b_roc, ax_b_delong,
    ax_c_roc, ax_c_perm, ax_d_roc, ax_d_delong,
    *axes_e,
]:
    for spine in axis.spines.values():
        spine.set_linewidth(0.6)

# Export the final panel geometry in PDF points for reproducibility and for
# direct comparison with the Illustrator Figure 4 canvas.
layout_axes = {
    "a_ROC": ax_a_roc,
    "a_importance": ax_a_perm,
    "b_ROC": ax_b_roc,
    "b_DeLong": ax_b_delong,
    "c_ROC": ax_c_roc,
    "c_importance": ax_c_perm,
    "d_ROC": ax_d_roc,
    "d_DeLong": ax_d_delong,
}
layout_axes.update({
    f"e_{index + 1}_{outcome.replace(' ', '_')}": axis
    for index, (axis, outcome) in enumerate(zip(axes_e, outcome_order))
})
layout_rows = []
for panel_name, axis in layout_axes.items():
    box = axis.get_position()
    layout_rows.append({
        "Panel": panel_name,
        "Left_pt": box.x0 * 510.24,
        "Bottom_pt": box.y0 * 467.667,
        "Width_pt": box.width * 510.24,
        "Height_pt": box.height * 467.667,
        "Top_from_PDF_pt": (1 - box.y1) * 467.667,
    })
pd.DataFrame(layout_rows).to_csv(
    FIGURE_DIR / "figure_geometry_v17_points.csv",
    index=False,
)

pdf_path = (
    FIGURE_DIR
    / "Extended_Data_Fig10_FINAL_v17.pdf"
)
png_path = (
    FIGURE_DIR
    / "Extended_Data_Fig10_FINAL_v17.png"
)
tiff_path = (
    FIGURE_DIR
    / "Extended_Data_Fig10_FINAL_v17.tiff"
)

# Fixed export box. Do not use bbox_inches="tight": it changes the final page
# dimensions when panel letters or headings extend beyond an axes object.
fig.savefig(
    pdf_path,
    bbox_inches=None,
    pad_inches=0,
)
fig.savefig(
    png_path,
    bbox_inches=None,
    pad_inches=0,
    dpi=600,
)
fig.savefig(
    tiff_path,
    bbox_inches=None,
    pad_inches=0,
    dpi=600,
    pil_kwargs={"compression": "tiff_lzw"},
)
# Confirm that all exports exist and are non-empty.
exported_files = [pdf_path, png_path, tiff_path]
missing_exports = [
    path for path in exported_files
    if (not path.exists()) or path.stat().st_size == 0
]
if missing_exports:
    raise RuntimeError(
        "Figure export did not complete. Missing or empty files: "
        + ", ".join(str(path) for path in missing_exports)
    )

success_file = OUTPUT_ROOT / "EXPORT_SUCCESSFUL.txt"
success_file.write_text(
    "Extended Data Fig. 10 v17 exported successfully.\n\n"
    + "\n".join(
        f"figures/{path.name} | {path.stat().st_size:,} bytes"
        for path in exported_files
    )
    + "\n",
    encoding="utf-8",
)

plt.show()

print("\nFIGURES EXPORTED SUCCESSFULLY")
for path in exported_files:
    print(f"  {path} ({path.stat().st_size:,} bytes)")
print("Status file:", success_file)
print(
    "Fixed figure size (inches):",
    FIGURE_WIDTH_IN,
    "×",
    FIGURE_HEIGHT_IN,
)


# ============================================================
# 10. SOURCE DATA, AUDIT AND MANUSCRIPT-CONSISTENCY CHECKS
# ============================================================

apoe_auc_summary = roc_auc_summary(apoe_roc)
fixed_ptau_auc_summary = roc_auc_summary(fixed_ptau_roc)
matched_auc_summary = roc_auc_summary(matched_roc)
primary_fixed_auc_summary = roc_auc_summary(primary_fixed_results["roc"])
matched_ptau_auc_summary = roc_auc_summary(matched_ptau_roc)

counts_matched = matched["SampleGroup"].value_counts()
counts_ptau = matched_ptau_results["cohort"].iloc[0]

canonical_model_audit_file = (
    MATCH_RESULTS_DIR / "canonical_191_model_audit.csv"
)
canonical_model_audit = pd.read_csv(
    canonical_model_audit_file
)
if int(canonical_model_audit.iloc[0]["N_total"]) != 382:
    raise RuntimeError(
        "Panel c model audit does not contain 382 participants."
    )
if int(canonical_model_audit.iloc[0]["N_unique_outer_test_IDs"]) != 382:
    raise RuntimeError(
        "Panel c outer-fold predictions do not cover all 382 participants."
    )

consistency_checks = pd.DataFrame([
    {
        "Check": "Matched CN",
        "Observed": int(counts_matched.get("CN", 0)),
        "Manuscript": TARGET_CN,
    },
    {
        "Check": "Matched AD",
        "Observed": int(counts_matched.get("AD", 0)),
        "Manuscript": TARGET_AD,
    },
    {
        "Check": "Matched p-tau CN",
        "Observed": int(counts_ptau["N_CN"]),
        "Manuscript": 158,
    },
    {
        "Check": "Matched p-tau AD",
        "Observed": int(counts_ptau["N_AD"]),
        "Manuscript": 157,
    },
])
consistency_checks["Matches"] = (
    consistency_checks["Observed"] == consistency_checks["Manuscript"]
)
consistency_checks["Interpretation"] = np.where(
    consistency_checks["Matches"],
    "Consistent",
    "Previous manuscript count requires correction",
)

source_xlsx = SOURCE_DATA_DIR / "Extended_Data_Fig10_Source_Data.xlsx"
with pd.ExcelWriter(source_xlsx, engine="openpyxl") as writer:
    apoe_auc_summary.to_excel(
        writer, sheet_name="a_APOE_AUC", index=False
    )
    apoe_delong.to_excel(
        writer, sheet_name="a_APOE_DeLong", index=False
    )
    standardize_permutation_columns(apoe_permutation).to_excel(
        writer, sheet_name="a_APOE_importance", index=False
    )
    fixed_ptau_auc_summary.to_excel(
        writer, sheet_name="b_fixed_panel_AUC", index=False
    )
    fixed_ptau_delong.to_excel(
        writer, sheet_name="b_fixed_panel_DeLong", index=False
    )
    matched_auc_summary.to_excel(
        writer, sheet_name="c_matched_AUC", index=False
    )
    standardize_permutation_columns(matched_permutation).to_excel(
        writer, sheet_name="c_matched_importance", index=False
    )
    primary_fixed_auc_summary.to_excel(
        writer, sheet_name="c_primary_fixed_AUC", index=False
    )
    primary_fixed_results["metrics"].to_excel(
        writer, sheet_name="c_primary_fixed_metrics", index=False
    )
    matched_panel_overlap.to_excel(
        writer, sheet_name="c_panel_overlap", index=False
    )
    matched_ptau_auc_summary.to_excel(
        writer, sheet_name="d_matched_ptau_AUC", index=False
    )
    matched_ptau_delong.to_excel(
        writer, sheet_name="d_matched_ptau_DeLong", index=False
    )
    regression_coefficients.to_excel(
        writer, sheet_name="e_coefficients", index=False
    )
    regression_models.to_excel(
        writer, sheet_name="e_model_summary", index=False
    )
    matching_audit.to_excel(
        writer, sheet_name="matching_balance", index=False
    )
    if not historical_reconstruction_audit.empty:
        historical_reconstruction_audit.to_excel(
            writer, sheet_name="historical_cohort_audit", index=False
        )
    consistency_checks.to_excel(
        writer, sheet_name="manuscript_checks", index=False
    )

consistency_note_file = SOURCE_DATA_DIR / (
    "MANUSCRIPT_CONSISTENCY_CONFIRMED.txt"
    if consistency_checks["Matches"].all()
    else "MANUSCRIPT_UPDATE_REQUIRED.txt"
)
with open(consistency_note_file, "w", encoding="utf-8") as handle:
    if consistency_checks["Matches"].all():
        handle.write(
            "MATCHED p-tau217 SAMPLE-SIZE AUDIT\n"
            "===================================\n\n"
            f"Observed complete-case subset: CN n={int(counts_ptau['N_CN'])}; "
            f"AD n={int(counts_ptau['N_AD'])}.\n"
            "These counts match the submitted manuscript/figure legend.\n"
            "No p-tau217 value was imputed and no participant was replaced.\n"
        )
    else:
        handle.write(
            "MATCHED p-tau217 SAMPLE-SIZE AUDIT\n"
            "===================================\n\n"
            f"Observed reproducible complete-case subset: CN n={int(counts_ptau['N_CN'])}; "
            f"AD n={int(counts_ptau['N_AD'])}.\n"
            "Submitted manuscript counts: CN n=158; AD n=157.\n"
            "Review the discrepancy before submission. No value was imputed and no "
            "participant was replaced to force a count.\n"
        )

interpretation_file = (
    SOURCE_DATA_DIR / "MATCHED_PANEL_INTERPRETATION.txt"
)
with open(interpretation_file, "w", encoding="utf-8") as handle:
    handle.write(
        "MATCHED NESTED VERSUS PRIMARY FIXED PANEL\n"
        "==========================================\n\n"
        "The matched nested model reruns differential-abundance filtering "
        "and elastic-net selection inside the matched cohort and may select "
        "a different recurrent protein set.\n\n"
        "The primary fixed-panel analysis retains SPC25, CPLX2, TCP11L1, "
        "ACHE, ODC1, SPON1 and RTN4RL1 without feature reselection. "
        "Hyperparameter tuning remains confined to training folds.\n\n"
        f"Shared proteins relative to the primary panel: "
        f"{len(shared_proteins)}/{len(PRIMARY_PANEL)} — "
        f"{', '.join(shared_proteins)}.\n\n"
        "Both analyses are internal sensitivity analyses and do not "
        "constitute independent validation.\n"
    )

manifest = {
    "project": "ReDLat plasma proteomics reproducibility package",
    "master_file": "outputs/ML/private/derived/ReDLat_ML_gene_master_RAW.csv",
    "matched_source": matched_source,
    "selected_match_ids_file": "outputs/ML/private/Result_matching_rebuild_v5/matched_ids_SELECTED.csv",
    "matched_reselected_panel": matched_panel,
    "primary_fixed_panel": PRIMARY_PANEL,
    "shared_primary_matched_proteins": shared_proteins,
    "settings": {
        "random_state": RANDOM_STATE,
        "matching_seed": MATCHING_SEED,
        "exact_country": False,
        "matching_performed_in_this_script": False,
        "primary_panel_feature_reselection": False,
        "matched_nested_feature_reselection": True,
        "match_caliper_sd": MATCH_CALIPER_SD,
        "outer_folds": OUTER_FOLDS,
        "inner_folds": INNER_FOLDS,
        "dep_fdr": DEP_FDR,
        "inner_stability": INNER_STABILITY,
        "final_panel_threshold": FINAL_PANEL_THRESHOLD,
        "final_panel_padding": False,
        "svm_search": "BayesSearchCV" if HAS_SKOPT else "GridSearchCV",
    },
    "outputs": {
        "pdf": f"figures/{pdf_path.name}",
        "png": f"figures/{png_path.name}",
        "tiff": f"figures/{tiff_path.name}",
        "source_data": f"source_data/{source_xlsx.name}",
        "manuscript_consistency_note": f"source_data/{consistency_note_file.name}",
        "matched_panel_interpretation": f"source_data/{interpretation_file.name}",
    },
}
with open(
    SOURCE_DATA_DIR / "Extended_Data_Fig10_manifest.json",
    "w",
    encoding="utf-8",
) as handle:
    json.dump(manifest, handle, indent=2, ensure_ascii=False)

print("\nAUC summaries")
print("\nPanel a")
print(apoe_auc_summary.to_string(index=False))
print("\nPanel b")
print(fixed_ptau_auc_summary.to_string(index=False))
print("\nPanel c")
print(matched_auc_summary.to_string(index=False))
print("\nPanel d")
print(matched_ptau_auc_summary.to_string(index=False))

print("\nManuscript consistency checks")
print(consistency_checks.to_string(index=False))

if not consistency_checks["Matches"].all():
    print(
        "\nWARNING: At least one cohort count differs from the manuscript. "
        "Do not update the manuscript or submit the figure until the selected "
        "matched-ID version has been confirmed."
    )

print("\nSource data workbook:", source_xlsx)
