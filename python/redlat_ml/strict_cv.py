from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, Any
import json
import os
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from scipy.stats import loguniform, norm
from sklearn.calibration import CalibratedClassifierCV
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, confusion_matrix, f1_score, precision_score,
    recall_score, roc_auc_score, roc_curve,
)
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

try:
    from skopt import BayesSearchCV
    from skopt.space import Categorical, Real
    HAS_SKOPT = True
except Exception:
    HAS_SKOPT = False


META_COLUMNS = {
    "SampleId", "SampleGroup", "Site", "Country", "Sex", "Age", "Education",
    "ApoE", "APOE_group", "APOE4_carrier", "cdr_global", "cdr_boxscore",
    "mmse_total", "udsfaq_total", "cog_benson", "cog_tmt_a", "cog_tmt_b",
    "cog_craft_verb_delayed", "NPI", "Mini.SEA", "T.ADLQ", "p.tau217",
    "p-tau217", "pTau217", "ptau217", "p.tau181", "NfL", "ratio.AB42.40",
    "GFAP_1", "distance", "weights", "subclass", "SampleGroup_bin",
    "pair_id", "match_pair_id", "propensity_score", "propensity_logit",
}


@dataclass(frozen=True)
class StrictCVSettings:
    random_state: int = 42
    inner_folds: int = 5
    stability_threshold: float = 1.0
    coefficient_threshold: float = 5e-3
    fallback_top_n: int = 10
    max_missing_fraction: float = 0.20
    min_variance: float = 1e-12
    selector_search_iter: int = 30
    svm_search_iter: int = 50
    permutation_repeats: int = 50
    max_iter: int = 50000
    tolerance: float = 1e-3
    n_jobs: int = -1

    @classmethod
    def from_env(cls) -> "StrictCVSettings":
        def geti(name: str, default: int) -> int:
            return int(os.getenv(name, str(default)))
        def getf(name: str, default: float) -> float:
            return float(os.getenv(name, str(default)))
        return cls(
            random_state=geti("REDLAT_ML_RANDOM_STATE", 42),
            inner_folds=geti("REDLAT_ML_INNER_FOLDS", 5),
            stability_threshold=getf("REDLAT_ML_INNER_STABILITY", 1.0),
            coefficient_threshold=getf("REDLAT_ML_COEFFICIENT_THRESHOLD", 5e-3),
            fallback_top_n=geti("REDLAT_ML_FALLBACK_TOP_N", 10),
            max_missing_fraction=getf("REDLAT_ML_MAX_MISSING_FRACTION", 0.20),
            min_variance=getf("REDLAT_ML_MIN_VARIANCE", 1e-12),
            selector_search_iter=geti("REDLAT_ML_SELECTOR_SEARCH_ITER", 30),
            svm_search_iter=geti("REDLAT_ML_SVM_SEARCH_ITER", 50),
            permutation_repeats=geti("REDLAT_ML_PERMUTATION_REPEATS", 50),
            max_iter=geti("REDLAT_ML_MAX_ITER", 50000),
            tolerance=getf("REDLAT_ML_TOLERANCE", 1e-3),
            n_jobs=geti("REDLAT_ML_N_JOBS", -1),
        )


def numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(pd.to_numeric, errors="coerce")


def resolve_column(frame: pd.DataFrame, candidates: Sequence[str]) -> str:
    direct = {str(c): str(c) for c in frame.columns}
    lower = {str(c).lower(): str(c) for c in frame.columns}
    for candidate in candidates:
        if candidate in direct:
            return direct[candidate]
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"None of the expected columns were found: {list(candidates)}")


def derive_apoe4_carrier(frame: pd.DataFrame) -> tuple[pd.Series, str]:
    for candidate in ("APOE4_carrier", "APOE_group"):
        if candidate in frame.columns:
            raw = frame[candidate]
            numeric = pd.to_numeric(raw, errors="coerce")
            if numeric.notna().any():
                return numeric.where(numeric.isin([0, 1])), candidate
            text = raw.astype(str).str.strip().str.lower()
            mapped = text.map({
                "carrier": 1, "e4 carrier": 1, "ε4 carrier": 1, "yes": 1,
                "true": 1, "non-carrier": 0, "noncarrier": 0, "no": 0,
                "false": 0,
            })
            if mapped.notna().any():
                return mapped, candidate
    genotype_col = resolve_column(frame, ("ApoE", "APOE", "apoe"))
    genotype = frame[genotype_col].astype(str).str.lower().str.replace("ε", "e", regex=False)
    genotype = genotype.str.replace(r"[^0-9e/]", "", regex=True)
    carrier = genotype.str.contains("4", regex=False).astype(float)
    carrier[frame[genotype_col].isna() | frame[genotype_col].astype(str).str.strip().eq("")] = np.nan
    return carrier, genotype_col


def safe_n_splits(y: Sequence[int] | pd.Series, requested: int) -> int:
    counts = pd.Series(y).value_counts()
    if len(counts) < 2:
        raise ValueError("Both diagnostic groups are required.")
    n = min(int(requested), int(counts.min()))
    if n < 2:
        raise ValueError("At least two observations per class are required for cross-validation.")
    return n


def read_ids(path: Path) -> list[str]:
    data = pd.read_csv(path)
    column = resolve_column(data, ("SampleId", "sample_id", "ID"))
    ids = data[column].astype(str).str.strip().tolist()
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicated IDs in {path}")
    return ids


def assert_fold_integrity(train_ids: Sequence[str], test_ids: Sequence[str], available_ids: Iterable[str], label: str) -> None:
    train, test, available = set(train_ids), set(test_ids), set(available_ids)
    overlap = train & test
    if overlap:
        raise ValueError(f"{label}: train/test overlap detected ({len(overlap)} IDs).")
    missing = (train | test) - available
    if missing:
        raise KeyError(f"{label}: {len(missing)} fold IDs are absent from the master matrix.")


def load_candidate_genes(path: Path) -> list[str]:
    data = pd.read_csv(path)
    column = resolve_column(data, ("EntrezGeneSymbol", "Protein", "Gene"))
    values = data[column].dropna().astype(str).str.strip()
    return list(dict.fromkeys(v for v in values if v))


def feature_audit(X_train: pd.DataFrame, candidates: Sequence[str], settings: StrictCVSettings) -> tuple[list[str], pd.DataFrame]:
    candidates = list(dict.fromkeys(str(x) for x in candidates if str(x) in X_train.columns))
    if not candidates:
        raise ValueError("No candidate proteins are present in the master matrix.")
    numeric = numeric_frame(X_train[candidates])
    rows = []
    eligible = []
    for feature in candidates:
        series = numeric[feature]
        n_observed = int(series.notna().sum())
        missing_fraction = float(series.isna().mean())
        variance = float(series.var(skipna=True)) if n_observed >= 2 else np.nan
        keep = (
            n_observed >= 2
            and missing_fraction <= settings.max_missing_fraction
            and np.isfinite(variance)
            and variance > settings.min_variance
        )
        if keep:
            eligible.append(feature)
        rows.append({
            "Feature": feature,
            "N_observed_training": n_observed,
            "Missing_fraction_training": missing_fraction,
            "Variance_training": variance,
            "Eligible_training_only": bool(keep),
        })
    audit = pd.DataFrame(rows)
    if not eligible:
        raise ValueError("No candidate protein passed the training-only missingness and variance filters.")
    return eligible, audit


def _elastic_net_search(X: pd.DataFrame, y: pd.Series, settings: StrictCVSettings, seed: int):
    n_splits = safe_n_splits(y, settings.inner_folds)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    estimator = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("classifier", LogisticRegression(
            penalty="elasticnet", solver="saga", max_iter=settings.max_iter,
            tol=settings.tolerance, random_state=seed,
        )),
    ])
    if HAS_SKOPT:
        search = BayesSearchCV(
            estimator=estimator,
            search_spaces={
                "classifier__C": Real(1e-4, 1.0, prior="log-uniform"),
                "classifier__l1_ratio": Categorical([0.85, 0.90, 0.95, 1.0]),
            },
            n_iter=settings.selector_search_iter,
            scoring="roc_auc", cv=cv, n_jobs=settings.n_jobs,
            random_state=seed, refit=True,
        )
    else:
        search = RandomizedSearchCV(
            estimator=estimator,
            param_distributions={
                "classifier__C": loguniform(1e-4, 1.0),
                "classifier__l1_ratio": [0.85, 0.90, 0.95, 1.0],
            },
            n_iter=settings.selector_search_iter,
            scoring="roc_auc", cv=cv, n_jobs=settings.n_jobs,
            random_state=seed, refit=True,
        )
    search.fit(numeric_frame(X), y)
    return search


def select_stable_panel(X_outer_train: pd.DataFrame, y_outer_train: pd.Series, candidates: Sequence[str], settings: StrictCVSettings, seed: int) -> dict[str, Any]:
    outer_eligible, outer_audit = feature_audit(X_outer_train, candidates, settings)
    selector_splits = safe_n_splits(y_outer_train, settings.inner_folds)
    outer_cv = StratifiedKFold(n_splits=selector_splits, shuffle=True, random_state=seed)
    selected_counts = {feature: 0 for feature in outer_eligible}
    coefficient_sums = {feature: 0.0 for feature in outer_eligible}
    coefficient_counts = {feature: 0 for feature in outer_eligible}
    selector_rows = []

    for split, (inner_train_idx, _) in enumerate(outer_cv.split(X_outer_train, y_outer_train), start=1):
        inner_X = X_outer_train.iloc[inner_train_idx]
        inner_y = y_outer_train.iloc[inner_train_idx]
        inner_eligible, _ = feature_audit(inner_X, outer_eligible, settings)
        search = _elastic_net_search(inner_X[inner_eligible], inner_y, settings, seed + split)
        coef = np.asarray(search.best_estimator_.named_steps["classifier"].coef_).reshape(-1)
        for feature, value in zip(inner_eligible, coef):
            coefficient_sums[feature] += abs(float(value))
            coefficient_counts[feature] += 1
            if abs(float(value)) > settings.coefficient_threshold:
                selected_counts[feature] += 1
        selector_rows.append({
            "Inner_split": split,
            "N_training": len(inner_y),
            "N_eligible": len(inner_eligible),
            "Best_C": search.best_params_.get("classifier__C"),
            "Best_l1_ratio": search.best_params_.get("classifier__l1_ratio"),
            "Preprocessing_fitted_inside_CV": True,
        })

    table = outer_audit.set_index("Feature").copy()
    table["Times_selected"] = pd.Series(selected_counts)
    table["Selection_frequency"] = table["Times_selected"] / selector_splits
    table["Mean_abs_coefficient"] = pd.Series({
        feature: coefficient_sums[feature] / max(coefficient_counts[feature], 1)
        for feature in outer_eligible
    })
    table = table.reset_index().sort_values(
        ["Selection_frequency", "Mean_abs_coefficient", "Feature"],
        ascending=[False, False, True],
    )
    panel = table.loc[
        table["Selection_frequency"] >= settings.stability_threshold, "Feature"
    ].tolist()
    fallback_used = False
    if not panel:
        panel = table.loc[table["Eligible_training_only"]].head(settings.fallback_top_n)["Feature"].tolist()
        fallback_used = True
    if not panel:
        raise RuntimeError("Elastic-net stability selection returned an empty panel.")
    return {
        "panel": panel,
        "frequency": table,
        "selector_audit": pd.DataFrame(selector_rows),
        "fallback_used": fallback_used,
        "n_outer_eligible": len(outer_eligible),
    }


def _svm_search(X: pd.DataFrame, y: pd.Series, settings: StrictCVSettings, seed: int):
    n_splits = safe_n_splits(y, settings.inner_folds)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    estimator = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("classifier", SVC(probability=False, random_state=seed)),
    ])
    if HAS_SKOPT:
        search = BayesSearchCV(
            estimator=estimator,
            search_spaces={
                "classifier__kernel": Categorical(["linear", "rbf"]),
                "classifier__C": Real(0.1, 100, prior="log-uniform"),
                "classifier__gamma": Real(1e-4, 1, prior="log-uniform"),
            },
            n_iter=settings.svm_search_iter,
            scoring="roc_auc", cv=cv, n_jobs=settings.n_jobs,
            random_state=seed, refit=True,
        )
    else:
        search = RandomizedSearchCV(
            estimator=estimator,
            param_distributions={
                "classifier__kernel": ["linear", "rbf"],
                "classifier__C": loguniform(0.1, 100),
                "classifier__gamma": loguniform(1e-4, 1),
            },
            n_iter=settings.svm_search_iter,
            scoring="roc_auc", cv=cv, n_jobs=settings.n_jobs,
            random_state=seed, refit=True,
        )
    search.fit(numeric_frame(X), y)
    return search


def fit_tuned_calibrated_svm(X_train: pd.DataFrame, y_train: pd.Series, settings: StrictCVSettings, seed: int):
    search = _svm_search(X_train, y_train, settings, seed)
    params = search.best_params_
    fixed = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("classifier", SVC(
            probability=False,
            kernel=params["classifier__kernel"],
            C=float(params["classifier__C"]),
            gamma=float(params["classifier__gamma"]),
            random_state=seed,
        )),
    ])
    calibration_cv = StratifiedKFold(
        n_splits=safe_n_splits(y_train, settings.inner_folds),
        shuffle=True, random_state=seed + 991,
    )
    calibrated = CalibratedClassifierCV(
        estimator=fixed, method="sigmoid", cv=calibration_cv,
        n_jobs=settings.n_jobs,
    )
    calibrated.fit(numeric_frame(X_train), y_train)
    clean_params = {
        "kernel": params["classifier__kernel"],
        "C": float(params["classifier__C"]),
        "gamma": float(params["classifier__gamma"]),
        "Search": "BayesSearchCV" if HAS_SKOPT else "RandomizedSearchCV",
        "Calibration": "sigmoid_inner_CV",
    }
    return calibrated, clean_params


def predict_probability(estimator, X: pd.DataFrame) -> np.ndarray:
    probabilities = estimator.predict_proba(numeric_frame(X))
    classes = list(estimator.classes_)
    return probabilities[:, classes.index(1)]


def metric_row(y_true: Sequence[int], probability: Sequence[float], **extra: Any) -> dict[str, Any]:
    y_true = np.asarray(y_true, dtype=int)
    probability = np.asarray(probability, dtype=float)
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    row = {
        "AUC": roc_auc_score(y_true, probability),
        "Accuracy": accuracy_score(y_true, prediction),
        "Sensitivity": recall_score(y_true, prediction, zero_division=0),
        "Specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "Precision": precision_score(y_true, prediction, zero_division=0),
        "F1": f1_score(y_true, prediction, zero_division=0),
        "N_test": len(y_true),
        "N_CN_test": int((y_true == 0).sum()),
        "N_AD_test": int((y_true == 1).sum()),
    }
    row.update(extra)
    return row


def roc_table(y_true: Sequence[int], probability: Sequence[float], **extra: Any) -> pd.DataFrame:
    fpr, tpr, thresholds = roc_curve(y_true, probability)
    table = pd.DataFrame({"FPR": fpr, "TPR": tpr, "Threshold": thresholds})
    table["AUC"] = roc_auc_score(y_true, probability)
    for key, value in extra.items():
        table[key] = value
    return table


def fold_permutation_importance(estimator, X_test: pd.DataFrame, y_test: pd.Series, features: Sequence[str], settings: StrictCVSettings, fold_label: str) -> pd.DataFrame:
    if pd.Series(y_test).nunique() < 2:
        return pd.DataFrame()
    result = permutation_importance(
        estimator, numeric_frame(X_test[list(features)]), y_test,
        scoring="roc_auc", n_repeats=settings.permutation_repeats,
        random_state=settings.random_state, n_jobs=settings.n_jobs,
    )
    return pd.DataFrame({
        "Proteins": list(features),
        "Permutation": result.importances_mean,
        "Permutation_sd_within_fold": result.importances_std,
        "Fold": fold_label,
        "Importance_scope": "outer_test_fold",
    })


def aggregate_importance(fold_tables: Sequence[pd.DataFrame], selection_frequency: pd.Series | None = None) -> pd.DataFrame:
    valid = [x for x in fold_tables if x is not None and not x.empty]
    if not valid:
        return pd.DataFrame(columns=["Proteins", "Permutation", "Permutation_sd", "N_folds_evaluated"])
    combined = pd.concat(valid, ignore_index=True)
    out = combined.groupby("Proteins", as_index=False).agg(
        Permutation=("Permutation", "mean"),
        Permutation_sd=("Permutation", "std"),
        N_folds_evaluated=("Fold", "nunique"),
    )
    out["Permutation_sd"] = out["Permutation_sd"].fillna(0.0)
    if selection_frequency is not None:
        out["Selection_frequency"] = out["Proteins"].map(selection_frequency).fillna(0.0)
    out["Importance_scope"] = "mean_outer_test_fold_importance"
    return out.sort_values("Permutation", ascending=False)


def aggregate_panel_frequency(panels: Sequence[Sequence[str]]) -> pd.Series:
    counts: dict[str, int] = {}
    for panel in panels:
        for feature in set(panel):
            counts[feature] = counts.get(feature, 0) + 1
    denominator = max(len(panels), 1)
    return pd.Series({feature: count / denominator for feature, count in counts.items()}, dtype=float).sort_values(ascending=False)


def compute_midrank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    sorted_values = values[order]
    ranks = np.zeros(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i
        while j < len(values) and sorted_values[j] == sorted_values[i]:
            j += 1
        ranks[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    output = np.empty(len(values), dtype=float)
    output[order] = ranks
    return output


def fast_delong(predictions: np.ndarray, positive_count: int):
    m = positive_count
    n = predictions.shape[1] - m
    k = predictions.shape[0]
    tx = np.empty((k, m + n))
    for row in range(k):
        tx[row] = compute_midrank(predictions[row])
    aucs = (tx[:, :m].mean(axis=1) - tx[:, m:].mean(axis=1)) / n
    v01 = (tx[:, :m] - tx[:, :m].mean(axis=1)[:, None]) / n
    v10 = (tx[:, m:] - tx[:, m:].mean(axis=1)[:, None]) / n
    covariance = np.cov(v01) / m + np.cov(v10) / n
    return aucs, covariance


def delong_pairwise(y_true: Sequence[int], score_a: Sequence[float], score_b: Sequence[float]) -> dict[str, float]:
    y = np.asarray(y_true, dtype=int)
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    valid = np.isfinite(y) & np.isfinite(a) & np.isfinite(b)
    y, a, b = y[valid], a[valid], b[valid]
    positive = np.where(y == 1)[0]
    negative = np.where(y == 0)[0]
    predictions = np.vstack([
        np.concatenate([a[positive], a[negative]]),
        np.concatenate([b[positive], b[negative]]),
    ])
    _, covariance = fast_delong(predictions, len(positive))
    auc_a = roc_auc_score(y, a)
    auc_b = roc_auc_score(y, b)
    variance = covariance[0, 0] + covariance[1, 1] - 2 * covariance[0, 1]
    if not np.isfinite(variance) or variance <= 1e-12:
        return {"AUC_1": auc_a, "AUC_2": auc_b, "Delta_AUC": auc_b - auc_a, "Z": np.nan, "P": np.nan}
    z_value = (auc_a - auc_b) / np.sqrt(variance)
    p_value = 2 * (1 - norm.cdf(abs(z_value)))
    return {"AUC_1": auc_a, "AUC_2": auc_b, "Delta_AUC": auc_b - auc_a, "Z": z_value, "P": p_value}


def settings_table(settings: StrictCVSettings) -> pd.DataFrame:
    return pd.DataFrame([{
        **settings.__dict__,
        "preprocessing_scope": "inside_each_training_CV_split",
        "oof_score_type": "calibrated_probability_without_foldwise_rescaling",
        "importance_scope": "outer_test_folds_only",
        "selector_search": "BayesSearchCV" if HAS_SKOPT else "RandomizedSearchCV",
    }])


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
