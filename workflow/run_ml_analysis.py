from pathlib import Path
import os, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("REDLAT_PROJECT_ROOT", str(ROOT))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

def find_rscript() -> str:
    env_r = os.environ.get("RSCRIPT_PATH")
    if env_r and Path(env_r).exists():
        return env_r
    which_r = shutil.which("Rscript")
    if which_r:
        return which_r
    prog_files = Path(r"C:\Program Files\R")
    if prog_files.exists():
        candidates = sorted(prog_files.glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates and candidates[0].exists():
            return str(candidates[0])
    return "Rscript"

RSCRIPT = find_rscript()

steps = [
    (sys.executable, ROOT / "scripts/ML/25_ML_audit_strict_pipeline.py"),
    (RSCRIPT, ROOT / "scripts/ML/01_ML_create_primary_folds.R"),
    (RSCRIPT, ROOT / "scripts/ML/02_ML_fit_primary_dep_folds.R"),
    (sys.executable, ROOT / "scripts/ML/03_ML_run_primary_nested_cv.py"),
    (RSCRIPT, ROOT / "scripts/ML/04_ML_create_ptau_folds.R"),
    (RSCRIPT, ROOT / "scripts/ML/05_ML_fit_ptau_dep_folds.R"),
    (sys.executable, ROOT / "scripts/ML/06_ML_run_ptau_nested_cv.py"),
    (RSCRIPT, ROOT / "scripts/ML/07_ML_create_apoe_folds.R"),
    (RSCRIPT, ROOT / "scripts/ML/08_ML_fit_apoe_dep_folds.R"),
    (sys.executable, ROOT / "scripts/ML/09_ML_run_apoe_nested_cv.py"),
    (RSCRIPT, ROOT / "scripts/ML/10_ML_create_loco_folds.R"),
    (RSCRIPT, ROOT / "scripts/ML/11_ML_fit_loco_dep_folds.R"),
    (sys.executable, ROOT / "scripts/ML/12_ML_run_loco.py"),
    (sys.executable, ROOT / "scripts/ML/13_ML_evaluate_fixed_panel_ptau.py"),
    (sys.executable, ROOT / "scripts/ML/14_ML_run_clinical_regression.py"),
    (RSCRIPT, ROOT / "scripts/ML/15_ML_rebuild_matching.R"),
    (sys.executable, ROOT / "scripts/ML/16_ML_run_matched_analysis.py"),
    (sys.executable, ROOT / "scripts/ML/17_ML_meta_analyze_loco.py"),
]
for executable, script in steps:
    print(f"\n>>> {script.name}", flush=True)
    subprocess.run([str(executable), str(script)], check=True, cwd=ROOT)
