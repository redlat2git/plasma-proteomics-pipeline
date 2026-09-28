from pathlib import Path
import os, subprocess, sys
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("REDLAT_PROJECT_ROOT", str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")
steps = [
    ROOT / "scripts/ML/25_ML_audit_strict_pipeline.py",
    ROOT / "scripts/ML/18_ML_plot_primary_results.py",
    ROOT / "scripts/ML/19_ML_plot_ptau_results.py",
    ROOT / "scripts/ML/20_ML_plot_delong_results.py",
    ROOT / "scripts/ML/21_ML_plot_loco_results.py",
    ROOT / "scripts/ML/22_ML_plot_loco_forest.py",
    ROOT / "scripts/ML/26_ML_generate_extended_data_figure10.py",
    ROOT / "scripts/ML/23_ML_generate_source_data.py",
    ROOT / "scripts/ML/24_ML_audit_publication_outputs.py",
]
for script in steps:
    print(f"\n>>> {script.name}", flush=True)
    subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)
