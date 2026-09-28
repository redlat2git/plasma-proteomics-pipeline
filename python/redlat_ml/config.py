from dataclasses import dataclass
from pathlib import Path
import os
import re

@dataclass(frozen=True)
class Config:
    project_root: Path
    data_dir: Path
    metadata_file: Path
    master_file: Path
    proteomics_file: Path
    annotation_file: Path
    adat_file: Path
    private_root: Path
    publication_root: Path
    public_root: Path


def _find_root(start: Path) -> Path:
    env = os.environ.get("REDLAT_PROJECT_ROOT", "").strip()
    if env:
        root = Path(env).expanduser().resolve()
        if not root.exists():
            raise FileNotFoundError(f"REDLAT_PROJECT_ROOT does not exist: {root}")
        return root
    for candidate in [start, *start.parents]:
        if (candidate / ".redlat-root").exists() or (candidate / ".here").exists():
            return candidate
    raise FileNotFoundError("Could not locate project root.")


def _reviewer_data_dir(root: Path) -> Path:
    env = os.environ.get("REDLAT_REVIEWER_DATA_DIR", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return root / "data_private" / "reviewer_inputs"


def load_config(script_file=None) -> Config:
    start = Path(script_file).resolve() if script_file else Path.cwd().resolve()
    root = _find_root(start)
    data = _reviewer_data_dir(root)
    private = root / "outputs" / "ML" / "private"
    publication = root / "outputs" / "ML" / "publication"
    return Config(
        project_root=root,
        data_dir=data,
        metadata_file=private / "derived" / "ReDLat_metadata_ML_compat.csv",
        master_file=private / "derived" / "ReDLat_ML_gene_master_RAW.csv",
        proteomics_file=data / "ReDLat_proteomics_somamer_log2.csv",
        annotation_file=data / "ReDLat_feature_annotation.csv",
        # Legacy alias; never parsed as an ADAT in the processed-data workflow.
        adat_file=data / "ReDLat_proteomics_somamer_log2.csv",
        private_root=private,
        publication_root=publication,
        public_root=publication,
    )


def require_files(*paths):
    resolved = []
    for item in paths:
        if isinstance(item, (list, tuple)) and not isinstance(item, Path):
            if len(item) == 2 and isinstance(item[1], str) and isinstance(item[0], (str, Path)):
                resolved.append(Path(item[0]))
            else:
                for sub in item:
                    if isinstance(sub, (list, tuple)) and len(sub) >= 1:
                        resolved.append(Path(sub[0]))
                    else:
                        resolved.append(Path(sub))
        else:
            resolved.append(Path(item))
    missing = [str(p) for p in resolved if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))
    return True


def assert_public_table(data, name="public_table"):
    """Reject participant identifiers, pseudonyms, and personal paths in public tables."""
    if data is None:
        return True

    forbidden_columns = {
        "sampleid", "sample_id", "participantid", "participant_id",
        "subjectid", "subject_id", "record_id", "match_pair_id",
        "subclass", "study_id", "reviewer_id",
    }
    observed = {str(c).strip().lower() for c in data.columns}
    bad_columns = sorted(observed & forbidden_columns)
    if bad_columns:
        raise ValueError(
            f"{name}: participant-level identifier columns are not allowed "
            f"in public Source Data: {bad_columns}"
        )

    path_patterns = [
        re.compile(r"[A-Za-z]:[/\\]Users[/\\]", re.I),
        re.compile(r"/Users/[^/]+/", re.I),
        re.compile(r"/home/[^/]+/", re.I),
    ]
    pseudonym_pattern = re.compile(r"\bRDLAT_[A-F0-9]{8,}\b", re.I)

    object_columns = data.select_dtypes(include=["object", "string"]).columns
    for column in object_columns:
        for value in data[column].dropna().astype(str):
            if any(pattern.search(value) for pattern in path_patterns):
                raise ValueError(
                    f"{name}: local personal filesystem path detected in column '{column}'."
                )
            if pseudonym_pattern.search(value):
                raise ValueError(
                    f"{name}: participant pseudonym detected in public Source Data "
                    f"column '{column}'."
                )
    return True
