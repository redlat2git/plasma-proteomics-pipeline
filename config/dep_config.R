dep_project_config <- function(project_root) {
  data_dir <- Sys.getenv(
    "REDLAT_REVIEWER_DATA_DIR",
    unset = file.path(project_root, "data_private", "reviewer_inputs")
  )
  list(
    project_root = project_root,
    data_dir = data_dir,
    metadata_file = file.path(data_dir, "ReDLat_metadata_deidentified.csv"),
    proteomics_file = file.path(data_dir, "ReDLat_proteomics_somamer_log2.csv"),
    annotation_file = file.path(data_dir, "ReDLat_feature_annotation.csv"),
    # Legacy compatibility alias only; this is NOT an ADAT.
    adat_file = file.path(data_dir, "ReDLat_proteomics_somamer_log2.csv"),
    result_root = file.path(project_root, "outputs", "DEP"),
    publication_root = file.path(project_root, "outputs", "publication", "DEP"),
    allow_participant_level_exports = tolower(trimws(Sys.getenv("REDLAT_ALLOW_PARTICIPANT_EXPORTS", unset = "FALSE"))) %in% c("1", "true", "t", "yes", "y")
  )
}
