# Run the WGCNA analysis workflow
rm(list = ls())
options(stringsAsFactors = FALSE)

project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) {
  normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Run renv::restore().", call. = FALSE)
}

scripts <- c(
  "01_WGCNA_prepare_input.R",
  "02_WGCNA_construct_network.R",
  "03_WGCNA_characterize_modules.R",
  "04_WGCNA_module_trait_models.R",
  "05_WGCNA_site_biomarker_sensitivity.R",
  "06_WGCNA_association_stability.R",
  "07_WGCNA_biomarker_fdr_correction.R",
  "08_WGCNA_structural_preservation.R",
  "09_WGCNA_network_quality.R"
)
source(file.path(project_root, "R", "wgcna_bootstrap.R"), local = FALSE)
WGCNA_CONFIG <- wgcna_load_config(project_root)
log_dir <- file.path(WGCNA_CONFIG$result_root, "workflow_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(REDLAT_PROJECT_ROOT = project_root)
rscript <- file.path(R.home("bin"), "Rscript")

args <- commandArgs(trailingOnly = TRUE)
force_rerun <- "--force" %in% args

for (script in scripts) {
  path <- file.path(project_root, "scripts", "WGCNA", script)
  if (!file.exists(path)) stop("Missing workflow script: ", path, call. = FALSE)
  completed_flag <- file.path(log_dir, paste0(script, ".completed.txt"))
  if (!force_rerun && file.exists(completed_flag)) {
    message("\n=== Skipping completed ", script, " ===")
    next
  }
  message("\n=== Running ", script, " ===")
  started <- Sys.time()
  status <- system2(
    rscript,
    c("--vanilla", shQuote(path))
  )
  if (!identical(status, 0L)) {
    writeLines(c(
      paste0("SCRIPT=", script),
      paste0("STARTED=", started),
      paste0("FAILED=", Sys.time()),
      paste0("STATUS=", status)
    ), file.path(log_dir, paste0(script, ".failed.txt")))
    stop("WGCNA workflow stopped in ", script, call. = FALSE)
  }
  writeLines(c(
    paste0("SCRIPT=", script),
    paste0("STARTED=", started),
    paste0("COMPLETED=", Sys.time())
  ), file.path(log_dir, paste0(script, ".completed.txt")))
}
