# Run the WGCNA reporting workflow
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
  "10_WGCNA_generate_main_figure3.R",
  "11_WGCNA_generate_extended_data_figures.R",
  "12_WGCNA_generate_supplementary_tables.R",
  "13_WGCNA_generate_manuscript_text.R",
  "14_WGCNA_build_submission_package.R",
  "15_WGCNA_audit_submission_package.R"
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
