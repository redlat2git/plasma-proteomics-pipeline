# Run the complete WGCNA workflow
rm(list = ls())
project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE) else here::here()
sys.source(file.path(project_root, "workflow", "run_wgcna_analysis.R"), envir = new.env(parent = globalenv()))
sys.source(file.path(project_root, "workflow", "run_wgcna_reporting.R"), envir = new.env(parent = globalenv()))
