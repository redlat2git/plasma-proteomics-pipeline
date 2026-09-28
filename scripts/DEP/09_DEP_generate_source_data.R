###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 09. DEP Source Data
# Requires: Completed DEP analysis and reporting outputs
# Produces: Figure-specific Source Data workbooks
# Data policy: participant-level data and intermediate outputs remain local.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# Project setup
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the project environment with renv::restore().", call. = FALSE)
}
source(file.path(project_root, "R", "dep_bootstrap.R"), local = FALSE)
DEP_CONFIG <- dep_load_config(project_root)
analysis_root <- DEP_CONFIG$result_root
publication_root <- DEP_CONFIG$publication_root
options(stringsAsFactors = FALSE)

required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "stringr", "tibble", "openxlsx")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


EXPECTED_MAP_N <- 9638L
MAIN_FDR <- 0.05

if (!isTRUE(DEP_CONFIG$allow_participant_level_exports)) {
  stop(
    "Source Data generation includes de-identified participant-level plotted values. ",
    "Set allow_participant_level_exports = TRUE in the untracked local config only after governance approval.",
    call. = FALSE
  )
}
OUT_ROOT <- file.path(publication_root, "source_data")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Helpers
###############################################################################

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

read_optional_csv <- function(paths) {
  path <- first_existing_file(paths)
  if (is.na(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

read_required_csv <- function(paths, label) {
  out <- read_optional_csv(paths)
  if (is.null(out)) stop("Required source not found: ", label, call. = FALSE)
  out
}

clean_text_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A")] <- NA_character_
  x
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_log2 <- function(x) {
  x <- safe_numeric(x)
  x[x <= 0] <- NA_real_
  log2(x)
}

pick_col <- function(df, candidates) {
  if (is.null(df) || length(names(df)) == 0) return(NA_character_)
  exact <- candidates[candidates %in% names(df)][1]
  if (length(exact) > 0 && !is.na(exact)) return(exact)
  clean <- function(z) tolower(gsub("[^a-z0-9]+", "", z))
  idx <- match(clean(candidates), clean(names(df)), nomatch = 0)
  idx <- idx[idx > 0][1]
  if (length(idx) == 0 || is.na(idx)) return(NA_character_)
  names(df)[idx]
}

assert_fixed_map <- function(tbl, label, expected_n = EXPECTED_MAP_N) {
  if (is.null(tbl)) stop(label, " is NULL.", call. = FALSE)
  if (nrow(tbl) != expected_n) {
    stop(label, " must contain ", expected_n, " rows; observed ", nrow(tbl), ".", call. = FALSE)
  }
  if ("EntrezGeneSymbol" %in% names(tbl) && dplyr::n_distinct(tbl$EntrezGeneSymbol) != expected_n) {
    stop(label, " does not contain ", expected_n, " unique gene symbols.", call. = FALSE)
  }
  if ("AptName" %in% names(tbl) && dplyr::n_distinct(tbl$AptName) != expected_n) {
    stop(label, " does not contain ", expected_n, " unique SOMAmer identifiers.", call. = FALSE)
  }
  invisible(TRUE)
}

sanitize_sheet <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(x, 1, 31)
}

make_styles <- function() {
  list(
    title = openxlsx::createStyle(
      fontName = "Arial", fontSize = 12, textDecoration = "bold",
      fontColour = "#111111", fgFill = "#CFC7B7",
      halign = "center", valign = "center",
      border = "TopBottomLeftRight", borderColour = "#4A4A4A"
    ),
    note = openxlsx::createStyle(
      fontName = "Arial", fontSize = 9, textDecoration = "italic",
      fontColour = "#555555", fgFill = "#F4F2EC",
      wrapText = TRUE, valign = "center",
      border = "TopBottomLeftRight", borderColour = "#4A4A4A"
    ),
    header = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10, textDecoration = "bold",
      fontColour = "#111111", fgFill = "#E9E5DC",
      halign = "center", valign = "center", wrapText = TRUE,
      border = "TopBottomLeftRight", borderColour = "#4A4A4A"
    ),
    body = openxlsx::createStyle(
      fontName = "Arial", fontSize = 9, valign = "top", wrapText = TRUE,
      border = "TopBottomLeftRight", borderColour = "#4A4A4A"
    )
  )
}

relative_source_paths <- function(paths) {
  root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  vapply(paths, function(path) {
    if (is.na(path) || !nzchar(path)) return(NA_character_)
    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    prefix <- paste0(root, "/")
    if (startsWith(normalized, prefix)) {
      substring(normalized, nchar(prefix) + 1L)
    } else {
      basename(normalized)
    }
  }, character(1))
}

source_table_counter <- 0L

write_source_sheet <- function(wb, sheet, title, note, data, styles) {
  sheet <- sanitize_sheet(sheet)
  data <- if (is.null(data)) NULL else as.data.frame(data, check.names = FALSE)
  if (!is.null(data)) names(data) <- make.unique(names(data), sep = "_dup")
  ncols <- if (is.null(data) || ncol(data) == 0) 2L else max(ncol(data), 2L)

  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)

  openxlsx::mergeCells(wb, sheet, cols = seq_len(ncols), rows = 1)
  openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet, styles$title, rows = 1, cols = seq_len(ncols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 25)

  openxlsx::mergeCells(wb, sheet, cols = seq_len(ncols), rows = 2)
  openxlsx::writeData(wb, sheet, note, startRow = 2, startCol = 1)
  openxlsx::addStyle(wb, sheet, styles$note, rows = 2, cols = seq_len(ncols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 34)

  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    openxlsx::writeData(wb, sheet, "No data available.", startRow = 4, startCol = 1)
    return(invisible(sheet))
  }

  source_table_counter <<- source_table_counter + 1L
  table_name <- paste0("SrcTbl_", source_table_counter)

  openxlsx::writeDataTable(
    wb, sheet, data,
    startRow = 4, startCol = 1,
    tableName = table_name,
    tableStyle = "TableStyleLight1",
    withFilter = TRUE
  )
  openxlsx::addStyle(
    wb, sheet, styles$header,
    rows = 4, cols = seq_len(ncol(data)),
    gridExpand = TRUE, stack = TRUE
  )
  openxlsx::addStyle(
    wb, sheet, styles$body,
    rows = 5:(4 + nrow(data)), cols = seq_len(ncol(data)),
    gridExpand = TRUE, stack = TRUE
  )
  openxlsx::freezePane(wb, sheet, firstActiveRow = 5)

  widths <- vapply(names(data), function(nm) {
    if (grepl("description|targetfullname|core_enrichment|leading_edge|notes?|content|formula|source", nm, ignore.case = TRUE)) {
      34
    } else {
      min(max(nchar(nm) + 2, 11), 22)
    }
  }, numeric(1))
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = widths)
  invisible(sheet)
}

write_readme <- function(wb, figure_name, status, notes, source_files, styles) {
  readme <- tibble::tibble(
    Field = c("Figure", "Status", "Purpose", "Notes", "Privacy", "Source files"),
    Content = c(
      figure_name,
      status,
      "Exact numerical values required to reconstruct the quantitative panels in this workbook.",
      notes,
      "No direct or stable study identifiers are included. Participant-level values are limited to variables required to reconstruct the displayed panels.",
      paste(relative_source_paths(source_files), collapse = "\n")
    )
  )
  write_source_sheet(
    wb, "README",
    paste0("Source Data — ", figure_name),
    "Prepared for Nature Aging submission.",
    readme, styles
  )
}

save_source_workbook <- function(filename, figure_name, sheets, status, notes, source_files) {
  for (sheet_spec in sheets) {
    dep_assert_no_direct_identifiers(sheet_spec$data, paste0(figure_name, " / ", sheet_spec$sheet))
  }
  wb <- openxlsx::createWorkbook()
  styles <- make_styles()
  write_readme(wb, figure_name, status, notes, source_files, styles)
  for (x in sheets) {
    write_source_sheet(wb, x$sheet, x$title, x$note, x$data, styles)
  }
  path <- file.path(OUT_ROOT, filename)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

spearman_long <- function(dat, protein_map, trait_map, fdr_scope = "all_displayed_pairs") {
  out <- purrr::map_dfr(seq_len(nrow(protein_map)), function(i) {
    protein <- protein_map[i, , drop = FALSE]
    apt <- protein$AptName[[1]]
    if (!apt %in% names(dat)) return(tibble::tibble())
    x <- safe_log2(dat[[apt]])

    purrr::map_dfr(seq_len(nrow(trait_map)), function(j) {
      trait <- trait_map[j, , drop = FALSE]
      col <- trait$column[[1]]
      if (is.na(col) || !col %in% names(dat)) return(tibble::tibble())
      y <- dat[[col]]
      if (is.factor(y) || is.character(y)) {
        y_chr <- trimws(as.character(y))
        if (trait$trait[[1]] == "Diagnosis") {
          y <- ifelse(y_chr == "AD", 1, ifelse(y_chr == "CN", 0, NA_real_))
        } else if (trait$trait[[1]] == "Sex") {
          y <- ifelse(tolower(y_chr) %in% c("f", "female", "1"), 1,
                      ifelse(tolower(y_chr) %in% c("m", "male", "2"), 0, NA_real_))
        } else if (trait$trait[[1]] == "APOE ε4") {
          y <- suppressWarnings(as.numeric(y_chr))
        } else {
          y <- suppressWarnings(as.numeric(y_chr))
        }
      } else {
        y <- safe_numeric(y)
      }
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      if (n < 6) return(tibble::tibble())
      test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      tibble::tibble(
        Protein = protein$EntrezGeneSymbol[[1]],
        Protein_Name = protein$Protein_Name[[1]],
        AptName = apt,
        Trait = trait$trait[[1]],
        Trait_column = col,
        n = n,
        rho = unname(test$estimate),
        P.Value = test$p.value,
        FDR_scope = fdr_scope
      )
    })
  })
  if (nrow(out) > 0) out$adj.P.Val <- p.adjust(out$P.Value, method = "BH")
  out
}

###############################################################################
# Load canonical workspace and principal outputs
###############################################################################

workspace_file <- first_existing_file(c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
))
if (is.na(workspace_file)) {
  stop("analysis_workspace.RData was not found. Run Script 01 first.", call. = FALSE)
}
load(workspace_file)

if (!exists("DEP_gene") || !exists("dep_df")) {
  stop("Workspace must contain DEP_gene and dep_df.", call. = FALSE)
}
assert_fixed_map(DEP_gene, "Primary DEP_gene")

primary_full_path <- first_existing_file(c(
  file.path(analysis_root, "03_dep", "gene_collapsed", "AD_vs_CN_full_limma_results_gene_collapsed.csv")
))
primary_full <- read_required_csv(primary_full_path, "primary gene-collapsed DEP")
assert_fixed_map(primary_full, "Primary gene-collapsed DEP")

output_registry <- list()
coverage <- list()

record_output <- function(figure, file, status, missing = NA_character_) {
  output_registry[[length(output_registry) + 1L]] <<- tibble::tibble(
    Figure = figure,
    Output_file = ifelse(is.na(file), NA_character_, basename(file)),
    Status = status,
    Missing_requirement = missing
  )
}

###############################################################################
# Figure 2
###############################################################################

reactome_path <- first_existing_file(c(file.path(analysis_root, "05_enrichment_corrected", "gsea", "main_dep_gsea_reactome_bh.csv")))
reactome <- read_required_csv(reactome_path, "primary Reactome GSEA")

pathway_display <- tibble::tribble(
  ~Display_label, ~search_pattern,
  "Assembly of collagen fibrils and other multimeric structures", "assembly of collagen fibrils",
  "Degradation of the extracellular matrix", "degradation of the extracellular matrix",
  "Extracellular matrix organization", "extracellular matrix organization",
  "Metabolism of RNA", "metabolism of rna",
  "Processing of capped intron-containing pre-mRNA", "processing of capped intron-containing pre-mrna",
  "mRNA 3'-end processing", "mrna 3'-end processing"
)
reactome_display <- purrr::map_dfr(seq_len(nrow(pathway_display)), function(i) {
  hit <- reactome %>%
    dplyr::filter(stringr::str_detect(tolower(Description), stringr::fixed(pathway_display$search_pattern[i]))) %>%
    dplyr::slice(1)
  if (nrow(hit) == 0) return(tibble::tibble(Display_label = pathway_display$Display_label[i], note = "Term not found"))
  hit %>% dplyr::mutate(Display_label = pathway_display$Display_label[i], .before = 1)
})

fig2_proteins <- c("SPC25", "LRRN1", "ODC1", "CPLX2", "SMOC1", "ACHE", "PLCH1", "RNPEP", "C3", "PUM2", "CAPN1", "GDI1")
fig2_map <- DEP_gene %>%
  dplyr::filter(EntrezGeneSymbol %in% fig2_proteins) %>%
  dplyr::mutate(order = match(EntrezGeneSymbol, fig2_proteins)) %>%
  dplyr::arrange(order) %>%
  dplyr::select(Protein_Name, EntrezGeneSymbol, AptName, dplyr::everything())

trait_candidates_fig2 <- list(
  `CDR-SB` = c("cdr_boxscore", "CDR_SB", "cdr_sb"),
  PFAQ = c("udsfaq_total", "PFAQ"),
  MMSE = c("mmse_total", "MMSE"),
  `Mini-SEA` = c("Mini.SEA", "Mini-SEA"),
  `NPI-Q` = c("NPI", "NPI_Q"),
  `p-tau181` = c("p.tau181", "p-tau181", "p_tau181"),
  `p-tau217` = c("p.tau217", "p-tau217", "p_tau217"),
  NfL = c("NfL", "NFL"),
  GFAP = c("GFAP"),
  `Aβ42/40` = c("ratio.AB42.40", "ratio AB42/40", "ratio_AB42_40")
)
trait_map_fig2 <- tibble::tibble(
  trait = names(trait_candidates_fig2),
  column = vapply(trait_candidates_fig2, function(x) pick_col(dep_df, x), character(1))
)
fig2_corr <- spearman_long(dep_df, fig2_map, trait_map_fig2)

apoe_comparison <- read_required_csv(c(file.path(analysis_root, "04_sensitivity", "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv")), "APOE comparison")
cdr_alignment <- read_required_csv(c(file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv")), "AD-only CDR-SB alignment")
fig2_annotations <- fig2_map %>%
  dplyr::select(Protein_Name, EntrezGeneSymbol, AptName) %>%
  dplyr::left_join(primary_full %>% dplyr::select(EntrezGeneSymbol, AptName, primary_logFC = logFC, primary_FDR = adj.P.Val, primary_Direction = Direction), by = c("EntrezGeneSymbol", "AptName")) %>%
  dplyr::left_join(apoe_comparison %>% dplyr::select(EntrezGeneSymbol, AptName, APOE_logFC = logFC_apoe, APOE_FDR = adj.P.Val_apoe, APOE_same_direction = same_direction, APOE_preserved_FDR005 = preserved_fdr005), by = c("EntrezGeneSymbol", "AptName")) %>%
  dplyr::left_join(cdr_alignment %>% dplyr::select(EntrezGeneSymbol, AptName, CDRSB_beta = beta_CDRSB_AD_only, CDRSB_FDR = adj.P.Val_CDRSB_AD_only, CDRSB_same_direction = same_direction), by = c("EntrezGeneSymbol", "AptName"))

atn_comparison_path <- first_existing_file(c(file.path(analysis_root, "04_sensitivity", "atn_adjusted", "primary_vs_ATN_adjusted_gene_comparison.csv")))
atn_comparison <- read_required_csv(atn_comparison_path, "primary versus AT(N)-adjusted comparison")
loco_summary_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "country_loco", "tables", "LOCO_summary_metrics_FIXED_PRIMARY_MAP.csv")))
loco_summary <- read_required_csv(loco_summary_path, "fixed-map LOCO summary")
loco_mean_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "country_loco", "tables", "main_vs_meanLOCO_table_FIXED_PRIMARY_MAP.csv")))
loco_mean <- read_required_csv(loco_mean_path, "fixed-map mean LOCO table")
assert_fixed_map(loco_mean, "Mean LOCO protein table")

fig2_file <- save_source_workbook(
  "SourceData_Fig2_DEP.xlsx",
  "Figure 2 — DEP architecture and multicountry stability",
  list(
    list(sheet = "Fig2a_volcano", title = "Fig. 2a primary differential-abundance source data", note = "One row per fixed gene–SOMAmer pair.", data = primary_full),
    list(sheet = "Fig2a_displayed_labels", title = "Proteins labeled in Fig. 2a", note = "Labels visible in the submitted figure.", data = primary_full %>% dplyr::filter(EntrezGeneSymbol %in% c("SPC25", "LRRN1", "GJB3", "CPLX2", "CAPN1", "PUM2", "AMBRA1", "C3"))),
    list(sheet = "Fig2b_Reactome_full", title = "Fig. 2b full Reactome GSEA", note = "Complete BH-adjusted pre-ranked Reactome results.", data = reactome),
    list(sheet = "Fig2b_displayed_paths", title = "Pathways displayed in Fig. 2b", note = "Six representative pathways shown in the panel.", data = reactome_display),
    list(sheet = "Fig2c_correlations", title = "Fig. 2c protein–trait correlations", note = "BH adjustment applied across all protein–trait pairs displayed in Fig. 2c.", data = fig2_corr),
    list(sheet = "Fig2c_annotations", title = "Fig. 2c protein side annotations", note = "Primary direction, APOE preservation and within-AD severity alignment.", data = fig2_annotations),
    list(sheet = "Fig2d_ATN_comparison", title = "Fig. 2d primary versus AT(N)-adjusted effects", note = "Fixed primary map; one row per gene–SOMAmer pair.", data = atn_comparison),
    list(sheet = "Fig2e_LOCO_summary", title = "Fig. 2e country-exclusion summary", note = "Effect-size concordance and primary-DEP preservation by excluded country.", data = loco_summary),
    list(sheet = "Fig2f_mean_LOCO", title = "Fig. 2f primary versus mean LOCO effects", note = "One row per fixed gene–SOMAmer pair.", data = loco_mean)
  ),
  status = "Complete",
  notes = "The workbook includes exact data for all six panels. Protein–trait correlations are recalculated from non-residualized log2 RFU values in the canonical workspace.",
  source_files = c(primary_full_path, reactome_path,
    file.path(analysis_root, "04_sensitivity", "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"),
    file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"),
    atn_comparison_path, loco_summary_path, loco_mean_path)
)
record_output("Fig. 2", fig2_file, "GENERATED")

###############################################################################
# Supplementary / Extended Data Fig. 1: PCA
###############################################################################

pca_scores_path <- first_existing_file(c(file.path(analysis_root, "02_pca", "pca_scores_privacy_safe.csv")))
pca_var_path <- first_existing_file(c(file.path(analysis_root, "02_pca", "pca_variance_explained.csv")))

pca_scores_tbl <- if (!is.na(pca_scores_path)) readr::read_csv(pca_scores_path, show_col_types = FALSE, guess_max = 100000) else if (exists("pca_df")) tibble::as_tibble(get("pca_df")) else NULL
pca_var_tbl <- if (!is.na(pca_var_path)) readr::read_csv(pca_var_path, show_col_types = FALSE) else if (exists("pca_var")) tibble::as_tibble(get("pca_var")) else NULL

if (!is.null(pca_scores_tbl) && !is.null(pca_var_tbl)) {
  if (!"PlateId" %in% names(pca_scores_tbl)) {
    pca_scores_tbl$PlateId <- if ("Plate" %in% names(pca_scores_tbl)) pca_scores_tbl$Plate else NA_character_
  }
  if (!"site" %in% names(pca_scores_tbl)) {
    pca_scores_tbl$site <- if ("Site" %in% names(pca_scores_tbl)) pca_scores_tbl$Site else if ("SiteId" %in% names(pca_scores_tbl)) pca_scores_tbl$SiteId else NA_character_
  }
  if (!"SiteId" %in% names(pca_scores_tbl)) {
    pca_scores_tbl$SiteId <- pca_scores_tbl$site
  }
  pca_keep <- c("PC1", "PC2", "PC3", "SampleGroup", "Sex", "Age", "APOE4_carrier", "PlateId", "Country", "site", "SiteId")
  pca_keep <- pca_keep[pca_keep %in% names(pca_scores_tbl)]
  pca_file <- save_source_workbook(
    "SourceData_DEP_PCA.xlsx",
    "Supplementary / Extended Data figure — PCA quality control",
    list(
      list(sheet = "PCA_scores", title = "PCA sample scores and coloring variables", note = "Participant-level PC scores used across all PCA panels.", data = pca_scores_tbl %>% dplyr::select(dplyr::all_of(pca_keep))),
      list(sheet = "PCA_variance", title = "Variance explained by principal components", note = "Percent and cumulative variance for the displayed PCs.", data = pca_var_tbl)
    ),
    status = "Complete",
    notes = "One PCA score table supports all coloring panels; no duplicate table is required per panel.",
    source_files = c(pca_scores_path, pca_var_path)
  )
  record_output("DEP PCA figure", pca_file, "GENERATED")
} else {
  record_output("DEP PCA figure", NA_character_, "NOT GENERATED", "result/02_pca/pca_scores_with_metadata.csv and pca_variance_explained.csv")
}

###############################################################################
# Supplementary / Extended Data Fig. 2: demographics
###############################################################################

source_participants <- if (exists("sample_data")) tibble::as_tibble(sample_data) else tibble::as_tibble(dep_df)
group_col <- pick_col(source_participants, c("SampleGroup", "SampleGroup.csv", "SampleGroup.adat"))
sample_type_col <- pick_col(source_participants, c("SampleType", "SampleType.csv", "SampleType.adat"))

participants_demog <- source_participants
if (!is.na(sample_type_col)) participants_demog <- participants_demog %>% dplyr::filter(.data[[sample_type_col]] == "Sample")
if (!is.na(group_col)) participants_demog <- participants_demog %>% dplyr::filter(.data[[group_col]] %in% c("CN", "AD"))

keep_demog <- c(group_col, "Age", "Sex", "APOE4_carrier")
keep_demog <- unique(keep_demog[!is.na(keep_demog) & keep_demog %in% names(participants_demog)])
participants_demog <- participants_demog %>% dplyr::select(dplyr::all_of(keep_demog)) %>% dplyr::mutate(Record = sprintf("D%04d", dplyr::row_number()), .before = 1)
if (!is.na(group_col) && group_col != "SampleGroup") names(participants_demog)[names(participants_demog) == group_col] <- "SampleGroup"

demog_summary <- participants_demog %>%
  dplyr::group_by(SampleGroup) %>%
  dplyr::summarise(
    n = dplyr::n(),
    age_mean = mean(safe_numeric(Age), na.rm = TRUE),
    age_sd = stats::sd(safe_numeric(Age), na.rm = TRUE),
    female_n = sum(tolower(as.character(Sex)) %in% c("f", "female", "1"), na.rm = TRUE),
    .groups = "drop"
  )
if ("APOE4_carrier" %in% names(participants_demog)) {
  apoe_counts <- participants_demog %>%
    dplyr::group_by(SampleGroup) %>%
    dplyr::summarise(APOE4_carrier_n = sum(safe_numeric(APOE4_carrier) == 1, na.rm = TRUE), .groups = "drop")
  demog_summary <- demog_summary %>% dplyr::left_join(apoe_counts, by = "SampleGroup")
} else {
  demog_summary$APOE4_carrier_n <- NA_integer_
}

demog_file <- save_source_workbook(
  "SourceData_DEP_DemographicDistributions.xlsx",
  "Supplementary / Extended Data figure — APOE, age and sex distributions",
  list(
    list(sheet = "Participant_data", title = "Participant-level demographic source data", note = "Rows provide the values used to construct the age, sex and APOE distributions.", data = participants_demog),
    list(sheet = "Group_summary", title = "Diagnostic-group summary", note = "Convenience summary; participant-level values remain the source of truth.", data = demog_summary)
  ),
  status = "Complete",
  notes = "Participant-level age, sex and APOE ε4 carrier-status values are provided without direct or stable study identifiers.",
  source_files = c(workspace_file)
)
record_output("DEP demographic distributions", demog_file, "GENERATED")

###############################################################################
# Supplementary / Extended Data Fig. 3: P-value calibration
###############################################################################

pval_source <- primary_full %>%
  dplyr::transmute(
    Protein_Name, EntrezGeneSymbol, AptName, logFC,
    P.Value = safe_numeric(P.Value),
    adj.P.Val = safe_numeric(adj.P.Val),
    minus_log10_P = -log10(P.Value)
  )
qq_source <- pval_source %>%
  dplyr::filter(is.finite(P.Value), P.Value > 0) %>%
  dplyr::arrange(P.Value) %>%
  dplyr::mutate(
    rank = dplyr::row_number(),
    expected_P = (rank - 0.5) / dplyr::n(),
    observed_minus_log10_P = -log10(P.Value),
    expected_minus_log10_P = -log10(expected_P)
  )

pval_file <- save_source_workbook(
  "SourceData_DEP_PvalueCalibration.xlsx",
  "Supplementary / Extended Data figure — P-value distribution and QQ plot",
  list(
    list(sheet = "Pvalue_distribution", title = "P-value distribution source data", note = "One row per fixed gene–SOMAmer pair.", data = pval_source),
    list(sheet = "QQ_plot", title = "QQ plot source data", note = "Expected P values use (rank − 0.5)/n.", data = qq_source)
  ),
  status = "Complete",
  notes = "Both raw P values and −log10 transformed values are supplied.",
  source_files = c(primary_full_path)
)
record_output("DEP P-value calibration", pval_file, "GENERATED")

###############################################################################
# Supplementary / Extended Data Fig. 4: representative proteins and correlations
###############################################################################

boxplot_proteins <- c(
  "AMBRA1", "ANKRD13C", "C3", "CAPN1", "CPLX2",
  "ENPP2", "GJB3", "LRRN1", "MENT", "MMP7",
  "NEFL", "ODC1", "PLCH1", "PODXL", "PUM2",
  "RNPEP", "S100A13", "SMOC1", "SPC25", "USP5"
)
boxplot_map <- DEP_gene %>%
  dplyr::filter(EntrezGeneSymbol %in% boxplot_proteins) %>%
  dplyr::mutate(order = match(EntrezGeneSymbol, boxplot_proteins)) %>%
  dplyr::arrange(order) %>%
  dplyr::select(Protein_Name, EntrezGeneSymbol, AptName)

boxplot_source <- purrr::map_dfr(seq_len(nrow(boxplot_map)), function(i) {
  p <- boxplot_map[i, , drop = FALSE]
  if (!p$AptName[[1]] %in% names(dep_df)) return(tibble::tibble())
  tibble::tibble(
    SampleGroup = as.character(dep_df$SampleGroup),
    Protein = p$EntrezGeneSymbol[[1]],
    Protein_Name = p$Protein_Name[[1]],
    AptName = p$AptName[[1]],
    log2_RFU = safe_log2(dep_df[[p$AptName[[1]]]])
  ) %>%
    dplyr::left_join(primary_full %>% dplyr::select(EntrezGeneSymbol, AptName, logFC, adj.P.Val, Direction), by = c("Protein" = "EntrezGeneSymbol", "AptName"))
})
boxplot_source <- boxplot_source %>% dplyr::mutate(Record = sprintf("B%06d", dplyr::row_number()), .before = 1)

heatmap_proteins <- primary_full %>%
  dplyr::filter(adj.P.Val < MAIN_FDR) %>%
  dplyr::group_by(Direction) %>%
  dplyr::arrange(adj.P.Val, .by_group = TRUE) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(Direction), adj.P.Val) %>%
  dplyr::select(Protein_Name, EntrezGeneSymbol, AptName, logFC, adj.P.Val, Direction)

trait_candidates_s4 <- list(
  Diagnosis = c("SampleGroup"),
  Age = c("Age"),
  Sex = c("Sex"),
  Education = c("Education"),
  CDR_global = c("cdr_global", "CDR_global"),
  MMSE = c("mmse_total", "MMSE"),
  NPI = c("NPI", "NPI_Q"),
  p_tau181 = c("p.tau181", "p-tau181", "p_tau181"),
  p_tau217 = c("p.tau217", "p-tau217", "p_tau217"),
  NfL = c("NfL", "NFL")
)
trait_map_s4 <- tibble::tibble(
  trait = names(trait_candidates_s4),
  column = vapply(trait_candidates_s4, function(x) pick_col(dep_df, x), character(1))
)
heatmap_corr <- spearman_long(dep_df, heatmap_proteins, trait_map_s4)

protein_trait_file <- save_source_workbook(
  "SourceData_DEP_RepresentativeProteins.xlsx",
  "Supplementary / Extended Data figure — representative proteins and protein–trait landscape",
  list(
    list(sheet = "Boxplot_participant_data", title = "Panel a participant-level protein abundance", note = "Non-residualized log2 RFU for the 20 displayed proteins.", data = boxplot_source),
    list(sheet = "Boxplot_protein_map", title = "Proteins displayed in panel a", note = "Fixed primary gene–SOMAmer mapping and primary effect statistics.", data = boxplot_map %>% dplyr::left_join(primary_full, by = c("EntrezGeneSymbol", "AptName"))),
    list(sheet = "Heatmap_protein_selection", title = "Proteins displayed in panel b", note = "Top 25 higher- and top 25 lower-abundance primary DEPs ranked by BH-FDR.", data = heatmap_proteins),
    list(sheet = "Heatmap_correlations", title = "Panel b protein–trait correlations", note = "Spearman correlations with BH adjustment across all displayed protein–trait pairs.", data = heatmap_corr)
  ),
  status = "Complete",
  notes = "The source data recreate the participant-level boxplots and the extended correlation heatmap.",
  source_files = c(workspace_file, primary_full_path)
)
record_output("DEP representative proteins", protein_trait_file, "GENERATED")

###############################################################################
# Supplementary / Extended Data Fig. 5: internal sensitivity
###############################################################################

apoe_full_path <- first_existing_file(c(file.path(analysis_root, "04_sensitivity", "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv")))
atn_full_path <- first_existing_file(c(file.path(analysis_root, "04_sensitivity", "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv")))
cdr_full_path <- first_existing_file(c(file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv")))
apoe_full <- read_required_csv(apoe_full_path, "APOE-adjusted full fixed-map table")
atn_full <- read_required_csv(atn_full_path, "AT(N)-adjusted full fixed-map table")
cdr_full <- read_required_csv(cdr_full_path, "AD-only CDR-SB full fixed-map table")
assert_fixed_map(apoe_full, "APOE-adjusted table")
assert_fixed_map(atn_full, "AT(N)-adjusted table")
assert_fixed_map(cdr_full, "AD-only CDR-SB table")

apoe_summary <- read_optional_csv(c(file.path(analysis_root, "04_sensitivity", "apoe", "equal_sample_baseline", "APOE_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv")))
atn_summary <- read_optional_csv(c(file.path(analysis_root, "04_sensitivity", "atn_adjusted", "equal_sample_baseline", "ATN_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv")))
cdr_counts <- read_optional_csv(c(file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_counts_gene_collapsed.csv")))

sens_file <- save_source_workbook(
  "SourceData_DEP_InternalSensitivity.xlsx",
  "Supplementary / Extended Data figure — APOE, AT(N) and CDR-SB sensitivity",
  list(
    list(sheet = "APOE_volcano", title = "Panel a APOE-adjusted differential abundance", note = "One row per fixed gene–SOMAmer pair.", data = apoe_full),
    list(sheet = "APOE_comparison", title = "Panel b primary versus APOE-adjusted effects", note = "Full fixed-map effect comparison.", data = apoe_comparison),
    list(sheet = "ATN_volcano", title = "Panel c AT(N)-adjusted differential abundance", note = "One row per fixed gene–SOMAmer pair.", data = atn_full),
    list(sheet = "ATN_comparison", title = "Panel d primary versus AT(N)-adjusted effects", note = "Full fixed-map effect comparison.", data = atn_comparison),
    list(sheet = "CDRSB_volcano", title = "Panel e within-AD CDR-SB associations", note = "One row per fixed gene–SOMAmer pair.", data = cdr_full),
    list(sheet = "CDRSB_comparison", title = "Panel f primary diagnostic versus within-AD CDR-SB effects", note = "Full fixed-map effect comparison.", data = cdr_alignment),
    list(sheet = "Model_summaries", title = "Sensitivity model summaries", note = "Equal-sample summaries and within-AD significance counts.", data = dplyr::bind_rows(
      if (!is.null(apoe_summary)) apoe_summary %>% dplyr::mutate(source = "APOE") else NULL,
      if (!is.null(atn_summary)) atn_summary %>% dplyr::mutate(source = "ATN") else NULL,
      if (!is.null(cdr_counts)) cdr_counts %>% dplyr::mutate(source = "CDRSB") else NULL
    ))
  ),
  status = "Complete",
  notes = "All volcano and effect-comparison panels use the fixed 9,638 gene–SOMAmer map.",
  source_files = c(apoe_full_path, atn_full_path, cdr_full_path,
    file.path(analysis_root, "04_sensitivity", "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"),
    atn_comparison_path,
    file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"))
)
record_output("DEP internal sensitivity", sens_file, "GENERATED")

###############################################################################
# Supplementary / Extended Data Fig. 6: recruitment-context robustness
###############################################################################

loso_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "site_robustness", "loso", "tables", "LOSO_summary_metrics_FIXED_PRIMARY_MAP.csv")))
balanced_summary_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_summary_metrics_FIXED_PRIMARY_MAP.csv")))
balanced_protein_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_protein_stability_FIXED_PRIMARY_MAP.csv")))
country_meta_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "country_meta", "tables", "country_meta_analysis_results_FIXED_PRIMARY_MAP.csv")))
country_counts_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "country_loco", "tables", "country_group_counts.csv")))
site_counts_path <- first_existing_file(c(file.path(analysis_root, "06_robustness", "site_robustness", "site_group_counts.csv")))

loso <- read_required_csv(loso_path, "LOSO summary")
balanced_summary <- read_required_csv(balanced_summary_path, "balanced-resampling iteration summary")
balanced_protein <- read_required_csv(balanced_protein_path, "balanced-resampling protein stability")
country_meta <- read_required_csv(country_meta_path, "country meta-analysis")
assert_fixed_map(balanced_protein, "Balanced-resampling protein stability")
assert_fixed_map(country_meta, "Country meta-analysis")

robust_file <- save_source_workbook(
  "SourceData_DEP_RecruitmentRobustness.xlsx",
  "Supplementary / Extended Data figure — recruitment-context robustness",
  list(
    list(sheet = "LOCO_summary", title = "Panels a-b LOCO robustness", note = "Summary metrics for each excluded country.", data = loco_summary),
    list(sheet = "LOSO_summary", title = "Panels c-d LOSO robustness", note = "Summary metrics for each eligible excluded recruitment site.", data = loso),
    list(sheet = "Balanced_iterations", title = "Panels e-g balanced-resampling distributions", note = "One row per resampling iteration; these values reconstruct the histograms.", data = balanced_summary),
    list(sheet = "Balanced_protein_stability", title = "Protein-level balanced-resampling stability", note = "Complete fixed-map protein-level stability output.", data = balanced_protein),
    list(sheet = "Country_meta", title = "Country-level meta-analysis", note = "Complete fixed-map meta-analytic results and heterogeneity statistics.", data = country_meta),
    list(sheet = "Country_counts", title = "Country diagnostic-group counts", note = "Counts underlying country eligibility and balancing.", data = read_optional_csv(country_counts_path)),
    list(sheet = "Site_counts", title = "Site diagnostic-group counts", note = "Counts underlying site eligibility.", data = read_optional_csv(site_counts_path))
  ),
  status = "Complete",
  notes = "The balanced-resampling iteration table is the direct source for the three histogram panels.",
  source_files = c(loco_summary_path, loso_path, balanced_summary_path, balanced_protein_path, country_meta_path, country_counts_path, site_counts_path)
)
record_output("DEP recruitment robustness", robust_file, "GENERATED")

###############################################################################
# p-tau217-enriched threshold-sweep figure (current Results cite separately)
###############################################################################

ptau_root <- file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables")
ptau_summary_path <- first_existing_file(c(file.path(ptau_root, "p_tau217_enriched_threshold_sweep_summary_v3.csv")))
ptau_comparison_path <- first_existing_file(c(file.path(ptau_root, "primary_vs_p_tau217_enriched_threshold_sweep_all_comparisons_v3.csv")))
ptau_counts_path <- first_existing_file(c(file.path(ptau_root, "p_tau217_enriched_threshold_sweep_sample_counts_v3.csv")))
ptau_thresholds_path <- first_existing_file(c(file.path(ptau_root, "cohort_defined_biomarker_thresholds_v3.csv")))
ptau_membership_path <- first_existing_file(c(file.path(ptau_root, "p_tau217_threshold_membership_audit.csv")))

ptau_summary <- read_required_csv(ptau_summary_path, "p-tau217 threshold summary")
ptau_comparison <- read_required_csv(ptau_comparison_path, "p-tau217 all-comparison table")
ptau_counts <- read_required_csv(ptau_counts_path, "p-tau217 sample counts")
ptau_thresholds <- read_required_csv(ptau_thresholds_path, "cohort-defined thresholds")
ptau_membership <- read_optional_csv(ptau_membership_path)

ptau_p90 <- ptau_comparison %>% dplyr::filter(comparison == "p_tau217_enriched_CNp90")
ptau_p90_stringent <- ptau_comparison %>% dplyr::filter(comparison == "p_tau217_enriched_low_abeta42_40_CNp90")
if (nrow(ptau_p90) != EXPECTED_MAP_N) warning("p-tau217 p90 comparison does not contain 9,638 rows.")
if (nrow(ptau_p90_stringent) != EXPECTED_MAP_N) warning("Stringent p90 comparison does not contain 9,638 rows.")

ptau_file <- save_source_workbook(
  "SourceData_DEP_pTau217ThresholdSweep.xlsx",
  "Supplementary / Extended Data figure — p-tau217-enriched threshold sensitivity",
  list(
    list(sheet = "Threshold_definitions", title = "Cohort-defined p-tau217 and Aβ42/40 thresholds", note = "These thresholds are exploratory enrichment anchors, not clinical diagnostic cutoffs.", data = ptau_thresholds),
    list(sheet = "Threshold_summary", title = "Threshold-sweep summary", note = "Sample sizes, DEP counts, effect concordance and preservation metrics for p80, p90 and p95 models.", data = ptau_summary),
    list(sheet = "Sample_counts", title = "Threshold-sweep sample counts", note = "Group counts used in each contrast.", data = ptau_counts),
    list(sheet = "Threshold_membership", title = "Threshold membership audit", note = "Counts above and below each cohort-defined p-tau217 threshold.", data = ptau_membership),
    list(sheet = "P90_comparison", title = "Primary versus p-tau217-enriched CN p90 effects", note = "Main-text threshold sensitivity; one row per fixed gene–SOMAmer pair.", data = ptau_p90),
    list(sheet = "P90_FDR005", title = "FDR-significant proteins in the p90 contrast", note = "Subset of the p90 comparison with secondary FDR < 0.05.", data = ptau_p90 %>% dplyr::filter(secondary_fdr005 %in% TRUE)),
    list(sheet = "P90_stringent_comparison", title = "Primary versus p90 plus low Aβ42/40 contrast", note = "More selective secondary contrast using the CN p10 Aβ42/40 threshold.", data = ptau_p90_stringent)
  ),
  status = "Complete",
  notes = "Figure numbering should be finalized after the Extended Data figures are consolidated to the journal limit.",
  source_files = c(ptau_summary_path, ptau_comparison_path, ptau_counts_path, ptau_thresholds_path, ptau_membership_path)
)
record_output("DEP p-tau217 threshold sweep", ptau_file, "GENERATED")

###############################################################################
# Global source-data index and audit
###############################################################################

source_index <- dplyr::bind_rows(output_registry)
index_csv <- file.path(OUT_ROOT, "DEP_SourceData_Index.csv")
readr::write_csv(source_index, index_csv)

index_wb <- openxlsx::createWorkbook()
styles <- make_styles()
write_source_sheet(index_wb, "Index", "DEP Source Data index", "One workbook is generated per quantitative figure. Figure numbering may be updated during final Extended Data consolidation.", source_index, styles)
openxlsx::saveWorkbook(index_wb, file.path(OUT_ROOT, "DEP_SourceData_Index.xlsx"), overwrite = TRUE)

manifest <- tibble::tibble(
  generated_at = as.character(Sys.time()),
  project_root = ".",
  workspace = relative_source_paths(workspace_file),
  expected_fixed_map_n = EXPECTED_MAP_N,
  n_generated = sum(source_index$Status == "GENERATED"),
  n_not_generated = sum(source_index$Status != "GENERATED")
)
readr::write_csv(manifest, file.path(OUT_ROOT, "DEP_SourceData_RunManifest.csv"))
writeLines(capture.output(utils::sessionInfo()), file.path(OUT_ROOT, "DEP_SourceData_sessionInfo.txt"))

message("\n============================================================")
message("Nature Aging DEP Source Data package generated in:")
message(OUT_ROOT)
message("Generated workbooks: ", sum(source_index$Status == "GENERATED"))
message("Not generated: ", sum(source_index$Status != "GENERATED"))
message("============================================================\n")
print(source_index)
###############################################################################
# END
###############################################################################

