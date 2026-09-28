###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 01. Prepare outcome-independent WGCNA input
# Requires: canonical DEP workspace and gene-level results
# Produces: gene-collapsed expression, metadata, annotation and QC
# Data policy: participant-level inputs and intermediate outputs remain local.
###############################################################################

rm(list = ls())


# -----------------------------------------------------------------------------
# Repository configuration
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the project environment with renv::restore().", call. = FALSE)
}
source(file.path(project_root, "R", "wgcna_bootstrap.R"), local = FALSE)
WGCNA_CONFIG <- wgcna_load_config(project_root)

# =============================================================================
# 0) CONFIGURATION
# =============================================================================

# Explicit project roots.
#
# WGCNA project:
WGCNA_PROJECT_ROOT <- WGCNA_CONFIG$project_root

# DEP project containing the cleaned DEP workspace and primary gene table:
DEP_PROJECT_ROOT <- WGCNA_CONFIG$dep_project_root

# Keep PROJECT_ROOT as an alias for the WGCNA project to preserve the original
# script structure and downstream QC fields.
PROJECT_ROOT <- WGCNA_PROJECT_ROOT

# Candidate DEP output locations. The first existing path is used.
DEP_WORKSPACE_CANDIDATES <- c(
  file.path(WGCNA_CONFIG$dep_result_root, "workspace", "analysis_workspace.RData")
)

DEP_GENE_TABLE_CANDIDATES <- c(
  file.path(WGCNA_CONFIG$dep_result_root, "03_dep", "gene_collapsed",
            "AD_vs_CN_full_limma_results_gene_collapsed.csv")
)

# WGCNA output folder for this script.
OUTDIR <- file.path(WGCNA_CONFIG$result_root, "01_input"
)

MAIN_GROUPS <- c("CN", "AD")

CORE_METADATA <- c(
  "SampleId",
  "SampleGroup",
  "Age",
  "Sex",
  "Country",
  "Education"
)

# Expression scale:
#   "auto"     = detect raw RFU versus already-log2 values;
#   "raw_rfu"  = force log2 transformation;
#   "log2_rfu" = do not transform again.
EXPRESSION_SCALE_MODE <- "auto"

REMOVE_ZERO_VARIANCE_FEATURES <- TRUE
MIN_NON_MISSING_PER_FEATURE <- 10

# The WGCNA gene universe should remain identical to the primary DEP gene
# universe. Set FALSE only for exploratory use.
REQUIRE_EXACT_PRIMARY_GENE_UNIVERSE <- TRUE

# Deterministic scale audit uses a subset of features to avoid unnecessary
# memory duplication. Selection does not use diagnosis.
N_FEATURES_FOR_SCALE_AUDIT <- 500

# =============================================================================
# 1) PACKAGES AND HELPERS
# =============================================================================

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "purrr",
  "stringr"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A", "nan")] <- NA_character_
  x
}

safe_log2_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x <= 0] <- NA_real_
  log2(x)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_mad <- function(x) {
  x <- safe_numeric(x)
  if (sum(is.finite(x)) < 2) return(NA_real_)
  stats::mad(x, center = stats::median(x, na.rm = TRUE), na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- safe_numeric(x)
  if (sum(is.finite(x)) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

stop_if_missing_file <- function(file, label, attempted = character(0)) {
  if (is.na(file) || !file.exists(file)) {
    stop(
      label,
      " not found.\nAttempted paths:\n",
      paste(attempted, collapse = "\n"),
      call. = FALSE
    )
  }
}

pick_first_col <- function(df, candidates, required = FALSE, label = NULL) {
  hit <- candidates[candidates %in% names(df)][1]

  if (length(hit) == 0 || is.na(hit)) {
    if (required) {
      stop(
        ifelse(is.null(label), "Required column", label),
        " was not found. Candidates: ",
        paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
    return(NA_character_)
  }

  hit
}

detect_expression_scale <- function(df, n_features = 500) {
  if (ncol(df) == 0) {
    stop("No expression features were available for scale detection.", call. = FALSE)
  }

  n_take <- min(n_features, ncol(df))
  idx <- unique(round(seq(1, ncol(df), length.out = n_take)))
  sampled <- df[, idx, drop = FALSE]

  values <- unlist(
    lapply(sampled, safe_numeric),
    use.names = FALSE
  )

  values <- values[is.finite(values)]

  if (length(values) == 0) {
    stop("Expression scale could not be detected because all audited values were missing.", call. = FALSE)
  }

  q <- stats::quantile(
    values,
    probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
    na.rm = TRUE,
    names = FALSE
  )

  q01 <- q[[1]]
  q05 <- q[[2]]
  q50 <- q[[3]]
  q95 <- q[[4]]
  q99 <- q[[5]]

  # SOMAscan raw RFU values are typically far above the log2 range.
  if (q99 > 100 || q50 > 40) {
    detected_scale <- "raw_rfu"
    action <- "apply_log2"
    confidence <- "high"
  } else if (q99 <= 40 && q50 <= 30) {
    detected_scale <- "log2_rfu"
    action <- "keep_as_is"
    confidence <- "high"
  } else {
    detected_scale <- "ambiguous"
    action <- "stop_and_review"
    confidence <- "low"
  }

  audit <- tibble::tibble(
    n_features_audited = length(idx),
    n_values_audited = length(values),
    minimum = min(values, na.rm = TRUE),
    q01 = q01,
    q05 = q05,
    median = q50,
    q95 = q95,
    q99 = q99,
    maximum = max(values, na.rm = TRUE),
    nonpositive_proportion = mean(values <= 0, na.rm = TRUE),
    detected_scale = detected_scale,
    recommended_action = action,
    confidence = confidence
  )

  list(
    detected_scale = detected_scale,
    action = action,
    audit = audit
  )
}

transform_expression_df <- function(df, scale_mode, scale_detection) {
  allowed <- c("auto", "raw_rfu", "log2_rfu")

  if (!scale_mode %in% allowed) {
    stop(
      "EXPRESSION_SCALE_MODE must be one of: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  resolved_mode <- scale_mode

  if (scale_mode == "auto") {
    resolved_mode <- scale_detection$detected_scale

    if (resolved_mode == "ambiguous") {
      stop(
        "Expression scale was ambiguous. Review qc/expression_scale_audit.csv ",
        "and set EXPRESSION_SCALE_MODE explicitly to 'raw_rfu' or 'log2_rfu'.",
        call. = FALSE
      )
    }
  }

  if (resolved_mode == "raw_rfu") {
    transformed <- as.data.frame(
      lapply(df, safe_log2_vector),
      check.names = FALSE
    )
    transformation <- "log2 applied once"
  } else {
    transformed <- as.data.frame(
      lapply(df, safe_numeric),
      check.names = FALSE
    )
    transformation <- "input retained; no additional log2"
  }

  list(
    data = transformed,
    resolved_mode = resolved_mode,
    transformation = transformation
  )
}

# =============================================================================
# 2) CREATE OUTPUT STRUCTURE
# =============================================================================

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "qc"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 3) LOCATE AND CHECK INPUTS
# =============================================================================

DEP_WORKSPACE_FILE <- first_existing_file(DEP_WORKSPACE_CANDIDATES)
DEP_GENE_TABLE <- first_existing_file(DEP_GENE_TABLE_CANDIDATES)

stop_if_missing_file(
  DEP_WORKSPACE_FILE,
  "DEP workspace",
  DEP_WORKSPACE_CANDIDATES
)

stop_if_missing_file(
  DEP_GENE_TABLE,
  "DEP gene-collapsed table",
  DEP_GENE_TABLE_CANDIDATES
)

message("Input files detected.")
message("WGCNA project root: ", WGCNA_PROJECT_ROOT)
message("DEP project root: ", DEP_PROJECT_ROOT)
message("DEP workspace: ", DEP_WORKSPACE_FILE)
message("DEP gene table: ", DEP_GENE_TABLE)
message("WGCNA output directory: ", OUTDIR)

# =============================================================================
# 4) LOAD DEP WORKSPACE AND PRIMARY GENE TABLE
# =============================================================================

load(DEP_WORKSPACE_FILE)

required_objects <- c(
  "sample_data",
  "annot_tbl",
  "protein_universe"
)

missing_objects <- required_objects[
  !vapply(required_objects, function(x) exists(x, envir = environment()), logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "The DEP workspace is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# Aptamer-level DEP output is needed only to annotate the independently selected
# WGCNA representative with its own model statistics.
if (exists("DEP_aptamer")) {
  dep_aptamer_tbl <- as.data.frame(DEP_aptamer)
} else if (exists("DEP") && "AptName" %in% names(DEP)) {
  dep_aptamer_tbl <- as.data.frame(DEP)
} else {
  stop(
    "The DEP workspace must contain DEP_aptamer or an aptamer-level DEP object named DEP.",
    call. = FALSE
  )
}

DEP_gene_from_file <- readr::read_csv(
  DEP_GENE_TABLE,
  show_col_types = FALSE
)

required_dep_cols <- c(
  "EntrezGeneSymbol",
  "AptName",
  "SeqId",
  "Protein_Name",
  "logFC",
  "adj.P.Val"
)

missing_dep_cols <- setdiff(required_dep_cols, names(DEP_gene_from_file))

if (length(missing_dep_cols) > 0) {
  stop(
    "DEP gene table is missing required columns: ",
    paste(missing_dep_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!"SampleId" %in% names(sample_data)) {
  stop("sample_data must contain SampleId.", call. = FALSE)
}

if (anyDuplicated(sample_data$SampleId) > 0) {
  warning(
    "Duplicated SampleId values were detected in sample_data. ",
    "The first row per SampleId will be retained after analytic filtering."
  )
}

# Standardize primary gene map.
primary_gene_map <- DEP_gene_from_file %>%
  dplyr::mutate(
    EntrezGeneSymbol = clean_text_na(EntrezGeneSymbol),
    Primary_DEP_AptName = clean_text_na(AptName),
    Primary_DEP_SeqId = as.character(SeqId),
    Primary_DEP_Protein_Name = as.character(Protein_Name),
    Primary_DEP_logFC = safe_numeric(logFC),
    Primary_DEP_P.Value = if ("P.Value" %in% names(.)) safe_numeric(P.Value) else NA_real_,
    Primary_DEP_adj.P.Val = safe_numeric(adj.P.Val),
    Primary_DEP_type = if ("type" %in% names(.)) as.character(type) else NA_character_,
    Primary_DEP_Direction = if ("Direction" %in% names(.)) as.character(Direction) else NA_character_
  ) %>%
  dplyr::filter(
    !is.na(EntrezGeneSymbol),
    !is.na(Primary_DEP_AptName)
  ) %>%
  dplyr::distinct(EntrezGeneSymbol, .keep_all = TRUE) %>%
  dplyr::select(
    EntrezGeneSymbol,
    dplyr::starts_with("Primary_DEP_")
  )

if (anyDuplicated(primary_gene_map$EntrezGeneSymbol) > 0) {
  stop("Primary DEP gene map contains duplicated gene symbols.", call. = FALSE)
}

if (anyDuplicated(primary_gene_map$Primary_DEP_AptName) > 0) {
  warning(
    "The primary DEP map contains duplicated SOMAmer identifiers across genes. ",
    "This will be documented in the comparison audit."
  )
}

primary_gene_symbols <- primary_gene_map$EntrezGeneSymbol
n_primary_genes <- length(primary_gene_symbols)

message("Primary DEP gene universe: ", n_primary_genes, " genes.")

# =============================================================================
# 5) BUILD WGCNA SAMPLE METADATA
# =============================================================================

metadata_cols_optional <- c(
  "SampleType",
  "ApoE", "APOE_group", "APOE4_carrier",
  "cdr_global", "cdr_boxscore", "mmse_total", "udsfaq_total",
  "NPI", "Mini.SEA", "Mini-SEA", "Mini_SEA",
  "T.ADLQ", "T-ADLQ", "T_ADLQ",
  "p.tau181", "p-tau181", "p_tau181",
  "p.tau217", "p-tau217", "p_tau217",
  "NfL",
  "ratio.AB42.40", "ratio AB42/40", "ratio_AB42_40",
  "Site", "site", "Center", "center",
  "Cohort", "cohort",
  "PlateId", "PlateId.adat"
)

metadata_cols <- unique(
  c(CORE_METADATA, intersect(metadata_cols_optional, names(sample_data)))
)

missing_core <- setdiff(CORE_METADATA, names(sample_data))

if (length(missing_core) > 0) {
  stop(
    "sample_data is missing core metadata columns: ",
    paste(missing_core, collapse = ", "),
    call. = FALSE
  )
}

wgcna_metadata <- sample_data %>%
  dplyr::filter(SampleGroup %in% MAIN_GROUPS) %>%
  {
    if ("SampleType" %in% names(.)) {
      dplyr::filter(., is.na(SampleType) | SampleType == "Sample")
    } else {
      .
    }
  } %>%
  dplyr::select(dplyr::all_of(metadata_cols)) %>%
  dplyr::mutate(
    SampleId = as.character(SampleId),
    SampleGroup = as.character(SampleGroup),
    Sex = as.character(Sex),
    Country = as.character(Country),
    Age = safe_numeric(Age),
    Education = safe_numeric(Education)
  ) %>%
  dplyr::filter(
    stats::complete.cases(
      dplyr::across(dplyr::all_of(CORE_METADATA))
    )
  ) %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

wgcna_metadata <- wgcna_metadata %>%
  dplyr::mutate(
    Mini_SEA = dplyr::coalesce(
      if ("Mini_SEA" %in% names(.)) safe_numeric(.data[["Mini_SEA"]]) else NA_real_,
      if ("Mini.SEA" %in% names(.)) safe_numeric(.data[["Mini.SEA"]]) else NA_real_,
      if ("Mini-SEA" %in% names(.)) safe_numeric(.data[["Mini-SEA"]]) else NA_real_
    ),
    T_ADLQ = dplyr::coalesce(
      if ("T_ADLQ" %in% names(.)) safe_numeric(.data[["T_ADLQ"]]) else NA_real_,
      if ("T.ADLQ" %in% names(.)) safe_numeric(.data[["T.ADLQ"]]) else NA_real_,
      if ("T-ADLQ" %in% names(.)) safe_numeric(.data[["T-ADLQ"]]) else NA_real_
    ),
    p_tau181 = dplyr::coalesce(
      if ("p_tau181" %in% names(.)) safe_numeric(.data[["p_tau181"]]) else NA_real_,
      if ("p.tau181" %in% names(.)) safe_numeric(.data[["p.tau181"]]) else NA_real_,
      if ("p-tau181" %in% names(.)) safe_numeric(.data[["p-tau181"]]) else NA_real_
    ),
    p_tau217 = dplyr::coalesce(
      if ("p_tau217" %in% names(.)) safe_numeric(.data[["p_tau217"]]) else NA_real_,
      if ("p.tau217" %in% names(.)) safe_numeric(.data[["p.tau217"]]) else NA_real_,
      if ("p-tau217" %in% names(.)) safe_numeric(.data[["p-tau217"]]) else NA_real_
    ),
    ratio_AB42_40 = dplyr::coalesce(
      if ("ratio_AB42_40" %in% names(.)) safe_numeric(.data[["ratio_AB42_40"]]) else NA_real_,
      if ("ratio.AB42.40" %in% names(.)) safe_numeric(.data[["ratio.AB42.40"]]) else NA_real_,
      if ("ratio AB42/40" %in% names(.)) safe_numeric(.data[["ratio AB42/40"]]) else NA_real_
    )
  )

if (!"APOE4_carrier" %in% names(wgcna_metadata) && "ApoE" %in% names(wgcna_metadata)) {
  wgcna_metadata <- wgcna_metadata %>%
    dplyr::mutate(
      ApoE = trimws(as.character(ApoE)),
      APOE4_carrier = dplyr::case_when(
        ApoE %in% c("e2/e4", "e3/e4", "e4/e4") ~ 1,
        ApoE %in% c("e2/e2", "e2/e3", "e3/e3") ~ 0,
        TRUE ~ NA_real_
      )
    )
}

if ("APOE4_carrier" %in% names(wgcna_metadata)) {
  wgcna_metadata <- wgcna_metadata %>%
    dplyr::mutate(APOE4_carrier = safe_numeric(APOE4_carrier))
}

if (nrow(wgcna_metadata) == 0) {
  stop("No samples remained after CN/AD and core metadata filtering.", call. = FALSE)
}

message(
  "WGCNA samples after CN/AD and core metadata filtering: ",
  nrow(wgcna_metadata)
)

# =============================================================================
# 6) BUILD OUTCOME-INDEPENDENT SOMAmer CANDIDATE TABLE
# =============================================================================

required_annot_cols <- c(
  "AptName",
  "EntrezGeneSymbol"
)

missing_annot_cols <- setdiff(required_annot_cols, names(annot_tbl))

if (length(missing_annot_cols) > 0) {
  stop(
    "annot_tbl is missing required columns: ",
    paste(missing_annot_cols, collapse = ", "),
    call. = FALSE
  )
}

candidate_annotation <- as.data.frame(annot_tbl) %>%
  dplyr::mutate(
    AptName = clean_text_na(AptName),
    EntrezGeneSymbol = clean_text_na(EntrezGeneSymbol),
    SeqId = if ("SeqId" %in% names(.)) as.character(SeqId) else NA_character_,
    Protein_Name = if ("Protein_Name" %in% names(.)) as.character(Protein_Name) else EntrezGeneSymbol,
    Organism = if ("Organism" %in% names(.)) as.character(Organism) else NA_character_,
    Type = if ("Type" %in% names(.)) as.character(Type) else NA_character_
  ) %>%
  dplyr::filter(
    !is.na(AptName),
    !is.na(EntrezGeneSymbol),
    EntrezGeneSymbol %in% primary_gene_symbols,
    AptName %in% names(sample_data),
    AptName %in% as.character(protein_universe)
  ) %>%
  dplyr::distinct(AptName, .keep_all = TRUE) %>%
  dplyr::select(
    -dplyr::any_of(c(
      "n_non_missing_primary",
      "missing_prop_primary",
      "median_log2_RFU_primary",
      "MAD_log2_RFU_primary",
      "SD_log2_RFU_primary",
      "eligible_for_selection",
      "selection_rank_within_gene",
      "n_eligible_candidates_for_gene",
      "Selected_for_gene_matrix",
      "Gene_matrix_column",
      "Selection_rule"
    ))
  )

if (nrow(candidate_annotation) == 0) {
  stop(
    "No eligible human protein SOMAmer candidates could be mapped to the primary gene universe.",
    call. = FALSE
  )
}

candidate_gene_coverage <- candidate_annotation %>%
  dplyr::distinct(EntrezGeneSymbol)

genes_without_candidate <- setdiff(
  primary_gene_symbols,
  candidate_gene_coverage$EntrezGeneSymbol
)

safe_write_csv(
  tibble::tibble(EntrezGeneSymbol = genes_without_candidate),
  file.path(OUTDIR, "qc", "primary_genes_without_candidate_somamer.csv")
)

if (
  REQUIRE_EXACT_PRIMARY_GENE_UNIVERSE &&
  length(genes_without_candidate) > 0
) {
  stop(
    "Outcome-independent candidates were unavailable for ",
    length(genes_without_candidate),
    " primary genes. See qc/primary_genes_without_candidate_somamer.csv.",
    call. = FALSE
  )
}

candidate_aptamers <- candidate_annotation$AptName

candidate_expr_df <- sample_data %>%
  dplyr::filter(SampleId %in% wgcna_metadata$SampleId) %>%
  dplyr::select(
    SampleId,
    dplyr::all_of(candidate_aptamers)
  ) %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

candidate_expr_df <- wgcna_metadata %>%
  dplyr::select(SampleId) %>%
  dplyr::left_join(
    candidate_expr_df,
    by = "SampleId"
  )

stopifnot(all(candidate_expr_df$SampleId == wgcna_metadata$SampleId))

candidate_expr_raw <- candidate_expr_df %>%
  dplyr::select(-SampleId) %>%
  as.data.frame(check.names = FALSE)

# =============================================================================
# 7) AUDIT AND RESOLVE EXPRESSION SCALE
# =============================================================================

scale_detection <- detect_expression_scale(
  candidate_expr_raw,
  n_features = N_FEATURES_FOR_SCALE_AUDIT
)

scale_transformation <- transform_expression_df(
  candidate_expr_raw,
  scale_mode = EXPRESSION_SCALE_MODE,
  scale_detection = scale_detection
)

candidate_expr_log2 <- scale_transformation$data

expression_scale_audit <- scale_detection$audit %>%
  dplyr::mutate(
    requested_mode = EXPRESSION_SCALE_MODE,
    resolved_mode = scale_transformation$resolved_mode,
    transformation_applied = scale_transformation$transformation
  )

safe_write_csv(
  expression_scale_audit,
  file.path(OUTDIR, "qc", "expression_scale_audit.csv")
)

message(
  "Expression scale resolved as: ",
  scale_transformation$resolved_mode
)
message(
  "Transformation: ",
  scale_transformation$transformation
)

# =============================================================================
# 8) OUTCOME-INDEPENDENT GENE-COLLAPSE MAP
# =============================================================================

aptamer_selection_qc <- tibble::tibble(
  AptName = colnames(candidate_expr_log2),
  n_non_missing = vapply(
    candidate_expr_log2,
    function(x) sum(is.finite(safe_numeric(x))),
    numeric(1)
  ),
  missing_prop = vapply(
    candidate_expr_log2,
    function(x) mean(!is.finite(safe_numeric(x))),
    numeric(1)
  ),
  median_log2_RFU = vapply(
    candidate_expr_log2,
    function(x) stats::median(safe_numeric(x), na.rm = TRUE),
    numeric(1)
  ),
  MAD_log2_RFU = vapply(
    candidate_expr_log2,
    safe_mad,
    numeric(1)
  ),
  SD_log2_RFU = vapply(
    candidate_expr_log2,
    safe_sd,
    numeric(1)
  )
) %>%
  dplyr::mutate(
    keep_non_missing = n_non_missing >= MIN_NON_MISSING_PER_FEATURE,
    keep_variance = is.finite(MAD_log2_RFU) & MAD_log2_RFU > 0,
    eligible_for_selection = keep_non_missing & keep_variance
  ) %>%
  dplyr::left_join(
    candidate_annotation,
    by = "AptName"
  ) %>%
  dplyr::group_by(EntrezGeneSymbol) %>%
  dplyr::arrange(
    missing_prop,
    dplyr::desc(MAD_log2_RFU),
    AptName,
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    selection_rank_within_gene = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

eligible_candidate_counts <- aptamer_selection_qc %>%
  dplyr::filter(eligible_for_selection) %>%
  dplyr::count(
    EntrezGeneSymbol,
    name = "n_eligible_candidates"
  )

genes_without_eligible_candidate <- setdiff(
  primary_gene_symbols,
  eligible_candidate_counts$EntrezGeneSymbol
)

safe_write_csv(
  tibble::tibble(
    EntrezGeneSymbol = genes_without_eligible_candidate
  ),
  file.path(OUTDIR, "qc", "primary_genes_without_eligible_candidate.csv")
)

if (
  REQUIRE_EXACT_PRIMARY_GENE_UNIVERSE &&
  length(genes_without_eligible_candidate) > 0
) {
  stop(
    "No eligible outcome-independent SOMAmer remained for ",
    length(genes_without_eligible_candidate),
    " primary genes after missingness and variance QC. ",
    "See qc/primary_genes_without_eligible_candidate.csv.",
    call. = FALSE
  )
}

gene_map <- aptamer_selection_qc %>%
  dplyr::filter(eligible_for_selection) %>%
  dplyr::group_by(EntrezGeneSymbol) %>%
  dplyr::arrange(
    missing_prop,
    dplyr::desc(MAD_log2_RFU),
    AptName,
    .by_group = TRUE
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    eligible_candidate_counts,
    by = "EntrezGeneSymbol"
  ) %>%
  dplyr::mutate(
    selection_rule = paste(
      "lowest missing_prop",
      "highest MAD_log2_RFU",
      "AptName tie-breaker",
      sep = " > "
    )
  )

if (anyDuplicated(gene_map$EntrezGeneSymbol) > 0) {
  stop("Outcome-independent gene map contains duplicated genes.", call. = FALSE)
}

if (anyDuplicated(gene_map$AptName) > 0) {
  stop(
    "Outcome-independent gene map contains duplicated SOMAmer identifiers.",
    call. = FALSE
  )
}

if (
  REQUIRE_EXACT_PRIMARY_GENE_UNIVERSE &&
  nrow(gene_map) != n_primary_genes
) {
  stop(
    "Outcome-independent gene map contains ",
    nrow(gene_map),
    " genes, but the primary DEP universe contains ",
    n_primary_genes,
    ".",
    call. = FALSE
  )
}

message(
  "Outcome-independent gene-collapsed features selected: ",
  nrow(gene_map)
)

# =============================================================================
# 9) ANNOTATE SELECTED WGCNA SOMAmers WITH THEIR OWN DEP STATISTICS
# =============================================================================

dep_aptamer_stats <- dep_aptamer_tbl %>%
  dplyr::mutate(
    AptName = clean_text_na(AptName),
    WGCNA_logFC = if ("logFC" %in% names(.)) safe_numeric(logFC) else NA_real_,
    WGCNA_P.Value = if ("P.Value" %in% names(.)) safe_numeric(P.Value) else NA_real_,
    WGCNA_adj.P.Val = if ("adj.P.Val" %in% names(.)) safe_numeric(adj.P.Val) else NA_real_,
    WGCNA_type = if ("type" %in% names(.)) as.character(type) else NA_character_,
    WGCNA_Direction = if ("Direction" %in% names(.)) as.character(Direction) else NA_character_
  ) %>%
  dplyr::filter(!is.na(AptName)) %>%
  dplyr::distinct(AptName, .keep_all = TRUE) %>%
  dplyr::select(
    AptName,
    WGCNA_logFC,
    WGCNA_P.Value,
    WGCNA_adj.P.Val,
    WGCNA_type,
    WGCNA_Direction
  )

gene_map <- gene_map %>%
  dplyr::left_join(
    dep_aptamer_stats,
    by = "AptName"
  ) %>%
  dplyr::left_join(
    primary_gene_map,
    by = "EntrezGeneSymbol"
  ) %>%
  dplyr::mutate(
    same_somamer_as_primary_DEP = AptName == Primary_DEP_AptName,
    logFC = WGCNA_logFC,
    P.Value = WGCNA_P.Value,
    adj.P.Val = WGCNA_adj.P.Val,
    type = WGCNA_type,
    Direction = WGCNA_Direction
  )

gene_map_comparison <- gene_map %>%
  dplyr::transmute(
    EntrezGeneSymbol,
    WGCNA_AptName = AptName,
    WGCNA_SeqId = SeqId,
    WGCNA_Protein_Name = Protein_Name,
    WGCNA_missing_prop = missing_prop,
    WGCNA_MAD_log2_RFU = MAD_log2_RFU,
    WGCNA_logFC,
    WGCNA_P.Value,
    WGCNA_adj.P.Val,
    Primary_DEP_AptName,
    Primary_DEP_SeqId,
    Primary_DEP_Protein_Name,
    Primary_DEP_logFC,
    Primary_DEP_P.Value,
    Primary_DEP_adj.P.Val,
    same_somamer_as_primary_DEP,
    n_eligible_candidates,
    selection_rule
  ) %>%
  dplyr::arrange(
    dplyr::desc(!same_somamer_as_primary_DEP),
    EntrezGeneSymbol
  )

n_same_map <- sum(
  gene_map_comparison$same_somamer_as_primary_DEP,
  na.rm = TRUE
)

n_changed_map <- sum(
  !gene_map_comparison$same_somamer_as_primary_DEP,
  na.rm = TRUE
)

gene_map_audit <- tibble::tibble(
  metric = c(
    "n_primary_DEP_genes",
    "n_outcome_independent_WGCNA_genes",
    "n_unique_WGCNA_SOMAmers",
    "n_genes_same_SOMAmer_as_primary_DEP",
    "n_genes_different_SOMAmer_from_primary_DEP",
    "proportion_same_SOMAmer_as_primary_DEP",
    "n_primary_genes_without_candidate",
    "n_primary_genes_without_eligible_candidate",
    "selection_used_diagnosis",
    "selection_used_logFC",
    "selection_used_P_value_or_FDR"
  ),
  value = c(
    n_primary_genes,
    nrow(gene_map),
    dplyr::n_distinct(gene_map$AptName),
    n_same_map,
    n_changed_map,
    n_same_map / nrow(gene_map),
    length(genes_without_candidate),
    length(genes_without_eligible_candidate),
    FALSE,
    FALSE,
    FALSE
  )
)

safe_write_csv(
  aptamer_selection_qc,
  file.path(OUTDIR, "tables", "aptamer_selection_qc.csv")
)

safe_write_csv(
  gene_map,
  file.path(OUTDIR, "tables", "outcome_independent_gene_somamer_map.csv")
)

safe_write_csv(
  gene_map_comparison,
  file.path(OUTDIR, "tables", "gene_map_comparison_vs_primary_DEP.csv")
)

safe_write_csv(
  gene_map_audit,
  file.path(OUTDIR, "qc", "gene_map_audit.csv")
)

# =============================================================================
# 10) BUILD GENE-COLLAPSED EXPRESSION MATRIX
# =============================================================================

expr_cols <- gene_map$AptName

# Reuse the already transformed candidate matrix to avoid applying log2 twice.
selected_col_idx <- match(
  expr_cols,
  colnames(candidate_expr_log2)
)

if (any(is.na(selected_col_idx))) {
  stop(
    "Some selected SOMAmers were not found in the transformed candidate matrix.",
    call. = FALSE
  )
}

expr_mat <- candidate_expr_log2[
  ,
  selected_col_idx,
  drop = FALSE
]

colnames(expr_mat) <- gene_map$EntrezGeneSymbol

if (anyDuplicated(colnames(expr_mat)) > 0) {
  stop(
    "Duplicated gene symbols remained after outcome-independent collapse.",
    call. = FALSE
  )
}

feature_qc <- tibble::tibble(
  EntrezGeneSymbol = colnames(expr_mat),
  AptName = gene_map$AptName,
  n_non_missing = vapply(
    expr_mat,
    function(x) sum(is.finite(safe_numeric(x))),
    numeric(1)
  ),
  missing_prop = vapply(
    expr_mat,
    function(x) mean(!is.finite(safe_numeric(x))),
    numeric(1)
  ),
  sd = vapply(
    expr_mat,
    safe_sd,
    numeric(1)
  ),
  mad = vapply(
    expr_mat,
    safe_mad,
    numeric(1)
  )
) %>%
  dplyr::mutate(
    keep_non_missing = n_non_missing >= MIN_NON_MISSING_PER_FEATURE,
    keep_variance = is.finite(sd) & sd > 0,
    keep_final = keep_non_missing & keep_variance
  )

if (REMOVE_ZERO_VARIANCE_FEATURES) {
  keep_genes <- feature_qc %>%
    dplyr::filter(keep_final) %>%
    dplyr::pull(EntrezGeneSymbol)

  if (
    REQUIRE_EXACT_PRIMARY_GENE_UNIVERSE &&
    length(keep_genes) != n_primary_genes
  ) {
    excluded_genes <- setdiff(
      colnames(expr_mat),
      keep_genes
    )

    safe_write_csv(
      tibble::tibble(
        EntrezGeneSymbol = excluded_genes
      ),
      file.path(OUTDIR, "qc", "genes_excluded_after_feature_qc.csv")
    )

    stop(
      "Feature QC would reduce the fixed primary gene universe from ",
      n_primary_genes,
      " to ",
      length(keep_genes),
      ". See qc/genes_excluded_after_feature_qc.csv.",
      call. = FALSE
    )
  }

  expr_mat <- expr_mat[
    ,
    keep_genes,
    drop = FALSE
  ]
}

missing_before_imputation <- sum(
  !is.finite(as.matrix(expr_mat))
)

for (j in seq_len(ncol(expr_mat))) {
  x <- safe_numeric(expr_mat[[j]])
  miss <- !is.finite(x)

  if (any(miss)) {
    med <- stats::median(x, na.rm = TRUE)

    if (!is.finite(med)) {
      stop(
        "Median imputation failed for gene: ",
        colnames(expr_mat)[[j]],
        call. = FALSE
      )
    }

    x[miss] <- med
  }

  expr_mat[[j]] <- x
}

missing_after_imputation <- sum(
  !is.finite(as.matrix(expr_mat))
)

if (missing_after_imputation != 0) {
  stop(
    "Non-finite values remained after median imputation.",
    call. = FALSE
  )
}

wgcna_expression <- dplyr::bind_cols(
  wgcna_metadata %>%
    dplyr::select(SampleId),
  tibble::as_tibble(expr_mat)
)

if (nrow(wgcna_expression) != nrow(wgcna_metadata)) {
  stop(
    "Expression and metadata sample counts are not aligned.",
    call. = FALSE
  )
}

if (!all(wgcna_expression$SampleId == wgcna_metadata$SampleId)) {
  stop(
    "Expression and metadata sample order is not aligned.",
    call. = FALSE
  )
}

# =============================================================================
# 11) BUILD GENE-LEVEL ANNOTATION
# =============================================================================

gene_level_annotation <- gene_map %>%
  dplyr::filter(
    EntrezGeneSymbol %in% colnames(expr_mat)
  ) %>%
  dplyr::select(
    EntrezGeneSymbol,
    Protein_Name,
    AptName,
    SeqId,
    dplyr::any_of(c(
      "TargetFullName",
      "Target",
      "UniProt",
      "EntrezGeneID",
      "Organism",
      "Type"
    )),
    missing_prop,
    MAD_log2_RFU,
    SD_log2_RFU,
    n_eligible_candidates,
    selection_rule,
    same_somamer_as_primary_DEP,
    logFC,
    P.Value,
    adj.P.Val,
    type,
    Direction,
    Primary_DEP_AptName,
    Primary_DEP_SeqId,
    Primary_DEP_Protein_Name,
    Primary_DEP_logFC,
    Primary_DEP_P.Value,
    Primary_DEP_adj.P.Val,
    Primary_DEP_type,
    Primary_DEP_Direction
  ) %>%
  dplyr::arrange(
    EntrezGeneSymbol
  )

if (nrow(gene_level_annotation) != ncol(expr_mat)) {
  stop(
    "Gene-level annotation does not match the exported expression feature count.",
    call. = FALSE
  )
}

# =============================================================================
# 12) QC SUMMARIES AND MANIFEST
# =============================================================================

sample_qc <- wgcna_metadata %>%
  dplyr::count(
    SampleGroup,
    Country,
    name = "n"
  ) %>%
  dplyr::arrange(
    Country,
    SampleGroup
  )

wgcna_input_qc_summary <- tibble::tibble(
  metric = c(
    "wgcna_project_root",
    "dep_project_root",
    "dep_workspace_file",
    "dep_gene_table",
    "wgcna_output_directory",
    "n_samples_wgcna",
    "n_CN",
    "n_AD",
    "n_countries",
    "n_primary_DEP_genes",
    "n_outcome_independent_gene_features",
    "n_unique_outcome_independent_SOMAmers",
    "n_same_SOMAmer_as_primary_DEP",
    "n_different_SOMAmer_from_primary_DEP",
    "proportion_same_SOMAmer_as_primary_DEP",
    "expression_scale_mode_requested",
    "expression_scale_resolved",
    "expression_transformation",
    "remove_zero_variance_features",
    "min_non_missing_per_feature",
    "missing_values_before_imputation",
    "missing_values_after_imputation",
    "wgcna_input_level",
    "collapse_rule",
    "gene_universe_rule"
  ),
  value = c(
    WGCNA_PROJECT_ROOT,
    DEP_PROJECT_ROOT,
    DEP_WORKSPACE_FILE,
    DEP_GENE_TABLE,
    OUTDIR,
    as.character(nrow(wgcna_metadata)),
    as.character(sum(wgcna_metadata$SampleGroup == "CN", na.rm = TRUE)),
    as.character(sum(wgcna_metadata$SampleGroup == "AD", na.rm = TRUE)),
    as.character(dplyr::n_distinct(wgcna_metadata$Country)),
    as.character(n_primary_genes),
    as.character(ncol(expr_mat)),
    as.character(dplyr::n_distinct(gene_map$AptName)),
    as.character(n_same_map),
    as.character(n_changed_map),
    as.character(n_same_map / nrow(gene_map)),
    EXPRESSION_SCALE_MODE,
    scale_transformation$resolved_mode,
    scale_transformation$transformation,
    as.character(REMOVE_ZERO_VARIANCE_FEATURES),
    as.character(MIN_NON_MISSING_PER_FEATURE),
    as.character(missing_before_imputation),
    as.character(missing_after_imputation),
    "GENE-COLLAPSED, not aptamer-level",
    paste(
      "One representative AptName per EntrezGeneSymbol:",
      "lowest missingness > highest log2 MAD > AptName tie-breaker"
    ),
    "Exact EntrezGeneSymbol universe from the primary gene-collapsed DEP table"
  )
)

wgcna_input_manifest <- tibble::tibble(
  output_file = c(
    "gene_collapsed_expression_matrix.csv",
    "wgcna_sample_metadata.csv",
    "gene_level_annotation.csv",
    "tables/outcome_independent_gene_somamer_map.csv",
    "tables/gene_map_comparison_vs_primary_DEP.csv",
    "tables/aptamer_selection_qc.csv",
    "tables/gene_level_feature_qc.csv",
    "tables/sample_counts_by_country_group.csv",
    "qc/expression_scale_audit.csv",
    "qc/gene_map_audit.csv",
    "qc/primary_genes_without_candidate_somamer.csv",
    "qc/primary_genes_without_eligible_candidate.csv",
    "wgcna_input_qc_summary.csv",
    "wgcna_input_manifest.csv",
    "sessionInfo.txt"
  ),
  description = c(
    "Sample-by-gene log2 RFU matrix for WGCNA using outcome-independent SOMAmer selection.",
    "Metadata aligned exactly to expression-matrix sample order.",
    "Gene-level annotation for the selected WGCNA SOMAmer, including its own aptamer-level DEP statistics and the canonical primary DEP representative for comparison.",
    "Final outcome-independent EntrezGeneSymbol-SOMAmer map.",
    "Gene-by-gene comparison between the WGCNA-selected SOMAmer and the primary DEP representative.",
    "Candidate SOMAmer missingness, log2 MAD, log2 SD and selection rank.",
    "Feature-level missingness and variance QC after gene collapse.",
    "Sample counts by country and diagnostic group.",
    "Automatic raw-RFU versus log2-RFU scale audit.",
    "Gene-map size, uniqueness and overlap audit.",
    "Primary genes without any eligible annotated candidate before QC.",
    "Primary genes without an eligible candidate after missingness and variance QC.",
    "Global QC and provenance summary.",
    "Manifest describing all exported files.",
    "R session information for reproducibility."
  )
)

# =============================================================================
# 13) EXPORT
# =============================================================================

safe_write_csv(
  wgcna_expression,
  file.path(
    OUTDIR,
    "gene_collapsed_expression_matrix.csv"
  )
)

safe_write_csv(
  wgcna_metadata,
  file.path(
    OUTDIR,
    "wgcna_sample_metadata.csv"
  )
)

safe_write_csv(
  gene_level_annotation,
  file.path(
    OUTDIR,
    "gene_level_annotation.csv"
  )
)

safe_write_csv(
  feature_qc,
  file.path(
    OUTDIR,
    "tables",
    "gene_level_feature_qc.csv"
  )
)

safe_write_csv(
  sample_qc,
  file.path(
    OUTDIR,
    "tables",
    "sample_counts_by_country_group.csv"
  )
)

safe_write_csv(
  wgcna_input_qc_summary,
  file.path(
    OUTDIR,
    "wgcna_input_qc_summary.csv"
  )
)

safe_write_csv(
  wgcna_input_manifest,
  file.path(
    OUTDIR,
    "wgcna_input_manifest.csv"
  )
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

message("WGCNA outcome-independent input builder complete.")
message("Output directory: ", OUTDIR)
message("Samples exported: ", nrow(wgcna_metadata))
message("Gene-level features exported: ", ncol(expr_mat))
message("Unique SOMAmers exported: ", dplyr::n_distinct(gene_map$AptName))
message("SOMAmers identical to primary DEP map: ", n_same_map)
message("SOMAmers different from primary DEP map: ", n_changed_map)
message("Missing values before imputation: ", missing_before_imputation)
message("Missing values after imputation: ", missing_after_imputation)
message("IMPORTANT: Script 11 must use gene_collapsed_expression_matrix.csv.")
message("IMPORTANT: WGCNA input level is GENE-COLLAPSED, not aptamer-level.")
message("IMPORTANT: SOMAmer selection did not use diagnosis, logFC, P value or FDR.")

###############################################################################
# END
###############################################################################

