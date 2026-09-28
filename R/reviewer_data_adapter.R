###############################################################################
# Processed-data compatibility layer
#
# Scientific principle:
#   The three reviewer CSVs are the only participant-level source inputs.
#   The proteomics CSV is already log2-transformed. Legacy R scripts were
#   originally written against normalized RFU from an ADAT and then applied
#   log2 internally. To preserve the submitted numerical workflow without
#   double-transforming, this adapter reconstructs normalized-RFU-compatible
#   values in memory as 2^(log2 RFU). No raw ADAT is required or created.
#
# Privacy principle:
#   Study_ID is a pseudonym. Some legacy scripts expect a column named SampleId.
#   We create SampleId <- Study_ID only as an in-memory compatibility alias.
#   It is NOT the original ReDLat SampleId.
###############################################################################

reviewer_read_metadata <- function(path) {
  if (!file.exists(path)) stop("Reviewer metadata file not found: ", path, call. = FALSE)
  x <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)

  # ML-derived metadata uses SampleId as the pseudonymous compatibility alias.
  # Accept it as Study_ID when Study_ID is not explicitly present.
  if (!"Study_ID" %in% names(x) && "SampleId" %in% names(x)) {
    x$Study_ID <- as.character(x$SampleId)
  }
  required <- c(
    "Study_ID", "SampleGroup", "Site", "Country", "Sex", "Age", "Education",
    "APOE4_carrier", "cdr_global", "cdr_boxscore", "mmse_total",
    "udsfaq_total", "NPI", "Mini.SEA", "T.ADLQ",
    "p.tau217", "p.tau181", "NfL", "ratio.AB42.40"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Reviewer metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(x$Study_ID)) stop("Duplicated Study_ID in reviewer metadata.", call. = FALSE)

  # Synthetic genotype labels encode carrier status only. They are compatibility
  # labels, not recovered APOE genotypes. All submitted APOE models use carrier
  # status, so this preserves the estimand without disclosing genotype detail.
  x <- x %>%
    dplyr::mutate(
      Study_ID = as.character(Study_ID),
      SampleId = Study_ID,
      SampleType = "Sample",
      RowCheck = NA_character_,
      PlateId = if ("AssayPlate" %in% names(.)) as.character(AssayPlate) else NA_character_,
      site = as.character(Site),
      SiteId = as.character(Site),
      ApoE = dplyr::case_when(
        APOE4_carrier == 1 ~ "e3/e4",
        APOE4_carrier == 0 ~ "e3/e3",
        TRUE ~ NA_character_
      ),
      APOE_group = dplyr::case_when(
        APOE4_carrier == 1 ~ "E4 carrier",
        APOE4_carrier == 0 ~ "Non-E4",
        TRUE ~ "Unknown"
      ),
      GFAP_1 = if ("GFAP_biomarker" %in% names(.)) suppressWarnings(as.numeric(GFAP_biomarker)) else NA_real_
    )
  x
}

reviewer_read_annotation <- function(path) {
  if (!file.exists(path)) stop("Reviewer annotation file not found: ", path, call. = FALSE)
  x <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  required <- c("AptName", "EntrezGeneSymbol", "Organism", "Type")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Reviewer annotation is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(x) != 10751L) warning("Expected 10,751 annotation rows; observed ", nrow(x))
  x
}

reviewer_read_proteomics_as_raw_compat <- function(path) {
  if (!file.exists(path)) stop("Reviewer proteomics file not found: ", path, call. = FALSE)
  x <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  if (!"Study_ID" %in% names(x)) stop("Proteomics file must contain Study_ID.", call. = FALSE)
  if (anyDuplicated(x$Study_ID)) stop("Duplicated Study_ID in reviewer proteomics.", call. = FALSE)
  seq_cols <- grep("^seq[._]", names(x), value = TRUE)
  if (length(seq_cols) != 10751L) {
    stop("Expected 10,751 SOMAmer columns; observed ", length(seq_cols), call. = FALSE)
  }

  # Invert the stored log2 transform in memory so legacy scripts can perform
  # their original single log2 transformation. This is a reversible numerical
  # interface adaptation only; no raw ADAT is reconstructed.
  mat_log2 <- as.matrix(x[, seq_cols, drop = FALSE])
  storage.mode(mat_log2) <- "double"
  mat_raw_compat <- 2^mat_log2
  mat_raw_compat[!is.finite(mat_log2)] <- NA_real_

  x <- data.frame(
    SampleId = as.character(x$Study_ID),
    mat_raw_compat,
    check.names = FALSE
  )
  x
}

reviewer_load_aptamer_data <- function(metadata_file, proteomics_file, annotation_file) {
  meta <- reviewer_read_metadata(metadata_file)
  prot_raw <- reviewer_read_proteomics_as_raw_compat(proteomics_file)
  annot <- reviewer_read_annotation(annotation_file)

  if (!all(meta$SampleId %in% prot_raw$SampleId)) {
    stop("Study_ID/SampleId sets differ between metadata and proteomics.", call. = FALSE)
  }

  sample_data <- meta %>%
    dplyr::left_join(prot_raw, by = "SampleId")

  seq_cols <- grep("^seq[._]", names(sample_data), value = TRUE)
  protein_universe <- annot %>%
    dplyr::filter(
      as.character(Organism) == "Human",
      as.character(Type) == "Protein",
      !is.na(EntrezGeneSymbol),
      EntrezGeneSymbol != "",
      AptName %in% seq_cols
    ) %>%
    dplyr::pull(AptName) %>%
    unique()

  if (length(protein_universe) != 10751L) {
    stop("Expected 10,751 eligible human protein SOMAmers; observed ",
         length(protein_universe), call. = FALSE)
  }

  list(
    metadata = meta,
    proteomics_raw_compat = prot_raw,
    sample_data = sample_data,
    annotation = annot,
    seq_cols = seq_cols,
    protein_universe = protein_universe
  )
}
