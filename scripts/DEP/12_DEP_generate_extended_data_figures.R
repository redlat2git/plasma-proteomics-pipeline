###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 12. DEP Extended Data figures
# Requires: Approved DEP Source Data workbooks
# Produces: Editable and submission-ready Extended Data figures
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

# The script is optimized for source(...), but it no longer aborts when opened
# interactively. This avoids the artificial startup error seen in v5.
running_through_source <- any(vapply(
  sys.calls(),
  function(call_i) {
    call_name <- tryCatch(
      as.character(call_i[[1]])[1],
      error = function(e) ""
    )
    call_name %in% c("source", "sys.source")
  },
  logical(1)
))

if (interactive() && !running_through_source) {
  message(
    "Interactive execution detected. The recommended command is: ",
    "source('scripts/DEP/12_DEP_generate_extended_data_figures.R', ",
    "encoding = 'UTF-8')."
  )
}

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "forcats", "tibble",
  "ggplot2", "ggrepel", "patchwork", "openxlsx", "scales",
  "viridis", "svglite"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(
  required_pkgs,
  function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
))

message("Extended Data Figure pipeline v6 initialized successfully.")


# -----------------------------------------------------------------------------
# 1. Project paths and submission settings
# -----------------------------------------------------------------------------


first_existing_dir <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[dir.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

source_root_env <- Sys.getenv("DEP_SOURCE_DATA_DIR", unset = "")
source_root <- first_existing_dir(c(
  file.path(publication_root, "source_data")
))
if (is.na(source_root)) {
  stop(
    "Could not locate the NatureAging_DEP Source Data directory.\n",
    "Generate Source Data first with scripts/DEP/09_DEP_generate_source_data.R.",
    call. = FALSE
  )
}

out_root <- file.path(publication_root, "figures", "extended_data")

submission_root <- file.path(out_root, "submission")
illustrator_root <- file.path(out_root, "Illustrator_ready")
master_vector_root <- file.path(illustrator_root, "master_vectors")
panel_vector_root <- file.path(illustrator_root, "panels_vector")
preview_root <- file.path(out_root, "preview")
panel_preview_root <- file.path(preview_root, "panels")
diagnostics_root <- file.path(out_root, "diagnostics")

for (d in c(
  submission_root,
  master_vector_root,
  panel_vector_root,
  preview_root,
  panel_preview_root,
  diagnostics_root
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ED_WIDTH_MM <- 180
ED1_HEIGHT_MM <- 220
ED2_HEIGHT_MM <- 200
ED3_HEIGHT_MM <- 220
ED4_HEIGHT_MM <- 195
ED_DPI <- 300

# Hard publication-dimension audit. Typography and label changes must be
# accommodated inside these artboards rather than enlarging the figures.
expected_dimensions <- c(
  ED_WIDTH_MM = 180,
  ED1_HEIGHT_MM = 220,
  ED2_HEIGHT_MM = 200,
  ED3_HEIGHT_MM = 220,
  ED4_HEIGHT_MM = 195
)
observed_dimensions <- c(
  ED_WIDTH_MM = ED_WIDTH_MM,
  ED1_HEIGHT_MM = ED1_HEIGHT_MM,
  ED2_HEIGHT_MM = ED2_HEIGHT_MM,
  ED3_HEIGHT_MM = ED3_HEIGHT_MM,
  ED4_HEIGHT_MM = ED4_HEIGHT_MM
)
if (!isTRUE(all.equal(observed_dimensions, expected_dimensions))) {
  stop(
    "Extended Data publication dimensions changed unexpectedly. ",
    "Do not enlarge the artboards to accommodate text.",
    call. = FALSE
  )
}

# Existing visual identity, adjusted for high contrast and colour-blind access.
COL_ORANGE <- "#F46D43"
COL_BLUE <- "#3288BD"
COL_DARK <- "#252525"
COL_MID <- "#8C8C8C"
COL_LIGHT <- "#D9D9D9"
COL_PALE <- "#F2F2F2"
COL_WHITE <- "#FFFFFF"
COL_PURPLE <- "#7B3294"
COL_GREEN <- "#008837"

FONT_ARIAL <- "Arial"
FONT_MYRIAD <- "Myriad Pro"

# Publication typography. Theme text sizes are specified directly in points.
# geom_text()/annotate() sizes are converted from points to millimetres with
# pt_to_geom(), so changing typography never changes the fixed artboard.
TEXT_PT <- list(
  panel_tag = 9,
  panel_title = 7,
  axis_title = 7,
  axis_tick = 6,
  protein = 6,
  legend = 6,
  legend_title = 6,
  heatmap_colourbar_tick = 5,
  compact_numeric = 5,
  micro = 4.08,
  heatmap_symbol = 6,
  exponent = 3.5
)

pt_to_geom <- function(pt) pt / ggplot2::.pt

BASE_SIZE <- TEXT_PT$axis_tick
TAG_SIZE <- TEXT_PT$panel_tag
PANEL_TITLE_SIZE <- TEXT_PT$panel_title
LINE_MM <- 0.30
POINT_SIZE <- 0.72

# Gene symbols are italicized in all displayed titles and axis labels.
APOE_PCA_TITLE <- expression(
  "PCA by " * italic(APOE) * " " * epsilon * 4
)
APOE_CARRIER_TITLE <- expression(
  italic(APOE) * " " * epsilon * 4 * " carrier status"
)
APOE_ADJUSTED_TITLE <- expression(
  italic(APOE) * " " * epsilon * 4 * "-adjusted model"
)
APOE_PRESERVATION_TITLE <- expression(
  "Effect preservation after " *
    italic(APOE) * " " * epsilon * 4 * " adjustment"
)
APOE_ADJUSTED_YLAB <- expression(
  italic(APOE) * "-adjusted " * log[2] * " fold-change"
)

# -----------------------------------------------------------------------------
# 2. Input helpers and strict checks
# -----------------------------------------------------------------------------

source_file <- function(filename) {
  path <- file.path(source_root, filename)
  if (!file.exists(path)) stop("Required Source Data workbook not found: ", path, call. = FALSE)
  path
}

read_source_sheet <- function(filename, sheet) {
  path <- source_file(filename)
  sheets <- openxlsx::getSheetNames(path)
  if (!sheet %in% sheets) {
    stop("Sheet '", sheet, "' not found in ", basename(path), call. = FALSE)
  }
  out <- openxlsx::read.xlsx(
    path,
    sheet = sheet,
    startRow = 4,
    colNames = TRUE,
    check.names = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = TRUE
  )
  tibble::as_tibble(out)
}

assert_cols <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_chr <- function(x) trimws(as.character(x))

normalize_site_key <- function(x) {
  z <- as_chr(x)
  z <- suppressWarnings(
    iconv(z, from = "", to = "ASCII//TRANSLIT")
  )
  z[is.na(z)] <- as_chr(x)[is.na(z)]
  z <- tolower(z)
  z <- gsub("[^a-z0-9]+", " ", z)
  trimws(z)
}

as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(as_chr(x)) %in% c("true", "t", "1", "yes", "y")
}
neglog10 <- function(x) -log10(pmax(as_num(x), .Machine$double.xmin))

recode_diagnosis <- function(x) {
  z <- toupper(as_chr(x))
  dplyr::case_when(
    z %in% c("CN", "CONTROL", "COGNITIVELY NORMAL") ~ "CN",
    z %in% c("AD", "ALZHEIMER", "ALZHEIMER'S DISEASE") ~ "AD",
    TRUE ~ as_chr(x)
  )
}

recode_sex <- function(x) {
  z <- tolower(as_chr(x))
  dplyr::case_when(
    z %in% c("1", "f", "female", "woman", "women") ~ "Female",
    z %in% c("2", "m", "male", "man", "men") ~ "Male",
    TRUE ~ as_chr(x)
  )
}

recode_apoe <- function(x) {
  z <- tolower(as_chr(x))
  dplyr::case_when(
    z %in% c("1", "true", "carrier", "yes") ~ "ε4 carrier",
    z %in% c("0", "false", "non-carrier", "noncarrier", "no") ~ "Non-carrier",
    TRUE ~ as_chr(x)
  )
}

format_p <- function(p) {
  p <- as_num(p)
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 1e-4 ~ format(p, scientific = TRUE, digits = 2),
    TRUE ~ sprintf("%.3f", p)
  )
}

significance_stars <- function(fdr) {
  fdr <- as_num(fdr)
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# -----------------------------------------------------------------------------
# 3. Nature-style theme and robust export functions
# -----------------------------------------------------------------------------

nature_theme <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = FONT_ARIAL
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = FONT_ARIAL,
        colour = COL_DARK
      ),
      axis.title = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$axis_title,
        face = "plain",
        colour = COL_DARK
      ),
      axis.text = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$axis_tick,
        face = "plain",
        colour = COL_DARK
      ),
      axis.line = ggplot2::element_line(
        linewidth = LINE_MM,
        colour = COL_DARK
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = LINE_MM,
        colour = COL_DARK
      ),
      axis.ticks.length = grid::unit(1.1, "mm"),
      plot.title = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$panel_title,
        face = "plain",
        hjust = 0,
        margin = ggplot2::margin(b = 1.8, unit = "pt")
      ),
      plot.subtitle = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$legend,
        face = "plain",
        hjust = 0
      ),
      legend.title = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$legend_title,
        face = "plain"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$legend,
        face = "plain"
      ),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(2.8, "mm"),
      legend.spacing.x = grid::unit(0.9, "mm"),
      legend.spacing.y = grid::unit(0.5, "mm"),
      legend.box.spacing = grid::unit(0.5, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      legend.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      # The small top margin reserves room for the elevated panel tag without
      # increasing any final figure dimension.
      plot.margin = ggplot2::margin(
        1.8, 1.6, 1.2, 1.6,
        unit = "mm"
      )
    )
}

tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    family = FONT_ARIAL,
    face = "bold",
    size = TEXT_PT$panel_tag,
    colour = COL_DARK,
    hjust = 0,
    vjust = 1
  ),
  # A subtle elevation above the plotting region. This is intentionally small
  # so letters remain visually attached to their panels.
  plot.tag.position = c(0, 1.006)
)

OUTPUT_REGISTRY <- tibble::tibble(
  display_item = character(),
  panel = character(),
  file = character(),
  path = character(),
  format = character(),
  width_mm = numeric(),
  height_mm = numeric(),
  intended_use = character()
)

register_output <- function(path, display_item, panel, width_mm, height_mm,
                            intended_use) {
  OUTPUT_REGISTRY <<- dplyr::bind_rows(
    OUTPUT_REGISTRY,
    tibble::tibble(
      display_item = display_item,
      panel = panel,
      file = basename(path),
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      format = toupper(tools::file_ext(path)),
      width_mm = width_mm,
      height_mm = height_mm,
      intended_use = intended_use
    )
  )
  invisible(path)
}

validate_graphics_file <- function(path, min_bytes = 5000L) {
  if (!file.exists(path)) {
    stop("Graphics file was not created: ", path, call. = FALSE)
  }
  size <- file.info(path)$size
  if (!is.finite(size) || size < min_bytes) {
    stop(
      "Graphics output is unexpectedly small or blank (",
      size, " bytes): ", path,
      call. = FALSE
    )
  }
  invisible(path)
}

validate_svg_file <- function(path, min_bytes = 1000L) {
  if (!file.exists(path)) {
    stop("SVG output was not created: ", path, call. = FALSE)
  }

  size <- file.info(path)$size
  if (!is.finite(size) || size < min_bytes) {
    stop(
      "SVG output is unexpectedly small or blank (",
      size, " bytes): ", path,
      call. = FALSE
    )
  }

  # Byte-level reading is stable on Windows and preserves UTF-8 labels.
  svg_raw <- readBin(
    con = path,
    what = "raw",
    n = as.integer(size)
  )
  svg_text <- rawToChar(svg_raw)
  svg_lower <- tolower(svg_text)

  has_opening_svg <- grepl(
    "<svg",
    svg_lower,
    fixed = TRUE
  )
  has_closing_svg <- grepl(
    "</svg>",
    svg_lower,
    fixed = TRUE
  )
  has_dimensions <- (
    grepl("viewbox=", svg_lower, fixed = TRUE) ||
      grepl("width=", svg_lower, fixed = TRUE) ||
      grepl("height=", svg_lower, fixed = TRUE)
  )

  graphics_tokens <- c(
    "<path", "<text", "<rect", "<circle", "<line",
    "<polyline", "<polygon", "<use", "<g"
  )

  has_graphics <- any(vapply(
    graphics_tokens,
    function(token) grepl(token, svg_lower, fixed = TRUE),
    logical(1)
  ))

  checks <- c(
    opening_svg = has_opening_svg,
    closing_svg = has_closing_svg,
    dimensions = has_dimensions,
    graphics = has_graphics
  )

  if (!all(checks)) {
    stop(
      "SVG failed structural validation: ",
      paste(names(checks)[!checks], collapse = ", "),
      ". File: ", path,
      call. = FALSE
    )
  }

  invisible(path)
}

atomic_replace <- function(tmp, final) {
  if (file.exists(final)) unlink(final)
  ok <- file.rename(tmp, final)
  if (!ok) {
    ok <- file.copy(tmp, final, overwrite = TRUE)
    unlink(tmp)
  }
  if (!ok) {
    stop(
      "Could not move temporary graphics file to: ",
      final,
      call. = FALSE
    )
  }
  invisible(final)
}

safe_pdf <- function(plot, filename, width_mm, height_mm) {
  tmp <- paste0(filename, ".tmp.pdf")
  if (file.exists(tmp)) unlink(tmp)

  grDevices::cairo_pdf(
    filename = tmp,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = FONT_ARIAL,
    bg = "white",
    onefile = TRUE
  )

  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop(
      "PDF drawing failed for ",
      basename(filename), ": ",
      conditionMessage(err),
      call. = FALSE
    )
  }

  validate_graphics_file(tmp, 5000L)
  atomic_replace(tmp, filename)
  validate_graphics_file(filename, 5000L)
  invisible(filename)
}

safe_svg <- function(plot, filename, width_mm, height_mm) {
  tmp <- paste0(filename, ".tmp.svg")
  if (file.exists(tmp)) unlink(tmp)

  svglite::svglite(
    file = tmp,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white"
  )

  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop(
      "SVG drawing failed for ",
      basename(filename), ": ",
      conditionMessage(err),
      call. = FALSE
    )
  }

  validate_svg_file(tmp, min_bytes = 1000L)
  atomic_replace(tmp, filename)
  validate_svg_file(filename, min_bytes = 1000L)
  invisible(filename)
}

safe_raster <- function(plot, filename, width_mm, height_mm,
                        dpi = 300L, type = c("png", "tiff")) {
  type <- match.arg(type)
  tmp <- paste0(filename, ".tmp.", type)
  if (file.exists(tmp)) unlink(tmp)

  width_px <- round(width_mm / 25.4 * dpi)
  height_px <- round(height_mm / 25.4 * dpi)

  if (requireNamespace("ragg", quietly = TRUE)) {
    if (type == "png") {
      ragg::agg_png(
        tmp,
        width = width_px,
        height = height_px,
        units = "px",
        res = dpi,
        background = "white"
      )
    } else {
      ragg::agg_tiff(
        tmp,
        width = width_px,
        height = height_px,
        units = "px",
        res = dpi,
        compression = "lzw",
        background = "white"
      )
    }
  } else if (type == "png") {
    grDevices::png(
      tmp,
      width = width_px,
      height = height_px,
      res = dpi,
      bg = "white",
      type = "cairo"
    )
  } else {
    grDevices::tiff(
      tmp,
      width = width_px,
      height = height_px,
      res = dpi,
      compression = "lzw",
      bg = "white",
      type = "cairo"
    )
  }

  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop(
      toupper(type), " drawing failed for ",
      basename(filename), ": ",
      conditionMessage(err),
      call. = FALSE
    )
  }

  validate_graphics_file(tmp, 10000L)
  atomic_replace(tmp, filename)
  validate_graphics_file(filename, 10000L)
  invisible(filename)
}

save_extended_figure <- function(plot, display_item, stem,
                                 width_mm, height_mm) {
  tiff_path <- file.path(submission_root, paste0(stem, ".tiff"))
  pdf_path <- file.path(master_vector_root, paste0(stem, ".pdf"))
  svg_path <- file.path(master_vector_root, paste0(stem, ".svg"))
  png_path <- file.path(preview_root, paste0(stem, "_preview.png"))

  safe_raster(
    plot, tiff_path, width_mm, height_mm,
    dpi = ED_DPI, type = "tiff"
  )
  safe_pdf(plot, pdf_path, width_mm, height_mm)
  safe_svg(plot, svg_path, width_mm, height_mm)
  safe_raster(
    plot, png_path, width_mm, height_mm,
    dpi = 300L, type = "png"
  )

  register_output(
    tiff_path, display_item, "complete figure",
    width_mm, height_mm,
    "Nature Aging submission TIFF"
  )
  register_output(
    pdf_path, display_item, "complete figure",
    width_mm, height_mm,
    "Editable vector PDF"
  )
  register_output(
    svg_path, display_item, "complete figure",
    width_mm, height_mm,
    "Preferred Adobe Illustrator SVG"
  )
  register_output(
    png_path, display_item, "complete figure",
    width_mm, height_mm,
    "Review preview"
  )

  message("Generated ", display_item, ": ", tiff_path)
  invisible(c(tiff_path, pdf_path, svg_path, png_path))
}

save_panel <- function(plot, display_item, panel, stem,
                       width_mm, height_mm) {
  pdf_path <- file.path(panel_vector_root, paste0(stem, ".pdf"))
  svg_path <- file.path(panel_vector_root, paste0(stem, ".svg"))
  png_path <- file.path(panel_preview_root, paste0(stem, "_preview.png"))

  safe_pdf(plot, pdf_path, width_mm, height_mm)
  safe_svg(plot, svg_path, width_mm, height_mm)
  safe_raster(
    plot, png_path, width_mm, height_mm,
    dpi = 300L, type = "png"
  )

  register_output(
    pdf_path, display_item, panel,
    width_mm, height_mm,
    "Editable panel PDF"
  )
  register_output(
    svg_path, display_item, panel,
    width_mm, height_mm,
    "Editable panel SVG"
  )
  register_output(
    png_path, display_item, panel,
    width_mm, height_mm,
    "Panel review preview"
  )

  message(
    "Panel ", display_item, panel,
    " SVG created (",
    file.info(svg_path)$size,
    " bytes)"
  )

  invisible(c(pdf_path, svg_path, png_path))
}

# -----------------------------------------------------------------------------
# 4. SVG export preflight
# -----------------------------------------------------------------------------

svg_preflight_path <- file.path(
  diagnostics_root,
  "SVG_export_preflight.svg"
)

svg_preflight_plot <- ggplot2::ggplot(
  data.frame(x = 1, y = 1),
  ggplot2::aes(x, y)
) +
  ggplot2::geom_point(size = 2) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = 1.12,
    label = "SVG preflight",
    family = FONT_ARIAL,
    size = 2
  ) +
  ggplot2::theme_void(base_family = FONT_ARIAL) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    )
  )

safe_svg(
  svg_preflight_plot,
  svg_preflight_path,
  width_mm = 35,
  height_mm = 25
)

message(
  "SVG structural preflight passed (",
  file.info(svg_preflight_path)$size,
  " bytes)."
)

# -----------------------------------------------------------------------------
# 5. Load approved Extended Data Source Data
# -----------------------------------------------------------------------------

pca_file <- "SourceData_DEP_PCA.xlsx"
demog_file <- "SourceData_DEP_DemographicDistributions.xlsx"
pvalue_file <- "SourceData_DEP_PvalueCalibration.xlsx"
protein_file <- "SourceData_DEP_RepresentativeProteins.xlsx"
sensitivity_file <- "SourceData_DEP_InternalSensitivity.xlsx"
robust_file <- "SourceData_DEP_RecruitmentRobustness.xlsx"
ptau_file <- "SourceData_DEP_pTau217ThresholdSweep.xlsx"

pca_scores <- read_source_sheet(pca_file, "PCA_scores")
pca_variance <- read_source_sheet(pca_file, "PCA_variance")
demog <- read_source_sheet(demog_file, "Participant_data")
pvalues <- read_source_sheet(pvalue_file, "Pvalue_distribution")
qq_values <- read_source_sheet(pvalue_file, "QQ_plot")

boxplot_data <- read_source_sheet(
  protein_file,
  "Boxplot_participant_data"
)
heatmap_selection <- read_source_sheet(
  protein_file,
  "Heatmap_protein_selection"
)
heatmap_corr <- read_source_sheet(
  protein_file,
  "Heatmap_correlations"
)

apoe_volcano <- read_source_sheet(sensitivity_file, "APOE_volcano")
apoe_compare <- read_source_sheet(sensitivity_file, "APOE_comparison")
atn_volcano <- read_source_sheet(sensitivity_file, "ATN_volcano")
atn_compare <- read_source_sheet(sensitivity_file, "ATN_comparison")
cdrsb_volcano <- read_source_sheet(sensitivity_file, "CDRSB_volcano")
cdrsb_compare <- read_source_sheet(sensitivity_file, "CDRSB_comparison")

ptau_summary <- read_source_sheet(ptau_file, "Threshold_summary")
ptau_counts <- read_source_sheet(ptau_file, "Sample_counts")
ptau_p90 <- read_source_sheet(ptau_file, "P90_comparison")

loco <- read_source_sheet(robust_file, "LOCO_summary")
loso <- read_source_sheet(robust_file, "LOSO_summary")
balanced <- read_source_sheet(robust_file, "Balanced_iterations")
country_meta <- read_source_sheet(robust_file, "Country_meta")

if (!"PlateId" %in% names(pca_scores)) pca_scores$PlateId <- NA_character_
if (!"site" %in% names(pca_scores)) pca_scores$site <- if ("Site" %in% names(pca_scores)) pca_scores$Site else if ("SiteId" %in% names(pca_scores)) pca_scores$SiteId else NA_character_
if (!"SiteId" %in% names(pca_scores)) pca_scores$SiteId <- pca_scores$site

assert_cols(
  pca_scores,
  c("PC1", "PC2", "Age", "SampleGroup", "Sex",
    "APOE4_carrier", "PlateId", "Country"),
  "PCA scores"
)
assert_cols(pca_variance, c("PC", "percent"), "PCA variance")
assert_cols(
  demog,
  c("SampleGroup", "Age", "Sex", "APOE4_carrier"),
  "Demographic participant data"
)
assert_cols(pvalues, "P.Value", "P-value distribution")
assert_cols(
  qq_values,
  c("expected_minus_log10_P", "observed_minus_log10_P"),
  "QQ data"
)
assert_cols(
  boxplot_data,
  c("Protein", "SampleGroup", "log2_RFU", "Direction"),
  "Representative-protein participant data"
)
assert_cols(
  heatmap_selection,
  c("EntrezGeneSymbol", "Direction"),
  "Extended heatmap selection"
)
assert_cols(
  heatmap_corr,
  c("Protein", "Trait", "rho", "adj.P.Val"),
  "Extended heatmap correlations"
)
assert_cols(
  apoe_volcano,
  c("EntrezGeneSymbol", "logFC", "adj.P.Val"),
  "APOE volcano"
)
assert_cols(
  apoe_compare,
  c("EntrezGeneSymbol", "logFC_primary", "logFC_apoe",
    "adj.P.Val_primary", "adj.P.Val_apoe"),
  "APOE comparison"
)
assert_cols(
  atn_volcano,
  c("EntrezGeneSymbol", "logFC", "adj.P.Val"),
  "AT(N) volcano"
)
assert_cols(
  atn_compare,
  c("EntrezGeneSymbol", "logFC_primary", "logFC_atn",
    "adj.P.Val_primary", "adj.P.Val_atn"),
  "AT(N) comparison"
)
assert_cols(
  cdrsb_volcano,
  c("EntrezGeneSymbol", "logFC", "adj.P.Val"),
  "CDR-SB volcano"
)
assert_cols(
  cdrsb_compare,
  c("EntrezGeneSymbol", "logFC_primary", "beta_CDRSB_AD_only",
    "adj.P.Val_primary", "adj.P.Val_CDRSB_AD_only"),
  "CDR-SB comparison"
)
assert_cols(
  ptau_summary,
  c(
    "threshold_family",
    "p_tau217_CN_quantile",
    "fdr005_total",
    "n_CN_initial",
    "n_AD_initial",
    "logFC_correlation_with_primary",
    "direction_consistency_all",
    "preservation_fraction_primary_fdr005"
  ),
  "p-tau217 threshold summary"
)
assert_cols(
  ptau_counts,
  c("threshold_family", "p_tau217_CN_quantile",
    "BioContrast", "n_initial"),
  "p-tau217 sample counts"
)
assert_cols(
  ptau_p90,
  c("EntrezGeneSymbol", "logFC_primary", "logFC_secondary",
    "adj.P.Val_primary", "adj.P.Val_secondary"),
  "p-tau217 p90 comparison"
)
assert_cols(
  loco,
  c("excluded_country", "logFC_correlation",
    "prop_main_sig_preserved"),
  "LOCO summary"
)
assert_cols(
  loso,
  c("excluded_site", "logFC_correlation",
    "prop_main_sig_preserved"),
  "LOSO summary"
)
assert_cols(
  balanced,
  c("logFC_correlation", "direction_consistency_all",
    "prop_main_sig_preserved"),
  "Balanced-resampling iterations"
)
assert_cols(
  country_meta,
  c("EntrezGeneSymbol", "meta_logFC", "meta_se",
    "meta_adj.P.Val", "I2"),
  "Country meta-analysis"
)

# -----------------------------------------------------------------------------
# 6. Reusable panel builders
# -----------------------------------------------------------------------------

volcano_panel <- function(data, effect_col, fdr_col, gene_col, title,
                          labels = character(), count_label = TRUE) {
  d <- data %>%
    mutate(
      effect = as_num(.data[[effect_col]]),
      fdr = as_num(.data[[fdr_col]]),
      gene = as_chr(.data[[gene_col]]),
      y = neglog10(fdr),
      class = case_when(
        fdr < 0.05 & effect > 0 ~ "Higher in AD",
        fdr < 0.05 & effect < 0 ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    )

  label_data <- d %>% filter(gene %in% labels)
  n_high <- sum(d$class == "Higher in AD", na.rm = TRUE)
  n_low <- sum(d$class == "Lower in AD", na.rm = TRUE)

  p <- ggplot(d, aes(effect, y)) +
    geom_point(
      data = d %>% filter(class == "Not significant"),
      aes(fill = class), shape = 21, size = 0.62, stroke = 0.18,
      colour = COL_DARK, alpha = 0.50
    ) +
    geom_point(
      data = d %>% filter(class != "Not significant"),
      aes(fill = class), shape = 21, size = 0.72, stroke = 0.20,
      colour = COL_DARK, alpha = 0.88
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = LINE_MM, colour = COL_MID) +
    geom_vline(xintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = gene),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), colour = COL_DARK,
      min.segment.length = 0, segment.size = 0.20,
      box.padding = 0.18, point.padding = 0.10,
      max.overlaps = Inf, seed = 1408
    ) +
    scale_fill_manual(
      values = c(
        "Higher in AD" = COL_ORANGE,
        "Lower in AD" = COL_BLUE,
        "Not significant" = COL_LIGHT
      ),
      breaks = c("Higher in AD", "Lower in AD")
    ) +
    labs(
      title = title,
      x = expression(Adjusted~log[2]~fold-change~("AD versus CN")),
      y = expression(-log[10]~("BH-FDR")),
      fill = NULL
    ) +
    nature_theme() +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.2, alpha = 1)
      )
    )

  if (count_label) {
    xrange <- range(d$effect, finite = TRUE)
    yrange <- range(d$y, finite = TRUE)
    p <- p +
      annotate(
        "text",
        x = xrange[1] + 0.05 * diff(xrange),
        y = yrange[2] * 0.93,
        label = paste0(n_low, " lower"),
        size = pt_to_geom(TEXT_PT$protein),
        family = FONT_ARIAL,
        colour = COL_BLUE,
        fontface = "bold",
        hjust = 0
      ) +
      annotate(
        "text",
        x = xrange[2] - 0.05 * diff(xrange),
        y = yrange[2] * 0.93,
        label = paste0(n_high, " higher"),
        size = pt_to_geom(TEXT_PT$protein),
        family = FONT_ARIAL,
        colour = COL_ORANGE,
        fontface = "bold",
        hjust = 1
      )
  }
  p
}

effect_scatter_panel <- function(data, x_col, y_col, gene_col,
                                 primary_fdr_col, title,
                                 labels = character(),
                                 secondary_fdr_col = NULL,
                                 x_lab = "Primary-model log2 fold-change",
                                 y_lab = "Secondary-model effect",
                                 show_legend = FALSE) {
  d <- data %>%
    mutate(
      x = as_num(.data[[x_col]]),
      y = as_num(.data[[y_col]]),
      gene = as_chr(.data[[gene_col]]),
      primary_fdr = as_num(.data[[primary_fdr_col]]),
      class = case_when(
        primary_fdr < 0.05 & x > 0 ~ "Higher in AD",
        primary_fdr < 0.05 & x < 0 ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    ) %>%
    filter(is.finite(x), is.finite(y))

  r <- suppressWarnings(stats::cor(d$x, d$y, method = "pearson", use = "complete.obs"))
  direction <- mean(sign(d$x) == sign(d$y), na.rm = TRUE)
  label_data <- d %>% filter(gene %in% labels)

  lim <- max(abs(c(d$x, d$y)), na.rm = TRUE) * 1.05
  if (!is.finite(lim) || lim == 0) lim <- 1

  p <- ggplot(d, aes(x, y)) +
    geom_point(
      aes(fill = class), shape = 21, size = 0.65, stroke = 0.18,
      colour = COL_DARK, alpha = 0.68
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = LINE_MM, colour = COL_MID) +
    geom_hline(yintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
    geom_vline(xintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
    ggrepel::geom_text_repel(
      data = label_data, aes(label = gene),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), colour = COL_DARK,
      min.segment.length = 0, segment.size = 0.18,
      box.padding = 0.15, point.padding = 0.08,
      max.overlaps = Inf, seed = 2205
    ) +
    annotate(
      "text", x = lim * 0.95, y = -lim * 0.92,
      label = sprintf("r = %.2f\nDirection = %.1f%%", r, 100 * direction),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), hjust = 1, vjust = 0, colour = COL_DARK
    ) +
    scale_fill_manual(
      values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE, "Not significant" = COL_LIGHT),
      breaks = c("Higher in AD", "Lower in AD")
    ) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = title, x = x_lab, y = y_lab, fill = NULL) +
    nature_theme() +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "left"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.2, alpha = 1)
      )
    )

  if (!is.null(secondary_fdr_col) && secondary_fdr_col %in% names(data)) {
    n_secondary <- sum(as_num(data[[secondary_fdr_col]]) < 0.05, na.rm = TRUE)
    p <- p + annotate(
      "text", x = -lim * 0.96, y = lim * 0.94,
      label = paste0("Secondary FDR < 0.05: ", n_secondary),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), hjust = 0, vjust = 1
    )
  }
  p
}

robust_lollipop <- function(data, label_col, value_col, title,
                              percent = FALSE,
                              reference = NULL) {
  d <- data %>%
    transmute(
      label = as_chr(.data[[label_col]]),
      value = as_num(.data[[value_col]]),
      highlight = stringr::str_detect(tolower(label), "colombia")
    ) %>%
    filter(is.finite(value)) %>%
    arrange(value) %>%
    mutate(label = factor(label, levels = label))

  formatter <- if (percent) scales::label_percent(accuracy = 1) else scales::label_number(accuracy = 0.01)
  upper <- max(d$value, na.rm = TRUE)
  lower <- min(0, min(d$value, na.rm = TRUE))
  d <- d %>%
    mutate(
      label_text = formatter(value),
      label_x = value + 0.02 * max(upper - lower, 0.1)
    )

  p <- ggplot(d, aes(value, label)) +
    geom_segment(
      aes(x = lower, xend = value, yend = label),
      linewidth = 0.45,
      colour = COL_LIGHT
    ) +
    geom_point(aes(fill = highlight), shape = 21, size = 1.55, stroke = 0.25, colour = COL_DARK) +
    geom_text(
      aes(x = label_x, label = label_text),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), hjust = 0, colour = COL_DARK
    ) +
    scale_fill_manual(values = c(`TRUE` = COL_ORANGE, `FALSE` = COL_BLUE), guide = "none") +
    scale_x_continuous(labels = formatter, expand = expansion(mult = c(0.02, 0.16))) +
    labs(
      title = title,
      x = if (percent) {
        "Preserved primary signals"
      } else {
        "Correlation with primary model"
      },
      y = NULL
    ) +
    nature_theme() +
    theme(axis.text.y = element_text(size = TEXT_PT$protein))

  if (!is.null(reference) && is.finite(reference)) {
    p <- p + geom_vline(
      xintercept = reference,
      linetype = "dashed",
      linewidth = 0.30,
      colour = COL_MID
    )
  }

  p
}

balanced_hist <- function(data, value_col, title, x_lab, percent = FALSE) {
  values <- as_num(data[[value_col]])
  median_value <- median(values, na.rm = TRUE)
  formatter <- if (percent) scales::label_percent(accuracy = 1) else scales::label_number(accuracy = 0.01)

  ggplot(tibble(value = values), aes(value)) +
    geom_histogram(bins = 25, fill = COL_LIGHT, colour = COL_DARK, linewidth = 0.22) +
    geom_vline(xintercept = median_value, linetype = "dashed", linewidth = 0.45, colour = COL_ORANGE) +
    annotate(
      "text", x = median_value, y = Inf, label = paste0("Median ", formatter(median_value)),
      family = FONT_ARIAL, size = pt_to_geom(TEXT_PT$protein), hjust = -0.05, vjust = 1.3
    ) +
    scale_x_continuous(labels = formatter) +
    labs(title = title, x = x_lab, y = "Iterations") +
    nature_theme()
}

# -----------------------------------------------------------------------------
# 7. Extended Data Figure 1: QC, cohort structure and calibration
# -----------------------------------------------------------------------------

# Source Data versions may expose recruitment site as `site`, `SiteId`, or both.
# Create both placeholders before mutation so the script remains schema-safe.
if (!"site" %in% names(pca_scores)) {
  pca_scores$site <- NA_character_
}
if (!"SiteId" %in% names(pca_scores)) {
  pca_scores$SiteId <- NA_character_
}

pca_scores <- pca_scores %>%
  dplyr::mutate(
    PC1 = as_num(PC1),
    PC2 = as_num(PC2),
    Age = as_num(Age),
    SampleGroup = factor(
      recode_diagnosis(SampleGroup),
      levels = c("CN", "AD")
    ),
    Sex_label = factor(recode_sex(Sex)),
    APOE_label = factor(
      recode_apoe(APOE4_carrier),
      levels = c("Non-carrier", "ε4 carrier")
    ),
    # Remove only the requested prefix from the displayed plate labels.
    Plate_text = stringr::str_remove(
      as_chr(PlateId),
      "^CHI-24-002_"
    ),
    Country_text = as_chr(Country),
    Site_raw = dplyr::if_else(
      !is.na(site) & nzchar(as_chr(site)),
      as_chr(site),
      as_chr(SiteId)
    )
  )

# Never display investigator/site names. Countries with one site retain only
# the country label; countries with multiple sites become Country 1, Country 2,
# etc. This yields labels such as Chile 1/2 and Colombia 1/2 from the data.
site_anonymization <- pca_scores %>%
  dplyr::transmute(
    Country_text = as_chr(Country_text),
    Site_raw = as_chr(Site_raw)
  ) %>%
  dplyr::filter(
    !is.na(Country_text),
    nzchar(Country_text),
    !is.na(Site_raw),
    nzchar(Site_raw)
  ) %>%
  dplyr::distinct(Country_text, Site_raw) %>%
  dplyr::arrange(Country_text, Site_raw) %>%
  dplyr::group_by(Country_text) %>%
  dplyr::mutate(
    n_sites_country = dplyr::n(),
    site_number = dplyr::row_number(),
    Site_public = dplyr::if_else(
      n_sites_country > 1,
      paste0(Country_text, " ", site_number),
      Country_text
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Site_key = normalize_site_key(Site_raw)
  )

pca_scores <- pca_scores %>%
  dplyr::left_join(
    site_anonymization %>%
      dplyr::select(
        Country_text,
        Site_raw,
        Site_public
      ),
    by = c("Country_text", "Site_raw")
  ) %>%
  dplyr::mutate(
    Plate_label = factor(Plate_text),
    Country_label = factor(Country_text),
    Site_label = factor(Site_public)
  )

pc1_pct <- pca_variance %>% filter(as_chr(PC) == "PC1") %>% pull(percent) %>% as_num()
pc2_pct <- pca_variance %>% filter(as_chr(PC) == "PC2") %>% pull(percent) %>% as_num()
if (length(pc1_pct) == 0) pc1_pct <- NA_real_
if (length(pc2_pct) == 0) pc2_pct <- NA_real_

pca_axis_x <- if (is.finite(pc1_pct[1])) sprintf("PC1 (%.1f%%)", pc1_pct[1]) else "PC1"
pca_axis_y <- if (is.finite(pc2_pct[1])) sprintf("PC2 (%.1f%%)", pc2_pct[1]) else "PC2"

pca_discrete <- function(var, title, palette = NULL) {
  n_levels <- dplyr::n_distinct(pca_scores[[var]], na.rm = TRUE)
  legend_cols <- max(1L, min(3L, ceiling(n_levels / 2)))

  p <- ggplot(pca_scores, aes(PC1, PC2, colour = .data[[var]])) +
    geom_point(size = 0.66, alpha = 0.72) +
    labs(
      title = title,
      x = pca_axis_x,
      y = pca_axis_y,
      colour = NULL
    ) +
    nature_theme(5.8) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left"
    ) +
    guides(
      colour = guide_legend(
        ncol = legend_cols,
        byrow = TRUE,
        override.aes = list(size = 2.1, alpha = 1)
      )
    )

  if (!is.null(palette)) {
    p + scale_colour_manual(
      values = palette,
      na.value = COL_LIGHT
    )
  } else {
    p + viridis::scale_colour_viridis(
      discrete = TRUE,
      option = "D",
      end = 0.90,
      na.value = COL_LIGHT
    )
  }
}

pca_continuous <- function(var, title) {
  ggplot(pca_scores, aes(PC1, PC2, colour = .data[[var]])) +
    geom_point(size = 0.66, alpha = 0.75) +
    viridis::scale_colour_viridis(
      option = "C",
      end = 0.92,
      na.value = COL_LIGHT,
      guide = guide_colourbar(
        title.position = "top",
        barwidth = grid::unit(18, "mm"),
        barheight = grid::unit(2.4, "mm")
      )
    ) +
    labs(
      title = title,
      x = pca_axis_x,
      y = pca_axis_y,
      colour = "Years"
    ) +
    nature_theme(5.8) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left"
    )
}

ed1a <- pca_discrete("Sex_label", "PCA by sex", c("Female" = "#CC79A7", "Male" = COL_BLUE))
ed1b <- pca_continuous("Age", "PCA by age")
ed1c <- pca_discrete(
  "APOE_label",
  APOE_PCA_TITLE,
  c(
    "Non-carrier" = COL_MID,
    "ε4 carrier" = COL_ORANGE
  )
)
ed1d <- pca_discrete("Plate_label", "PCA by plate")
ed1e <- pca_discrete("Country_label", "PCA by country")
ed1f <- pca_discrete("Site_label", "PCA by recruitment site")
ed1g <- pca_discrete(
  "SampleGroup",
  "PCA by diagnosis",
  c("CN" = COL_BLUE, "AD" = COL_ORANGE)
)

ed1_gh_audit <- tibble::tibble(
  check = c(
    "Panel g contains both diagnostic groups",
    "Panel h contains positive finite ages"
  ),
  passed = c(
    all(c("CN", "AD") %in% unique(as_chr(pca_scores$SampleGroup))),
    any(is.finite(as_num(demog$Age)) & as_num(demog$Age) > 0)
  )
)

if (!all(ed1_gh_audit$passed)) {
  stop(
    "Extended Data Fig. 1 panels g/h failed their data audit: ",
    paste(ed1_gh_audit$check[!ed1_gh_audit$passed], collapse = "; "),
    call. = FALSE
  )
}

ed1h <- demog %>%
  transmute(
    Diagnosis = factor(recode_diagnosis(SampleGroup), levels = c("CN", "AD")),
    Age = as_num(Age)
  ) %>%
  ggplot(aes(Diagnosis, Age, fill = Diagnosis)) +
  geom_violin(width = 0.82, trim = FALSE, alpha = 0.40, colour = COL_DARK, linewidth = 0.25) +
  geom_boxplot(width = 0.23, outlier.shape = NA, colour = COL_DARK, linewidth = 0.30) +
  scale_fill_manual(values = c("CN" = COL_BLUE, "AD" = COL_ORANGE), guide = "none") +
  labs(title = "Age distribution", x = NULL, y = "Age (years)") +
  nature_theme(5.8)

sex_prop <- demog %>%
  transmute(
    Diagnosis = factor(recode_diagnosis(SampleGroup), levels = c("CN", "AD")),
    Sex = factor(recode_sex(Sex))
  ) %>%
  count(Diagnosis, Sex, name = "n") %>%
  group_by(Diagnosis) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ed1i <- ggplot(sex_prop, aes(Diagnosis, prop, fill = Sex)) +
  geom_col(width = 0.65, colour = COL_DARK, linewidth = 0.25) +
  scale_fill_manual(values = c("Female" = "#CC79A7", "Male" = COL_BLUE)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Sex distribution", x = NULL, y = "Participants", fill = NULL) +
  nature_theme(5.8) +
  theme(
    legend.position = "bottom",
    legend.justification = "left"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  )

apoe_prop <- demog %>%
  transmute(
    Diagnosis = factor(recode_diagnosis(SampleGroup), levels = c("CN", "AD")),
    APOE = factor(recode_apoe(APOE4_carrier), levels = c("Non-carrier", "ε4 carrier"))
  ) %>%
  filter(!is.na(APOE), APOE %in% c("Non-carrier", "ε4 carrier")) %>%
  count(Diagnosis, APOE, name = "n") %>%
  group_by(Diagnosis) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ed1j <- ggplot(apoe_prop, aes(Diagnosis, prop, fill = APOE)) +
  geom_col(width = 0.65, colour = COL_DARK, linewidth = 0.25) +
  scale_fill_manual(values = c("Non-carrier" = COL_MID, "ε4 carrier" = COL_ORANGE)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = APOE_CARRIER_TITLE,
    x = NULL,
    y = "Genotyped participants",
    fill = NULL
  ) +
  nature_theme(5.8) +
  theme(
    legend.position = "bottom",
    legend.justification = "left"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  )

ed1k <- pvalues %>%
  transmute(P = as_num(P.Value)) %>%
  filter(is.finite(P), P >= 0, P <= 1) %>%
  ggplot(aes(P)) +
  geom_histogram(bins = 40, fill = COL_MID, colour = COL_DARK, linewidth = 0.20) +
  labs(title = "Association P-value distribution", x = "Raw P value", y = "Proteins") +
  nature_theme(5.8)

ed1l <- qq_values %>%
  transmute(
    expected = as_num(expected_minus_log10_P),
    observed = as_num(observed_minus_log10_P)
  ) %>%
  filter(is.finite(expected), is.finite(observed)) %>%
  ggplot(aes(expected, observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = COL_ORANGE, linewidth = 0.40) +
  geom_point(size = 0.55, alpha = 0.72, colour = COL_DARK) +
  labs(
    title = "Quantile–quantile calibration",
    x = expression(Expected~-log[10](P)),
    y = expression(Observed~-log[10](P))
  ) +
  nature_theme(5.8)

ed1 <- patchwork::wrap_plots(
  ed1a, ed1b, ed1c, ed1d,
  ed1e, ed1f, ed1g, ed1h,
  ed1i, ed1j, ed1k, ed1l,
  ncol = 4,
  widths = rep(1, 4),
  heights = rep(1, 3),
  guides = "keep"
) +
  patchwork::plot_annotation(tag_levels = "a") &
  tag_theme

ed1_panel_list <- list(
  a = ed1a, b = ed1b, c = ed1c, d = ed1d,
  e = ed1e, f = ed1f, g = ed1g, h = ed1h,
  i = ed1i, j = ed1j, k = ed1k, l = ed1l
)

ed1_files <- save_extended_figure(
  ed1,
  "Extended Data Fig. 1",
  "ExtendedDataFig1_DEP_QC_CohortCalibration",
  ED_WIDTH_MM,
  ED1_HEIGHT_MM
)

purrr::iwalk(
  ed1_panel_list,
  ~ save_panel(
    .x,
    "Extended Data Fig. 1",
    .y,
    paste0("EDFig1", .y, "_DEP_QC_CohortCalibration"),
    45, 70
  )
)

# -----------------------------------------------------------------------------
# 8. Extended Data Figure 2: representative proteins and extended heatmap
# -----------------------------------------------------------------------------

protein_order_box <- unique(as_chr(boxplot_data$Protein))
boxplot_data2 <- boxplot_data %>%
  mutate(
    Diagnosis = factor(recode_diagnosis(SampleGroup), levels = c("CN", "AD")),
    Protein = factor(as_chr(Protein), levels = protein_order_box),
    log2_RFU = as_num(log2_RFU),
    Direction = as_chr(Direction),
    facet_label = paste0(as_chr(Protein), "\n", Direction)
  )
box_facet_order <- unique(boxplot_data2$facet_label)
boxplot_data2$facet_label <- factor(boxplot_data2$facet_label, levels = box_facet_order)

ed2a <- ggplot(boxplot_data2, aes(Diagnosis, log2_RFU, fill = Diagnosis)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.25, alpha = 0.60) +
  geom_jitter(width = 0.12, height = 0, size = 0.22, alpha = 0.22, colour = COL_DARK) +
  facet_wrap(~facet_label, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c("CN" = COL_BLUE, "AD" = COL_ORANGE), guide = "none") +
  labs(title = "Representative differentially abundant proteins", x = NULL, y = expression(Non-residualized~log[2]~RFU)) +
  nature_theme(5.2) +
  theme(
    strip.background = element_rect(fill = "white", colour = COL_DARK, linewidth = 0.22),
    strip.text = element_text(size = TEXT_PT$protein, family = FONT_ARIAL),
    axis.text.x = element_text(size = TEXT_PT$protein),
    axis.text.y = element_text(size = TEXT_PT$protein),
    panel.spacing = grid::unit(1.0, "mm")
  )

heat_protein_order <- unique(as_chr(heatmap_selection$EntrezGeneSymbol))
heat_trait_order <- unique(as_chr(heatmap_corr$Trait))

heat2 <- heatmap_corr %>%
  transmute(
    Protein = factor(as_chr(Protein), levels = heat_protein_order),
    Trait = factor(as_chr(Trait), levels = heat_trait_order),
    rho = as_num(rho),
    star = significance_stars(adj.P.Val)
  )

strip2 <- heatmap_selection %>%
  transmute(
    Protein = factor(as_chr(EntrezGeneSymbol), levels = heat_protein_order),
    Direction = as_chr(Direction),
    y = "DEP direction"
  )

ed2b_strip <- ggplot(strip2, aes(Protein, y, fill = Direction)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE), name = NULL) +
  labs(x = NULL, y = NULL) +
  nature_theme(5.0) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    plot.margin = margin(0, 2, 0, 12, "mm")
  )

ed2b_heat <- ggplot(heat2, aes(Protein, Trait, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  geom_text(
    aes(label = star),
    family = FONT_MYRIAD,
    size = pt_to_geom(TEXT_PT$heatmap_symbol),
    colour = COL_DARK
  ) +
  scale_fill_gradient2(
    low = COL_PURPLE, mid = "white", high = COL_GREEN, midpoint = 0,
    limits = c(-0.55, 0.55), oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  labs(title = "Extended protein–trait correlation landscape", x = NULL, y = NULL) +
  nature_theme(5.0) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = TEXT_PT$protein),
    axis.text.y = element_text(size = TEXT_PT$protein),
    axis.ticks = element_blank(), axis.line = element_blank(),
    legend.position = "right",
    legend.text = element_text(
      family = FONT_ARIAL,
      size = TEXT_PT$heatmap_colourbar_tick,
      face = "plain"
    ),
    plot.margin = margin(0, 2, 0, 2, "mm")
  )

ed2b <- patchwork::wrap_plots(
  ed2b_strip,
  ed2b_heat,
  ncol = 1,
  heights = c(0.13, 1),
  guides = "keep"
)

ed2 <- patchwork::wrap_plots(
  ed2a,
  ed2b,
  ncol = 1,
  heights = c(0.95, 1.05),
  guides = "keep"
) +
  patchwork::plot_annotation(tag_levels = "a") &
  tag_theme

ed2_files <- save_extended_figure(
  ed2,
  "Extended Data Fig. 2",
  "ExtendedDataFig2_DEP_RepresentativeProteins",
  ED_WIDTH_MM,
  ED2_HEIGHT_MM
)

save_panel(
  ed2a,
  "Extended Data Fig. 2",
  "a",
  "EDFig2a_representative_boxplots",
  180, 90
)
save_panel(
  ed2b,
  "Extended Data Fig. 2",
  "b",
  "EDFig2b_extended_heatmap",
  180, 100
)

# -----------------------------------------------------------------------------
# 9. Extended Data Figure 3: internal sensitivity and p-tau217 enrichment
# -----------------------------------------------------------------------------

apoe_labels <- c(
  "SPC25", "LRRN1", "AMBRA1", "PUM2",
  "NEFL", "PHGDH", "CAPN1", "USP5"
)

atn_labels <- c(
  "GDI1", "PLCH1", "C3", "PUM2", "PHGDH", "CAPN1"
)

cdrsb_labels <- c(
  "SPC25", "CPLX2", "SMOC1", "C3", "PUM2", "AMBRA1"
)

ptau_labels <- c(
  "SPC25", "LRRN1", "CPLX2", "SMOC1",
  "PLCH1", "C3", "PUM2", "GDI1"
)

ed3a <- volcano_panel(
  apoe_volcano,
  "logFC", "adj.P.Val", "EntrezGeneSymbol",
  APOE_ADJUSTED_TITLE,
  apoe_labels
)

ed3b <- effect_scatter_panel(
  apoe_compare,
  "logFC_primary", "logFC_apoe",
  "EntrezGeneSymbol", "adj.P.Val_primary",
  APOE_PRESERVATION_TITLE,
  apoe_labels,
  secondary_fdr_col = "adj.P.Val_apoe",
  y_lab = APOE_ADJUSTED_YLAB
)

ed3c <- volcano_panel(
  atn_volcano,
  "logFC", "adj.P.Val", "EntrezGeneSymbol",
  "Plasma AT(N)-adjusted model",
  c("GDI1")
)

ed3d <- effect_scatter_panel(
  atn_compare,
  "logFC_primary", "logFC_atn",
  "EntrezGeneSymbol", "adj.P.Val_primary",
  "Effect attenuation after plasma AT(N) adjustment",
  atn_labels,
  secondary_fdr_col = "adj.P.Val_atn",
  y_lab = "AT(N)-adjusted log2 fold-change"
)

ed3e <- volcano_panel(
  cdrsb_volcano,
  "logFC", "adj.P.Val", "EntrezGeneSymbol",
  "Within-AD CDR-SB severity model",
  character()
) +
  labs(x = "Within-AD CDR-SB slope")

ed3f <- effect_scatter_panel(
  cdrsb_compare,
  "logFC_primary", "beta_CDRSB_AD_only",
  "EntrezGeneSymbol", "adj.P.Val_primary",
  "Diagnostic effects versus within-AD severity",
  cdrsb_labels,
  secondary_fdr_col = "adj.P.Val_CDRSB_AD_only",
  y_lab = "Within-AD CDR-SB slope"
)

normalize_quantile <- function(x) {
  x <- as_num(x)
  dplyr::if_else(
    is.finite(x) & x > 1 & x <= 100,
    x / 100,
    x
  )
}

normalize_fraction_metric <- function(x) {
  x <- as_num(x)
  dplyr::if_else(
    is.finite(x) & x > 1 & x <= 100,
    x / 100,
    x
  )
}

ptau_main <- ptau_summary %>%
  dplyr::filter(
    as_chr(threshold_family) == "p_tau217_enriched"
  ) %>%
  dplyr::mutate(
    quantile = normalize_quantile(p_tau217_CN_quantile),
    Threshold = factor(
      paste0("p", round(100 * quantile)),
      levels = c("p80", "p90", "p95")
    ),
    fdr005_total = as_num(fdr005_total),
    n_CN_initial = as_num(n_CN_initial),
    n_AD_initial = as_num(n_AD_initial),
    correlation = normalize_fraction_metric(
      logFC_correlation_with_primary
    ),
    direction = normalize_fraction_metric(
      direction_consistency_all
    ),
    preservation = normalize_fraction_metric(
      preservation_fraction_primary_fdr005
    )
  ) %>%
  dplyr::filter(
    !is.na(Threshold),
    is.finite(quantile)
  ) %>%
  dplyr::arrange(Threshold)

# Panel g uses the model-summary counts directly.
ptau_counts_main <- ptau_main %>%
  dplyr::select(
    Threshold,
    n_CN_initial,
    n_AD_initial
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_CN_initial, n_AD_initial),
    names_to = "Group",
    values_to = "n"
  ) %>%
  dplyr::mutate(
    Group = dplyr::recode(
      Group,
      n_CN_initial = "CN",
      n_AD_initial = "AD"
    ),
    Group = factor(Group, levels = c("CN", "AD")),
    n = as_num(n)
  )

ptau_metrics <- ptau_main %>%
  dplyr::select(
    Threshold,
    correlation,
    direction,
    preservation
  ) %>%
  tidyr::pivot_longer(
    -Threshold,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Value = normalize_fraction_metric(Value),
    Metric = dplyr::recode(
      Metric,
      correlation = "Effect-size correlation",
      direction = "Direction consistency",
      preservation = "Primary FDR preservation"
    )
  )

ptau_gh_audit <- tibble::tibble(
  check = c(
    "Three p-tau217 thresholds are available",
    "Panel g sample counts are positive",
    "Panel h metrics are finite and non-zero",
    "Panel h metrics lie between 0 and 1"
  ),
  passed = c(
    dplyr::n_distinct(ptau_main$Threshold, na.rm = TRUE) == 3,
    nrow(ptau_counts_main) > 0 &&
      all(is.finite(ptau_counts_main$n)) &&
      all(ptau_counts_main$n > 0),
    nrow(ptau_metrics) > 0 &&
      all(is.finite(ptau_metrics$Value)) &&
      any(ptau_metrics$Value > 0),
    all(
      ptau_metrics$Value >= 0 &
        ptau_metrics$Value <= 1,
      na.rm = TRUE
    )
  )
)

if (!all(ptau_gh_audit$passed)) {
  stop(
    "Extended Data Fig. 3 panels g/h failed their data audit: ",
    paste(
      ptau_gh_audit$check[!ptau_gh_audit$passed],
      collapse = "; "
    ),
    call. = FALSE
  )
}

ed3g <- ggplot(
  ptau_counts_main,
  aes(Threshold, n, fill = Group)
) +
  geom_col(
    position = position_dodge(width = 0.70),
    width = 0.62,
    colour = COL_DARK,
    linewidth = 0.24
  ) +
  geom_text(
    aes(label = round(n)),
    position = position_dodge(width = 0.70),
    vjust = -0.35,
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    colour = COL_DARK
  ) +
  geom_text(
    data = ptau_main,
    aes(
      x = Threshold,
      y = Inf,
      label = paste0("FDR < 0.05: ", round(fdr005_total))
    ),
    inherit.aes = FALSE,
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    vjust = 1.25
  ) +
  scale_fill_manual(
    values = c("CN" = COL_BLUE, "AD" = COL_ORANGE),
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = "p-tau217 threshold-sweep sample sizes",
    x = "CN reference quantile",
    y = "Participants"
  ) +
  nature_theme(5.8) +
  theme(
    legend.position = "top",
    legend.justification = "left"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

ed3h <- ggplot(
  ptau_metrics,
  aes(
    Threshold,
    Value,
    group = Metric,
    colour = Metric,
    shape = Metric
  )
) +
  geom_line(linewidth = 0.50) +
  geom_point(size = 1.55, stroke = 0.30) +
  geom_text(
    aes(label = scales::percent(Value, accuracy = 1)),
    position = position_nudge(y = 0.035),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$compact_numeric),
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Effect-size correlation" = COL_DARK,
      "Direction consistency" = COL_BLUE,
      "Primary FDR preservation" = COL_ORANGE
    )
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.10))
  ) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(
    title = "p-tau217 threshold-sweep preservation",
    x = "CN reference quantile",
    y = "Metric",
    colour = NULL,
    shape = NULL
  ) +
  nature_theme(5.8) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(1.2, 2.0, 1.4, 1.8, "mm")
  ) +
  guides(
    colour = guide_legend(ncol = 1, byrow = TRUE),
    shape = guide_legend(ncol = 1, byrow = TRUE)
  )

ed3i <- effect_scatter_panel(
  ptau_p90, "logFC_primary", "logFC_secondary", "EntrezGeneSymbol", "adj.P.Val_primary",
  "Primary versus p-tau217-enriched p90 effects",
  labels = ptau_labels,
  secondary_fdr_col = "adj.P.Val_secondary",
  y_lab = "p-tau217-enriched p90 log2 fold-change",
  show_legend = TRUE
)

ed3_design <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 3),
  patchwork::area(t = 1, l = 4, b = 1, r = 6),
  patchwork::area(t = 2, l = 1, b = 2, r = 3),
  patchwork::area(t = 2, l = 4, b = 2, r = 6),
  patchwork::area(t = 3, l = 1, b = 3, r = 3),
  patchwork::area(t = 3, l = 4, b = 3, r = 6),
  patchwork::area(t = 4, l = 1, b = 4, r = 2),
  patchwork::area(t = 4, l = 3, b = 4, r = 4),
  patchwork::area(t = 4, l = 5, b = 4, r = 6)
)

ed3 <- patchwork::wrap_plots(
  ed3a, ed3b, ed3c, ed3d, ed3e, ed3f, ed3g, ed3h, ed3i,
  design = ed3_design,
  heights = c(1, 1, 1, 1.10),
  guides = "keep"
) +
  patchwork::plot_annotation(tag_levels = "a") &
  tag_theme

ed3_panel_list <- list(
  a = ed3a, b = ed3b, c = ed3c,
  d = ed3d, e = ed3e, f = ed3f,
  g = ed3g, h = ed3h, i = ed3i
)

ed3_files <- save_extended_figure(
  ed3,
  "Extended Data Fig. 3",
  "ExtendedDataFig3_DEP_InternalSensitivity_pTau217",
  ED_WIDTH_MM,
  ED3_HEIGHT_MM
)

purrr::iwalk(
  ed3_panel_list,
  ~ save_panel(
    .x,
    "Extended Data Fig. 3",
    .y,
    paste0("EDFig3", .y, "_internal_sensitivity_pTau217"),
    if (.y %in% letters[1:6]) 90 else 60,
    if (.y %in% letters[1:6]) 54 else 60
  )
)

# -----------------------------------------------------------------------------
# 10. Extended Data Figure 4: recruitment-context robustness and meta-analysis
# -----------------------------------------------------------------------------

# Apply the same public country-number labels used in Extended Data Fig. 1 to
# the LOSO panels. The raw site names remain available only in Source Data and
# are never used as plotting labels.
site_public_lookup <- site_anonymization %>%
  dplyr::select(Site_key, Site_public) %>%
  dplyr::filter(
    !is.na(Site_key),
    nzchar(Site_key),
    !is.na(Site_public),
    nzchar(Site_public)
  ) %>%
  dplyr::distinct(Site_key, .keep_all = TRUE)

loso_display <- loso %>%
  dplyr::mutate(
    excluded_site_raw = as_chr(excluded_site),
    Site_key = normalize_site_key(excluded_site_raw)
  ) %>%
  dplyr::left_join(
    site_public_lookup,
    by = "Site_key"
  ) %>%
  dplyr::mutate(
    excluded_site_public = dplyr::if_else(
      !is.na(Site_public) & nzchar(Site_public),
      Site_public,
      NA_character_
    )
  )

unmatched_loso_sites <- loso_display %>%
  dplyr::filter(
    is.na(excluded_site_public) |
      !nzchar(excluded_site_public)
  ) %>%
  dplyr::pull(excluded_site_raw) %>%
  unique()

if (length(unmatched_loso_sites) > 0) {
  stop(
    "Could not assign public country-number labels to LOSO sites: ",
    paste(unmatched_loso_sites, collapse = ", "),
    call. = FALSE
  )
}

ed4a <- robust_lollipop(
  loco,
  "excluded_country",
  "logFC_correlation",
  "LOCO effect-size concordance",
  reference = 0.90
)

ed4b <- robust_lollipop(
  loco,
  "excluded_country",
  "prop_main_sig_preserved",
  "LOCO preservation",
  percent = TRUE
)

ed4c <- robust_lollipop(
  loso_display,
  "excluded_site_public",
  "logFC_correlation",
  "LOSO effect-size concordance",
  reference = 0.90
)

ed4d <- robust_lollipop(
  loso_display,
  "excluded_site_public",
  "prop_main_sig_preserved",
  "LOSO preservation",
  percent = TRUE
)

ed4e <- balanced_hist(
  balanced,
  "logFC_correlation",
  "Balanced resampling: effect concordance",
  "Correlation with primary model"
)

ed4f <- balanced_hist(
  balanced,
  "direction_consistency_all",
  "Balanced resampling: direction consistency",
  "Direction consistency",
  percent = TRUE
)

# FDR-level preservation is retained as a numerical sensitivity result, but is
# not plotted because its near-zero distribution primarily reflects the severe
# loss of power under 8 CN + 8 AD per country rather than effect instability.
balanced_preservation_values <- as_num(
  balanced$prop_main_sig_preserved
)
balanced_preservation_values <- balanced_preservation_values[
  is.finite(balanced_preservation_values)
]

if (length(balanced_preservation_values) == 0) {
  stop(
    "Balanced-resampling FDR-preservation values are unavailable.",
    call. = FALSE
  )
}

balanced_preservation_summary <- tibble::tibble(
  metric = c(
    "iterations",
    "median_fraction",
    "median_percent",
    "nonzero_iterations",
    "nonzero_iteration_percent",
    "maximum_fraction",
    "maximum_percent"
  ),
  value = c(
    length(balanced_preservation_values),
    stats::median(balanced_preservation_values),
    100 * stats::median(balanced_preservation_values),
    sum(balanced_preservation_values > 0),
    100 * mean(balanced_preservation_values > 0),
    max(balanced_preservation_values),
    100 * max(balanced_preservation_values)
  )
)

balanced_preservation_summary_path <- file.path(
  diagnostics_root,
  "EDFig4_balanced_FDR_preservation_not_plotted.csv"
)
readr::write_csv(
  balanced_preservation_summary,
  balanced_preservation_summary_path
)

meta_significant <- country_meta %>%
  dplyr::mutate(
    meta_logFC = as_num(meta_logFC),
    meta_se = as_num(meta_se),
    meta_adj.P.Val = as_num(meta_adj.P.Val),
    I2 = as_num(I2),
    gene = as_chr(EntrezGeneSymbol),
    direction = dplyr::if_else(
      meta_logFC >= 0,
      "Higher in AD",
      "Lower in AD"
    )
  ) %>%
  dplyr::filter(
    meta_adj.P.Val < 0.05,
    is.finite(meta_logFC),
    is.finite(meta_se)
  )

if (nrow(meta_significant) == 0) {
  stop(
    "No BH-FDR-significant country meta-analysis signals were found.",
    call. = FALSE
  )
}

meta_i2_median <- stats::median(
  meta_significant$I2,
  na.rm = TRUE
)

if (!is.finite(meta_i2_median)) {
  meta_i2_median <- NA_real_
}

meta_selected <- meta_significant %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_min(
    meta_adj.P.Val,
    n = 7,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    lower = meta_logFC - 1.96 * meta_se,
    upper = meta_logFC + 1.96 * meta_se
  ) %>%
  dplyr::arrange(meta_logFC) %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene)
  )

meta_i2_label <- if (is.finite(meta_i2_median)) {
  sprintf("Median I² = %.0f%%", meta_i2_median)
} else {
  "Median I² unavailable"
}

ed4g <- ggplot(
  meta_selected,
  aes(meta_logFC, gene, colour = direction)
) +
  geom_vline(
    xintercept = 0,
    linewidth = LINE_MM,
    colour = COL_DARK
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0,
    orientation = "y",
    linewidth = 0.38
  ) +
  geom_point(size = 1.35) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = meta_i2_label,
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    colour = COL_DARK,
    hjust = 1.05,
    vjust = -0.65
  ) +
  scale_colour_manual(
    values = c(
      "Higher in AD" = COL_ORANGE,
      "Lower in AD" = COL_BLUE
    ),
    name = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Selected multicountry meta-analysis signals",
    x = "Random-effects meta-analytic log2 fold-change",
    y = NULL
  ) +
  nature_theme(5.8) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.y = element_text(size = TEXT_PT$protein),
    plot.margin = margin(1.5, 3.0, 3.5, 2.0, "mm")
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 0.8)
    )
  )

ed4_design <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 1),
  patchwork::area(t = 1, l = 2, b = 1, r = 2),
  patchwork::area(t = 2, l = 1, b = 2, r = 1),
  patchwork::area(t = 2, l = 2, b = 2, r = 2),
  patchwork::area(t = 3, l = 1, b = 3, r = 1),
  patchwork::area(t = 3, l = 2, b = 3, r = 2),
  patchwork::area(t = 4, l = 1, b = 4, r = 2)
)

ed4 <- patchwork::wrap_plots(
  ed4a,
  ed4b,
  ed4c,
  ed4d,
  ed4e,
  ed4f,
  ed4g,
  design = ed4_design,
  heights = c(1, 1, 1, 1.15),
  guides = "keep"
) +
  patchwork::plot_annotation(tag_levels = "a") &
  tag_theme

ed4_panel_list <- list(
  a = ed4a,
  b = ed4b,
  c = ed4c,
  d = ed4d,
  e = ed4e,
  f = ed4f,
  g = ed4g
)

ed4_files <- save_extended_figure(
  ed4,
  "Extended Data Fig. 4",
  "ExtendedDataFig4_DEP_RecruitmentRobustness",
  ED_WIDTH_MM,
  ED4_HEIGHT_MM
)

purrr::iwalk(
  ed4_panel_list,
  ~ save_panel(
    .x,
    "Extended Data Fig. 4",
    .y,
    paste0(
      "EDFig4",
      .y,
      "_recruitment_robustness"
    ),
    if (.y == "g") 180 else 90,
    if (.y == "g") 58 else 43
  )
)

# -----------------------------------------------------------------------------
# 11. Legends, Source Data crosswalk and submission audit
# -----------------------------------------------------------------------------

n_meta <- sum(
  as_num(country_meta$meta_adj.P.Val) < 0.05,
  na.rm = TRUE
)
n_meta_higher <- sum(
  as_num(country_meta$meta_adj.P.Val) < 0.05 &
    as_num(country_meta$meta_logFC) > 0,
  na.rm = TRUE
)
n_meta_lower <- sum(
  as_num(country_meta$meta_adj.P.Val) < 0.05 &
    as_num(country_meta$meta_logFC) < 0,
  na.rm = TRUE
)

legend_text <- c(
  paste0(
    "Extended Data Fig. 1 | Proteomic quality control, cohort structure ",
    "and statistical calibration. a–g, Principal component projections ",
    "coloured by sex, age, *APOE* ε4 carrier status, assay plate, country, ",
    "recruitment site and diagnostic group. Each legend is retained within ",
    "the corresponding PCA panel. h–j, Age, sex and *APOE* ε4 carrier-status ",
    "distributions by diagnostic group. k, Distribution of unadjusted ",
    "association P values. l, Quantile–quantile plot of observed and ",
    "expected −log10(P) values."
  ),
  paste0(
    "Extended Data Fig. 2 | Representative clinical AD-associated proteins ",
    "and the extended protein–trait landscape. a, Participant-level ",
    "non-residualized log2 RFU distributions for representative higher- and ",
    "lower-abundance proteins. Boxes show the interquartile range, horizontal ",
    "lines the median and whiskers 1.5 times the interquartile range. ",
    "b, Extended Spearman correlation heatmap across representative proteins ",
    "and diagnostic, demographic, clinical and plasma biomarker variables. ",
    "The upper strip denotes primary differential-abundance direction, the ",
    "local colour bar reports Spearman rho and asterisks denote BH-FDR-adjusted ",
    "significance."
  ),
  paste0(
    "Extended Data Fig. 3 | Internal sensitivity and plasma p-tau217 ",
    "enrichment analyses. a,b, *APOE* ε4-adjusted differential abundance and ",
    "effect preservation. c,d, Plasma AT(N)-adjusted differential abundance ",
    "and effect attenuation. e,f, Within-AD CDR-SB associations and their ",
    "alignment with primary diagnostic effects. g–i, Sample composition, ",
    "preservation metrics and primary-versus-secondary effect concordance ",
    "across cohort-defined p-tau217 enrichment thresholds. Legends are ",
    "retained within the sensitivity panel or paired analysis to which they ",
    "apply. Thresholds are exploratory enrichment anchors and not clinical ",
    "diagnostic cutoffs."
  ),
  sprintf(
    paste0(
      "Extended Data Fig. 4 | Recruitment-context robustness of the ",
      "clinical AD-associated plasma proteomic profile. a,b, Leave-one-country-",
      "out effect-size concordance and preservation of primary FDR-significant ",
      "proteins. c,d, Corresponding leave-one-site-out metrics using anonymized ",
      "country-number site labels. Dashed vertical lines in a and c mark ",
      "r = 0.90. e,f, Balanced country-resampling distributions for effect-size ",
      "correlation and directional consistency; dashed lines indicate medians. ",
      "The near-zero FDR-level preservation distribution from this severe ",
      "8 CN + 8 AD per-country stress test is reported in Source Data and the ",
      "diagnostic summary rather than plotted, because it primarily reflects ",
      "limited statistical power. g, Selected proteins from the country-level ",
      "random-effects meta-analysis with 95%% confidence intervals. A single ",
      "annotation reports the median I² across all BH-FDR-significant ",
      "meta-analytic signals. The full meta-analysis identified %d proteins at ",
      "BH-FDR < 0.05 (%d higher and %d lower)."
    ),
    n_meta,
    n_meta_higher,
    n_meta_lower
  )
)

legend_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_figure_legends_FINAL_DRAFT.txt"
)
writeLines(legend_text, legend_path, useBytes = TRUE)

crosswalk <- tibble::tribble(
  ~Display_item, ~Panel, ~Source_workbook, ~Source_sheet,
  "Extended Data Fig. 1", "a–g", pca_file,
  "PCA_scores; PCA_variance",
  "Extended Data Fig. 1", "h–j", demog_file,
  "Participant_data",
  "Extended Data Fig. 1", "k–l", pvalue_file,
  "Pvalue_distribution; QQ_plot",
  "Extended Data Fig. 2", "a", protein_file,
  "Boxplot_participant_data",
  "Extended Data Fig. 2", "b", protein_file,
  "Heatmap_protein_selection; Heatmap_correlations",
  "Extended Data Fig. 3", "a–f", sensitivity_file,
  paste0(
    "APOE_volcano; APOE_comparison; ATN_volcano; ATN_comparison; ",
    "CDRSB_volcano; CDRSB_comparison"
  ),
  "Extended Data Fig. 3", "g–i", ptau_file,
  "Threshold_summary; Sample_counts; P90_comparison",
  "Extended Data Fig. 4", "a–g", robust_file,
  "LOCO_summary; LOSO_summary; Balanced_iterations; Country_meta"
)

crosswalk_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_SourceData_crosswalk.csv"
)
readr::write_csv(crosswalk, crosswalk_path)

display_label_audit <- tibble::tibble(
  check = c(
    "CHI-24-002_ prefix absent from displayed plate labels",
    "Investigator/site names absent from all displayed site labels",
    "Multisite countries use country-number labels",
    "All LOSO sites map to public labels",
    "APOE is italicized in all displayed plot labels",
    "ED Fig. 1 panels g/h contain non-zero source values",
    "ED Fig. 3 panels g/h contain non-zero source values",
    "Near-zero balanced FDR preservation is not plotted",
    "Meta-analysis uses one global I2 annotation"
  ),
  passed = c(
    !any(
      stringr::str_detect(
        levels(pca_scores$Plate_label),
        "^CHI-24-002_"
      ),
      na.rm = TRUE
    ),
    !any(
      stringr::str_detect(
        c(
          levels(pca_scores$Site_label),
          as_chr(loso_display$excluded_site_public)
        ),
        stringr::regex(
          paste(
            c(
              "Avila", "Behrens", "Bruno", "Custodio",
              "Lopera", "Matallana", "Slachevsky"
            ),
            collapse = "|"
          ),
          ignore_case = TRUE
        )
      ),
      na.rm = TRUE
    ),
    all(
      site_anonymization$Site_public[
        site_anonymization$n_sites_country > 1
      ] == paste0(
        site_anonymization$Country_text[
          site_anonymization$n_sites_country > 1
        ],
        " ",
        site_anonymization$site_number[
          site_anonymization$n_sites_country > 1
        ]
      )
    ),
    all(
      !is.na(loso_display$excluded_site_public) &
        nzchar(loso_display$excluded_site_public)
    ),
    all(vapply(
      list(
        APOE_PCA_TITLE,
        APOE_CARRIER_TITLE,
        APOE_ADJUSTED_TITLE,
        APOE_PRESERVATION_TITLE,
        APOE_ADJUSTED_YLAB
      ),
      is.expression,
      logical(1)
    )),
    all(ed1_gh_audit$passed),
    all(ptau_gh_audit$passed),
    !exists("ed4_primary_preservation_plot", inherits = FALSE) &&
      identical(names(ed4_panel_list), letters[1:7]),
    length(meta_i2_label) == 1 &&
      !any(grepl("I²=", unlist(lapply(ed4g$layers, function(l) c(l$aes_params$label, if (!is.null(l$data) && "label" %in% names(l$data)) as.character(l$data$label))))))
  )
)

if (!all(display_label_audit$passed)) {
  stop(
    "Display/data audit failed: ",
    paste(
      display_label_audit$check[!display_label_audit$passed],
      collapse = "; "
    ),
    call. = FALSE
  )
}

display_label_audit_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_display_and_zero_audit.csv"
)
readr::write_csv(
  display_label_audit,
  display_label_audit_path
)

illustrator_guide_path <- file.path(
  illustrator_root,
  "README_OPEN_EXTENDED_DATA_IN_ADOBE_ILLUSTRATOR.txt"
)

writeLines(
  c(
    "DEP EXTENDED DATA FIGURES — ADOBE ILLUSTRATOR WORKFLOW",
    "",
    "Preferred files:",
    "- Open complete figures from Illustrator_ready/master_vectors/*.svg.",
    "- Open individual panels from Illustrator_ready/panels_vector/*.svg.",
    "- Save the opened SVG immediately as a native Adobe Illustrator .ai file.",
    "",
    "Checks after opening:",
    "1. Confirm the complete-figure artboard remains 180 mm wide.",
    "2. Confirm Arial is available before opening, so live text is preserved.",
    "3. Confirm every legend remains inside its corresponding panel.",
    "4. Do not globally release clipping masks; release only a local group when",
    "   an intentional manual adjustment is required.",
    "5. Keep panel letters at 7 pt and panel titles at 6.5 pt.",
    "6. Do not reduce any body, axis or legend text below 5 pt.",
    "9. Keep orange #F46D43 and blue #3288BD unchanged.",
    "10. Export the final submission TIFF at 300 dpi in RGB with LZW compression.",
    "11. Do not reinsert the near-zero balanced FDR-preservation histogram;",
    "   that stress-test result is retained numerically in diagnostics.",
    "12. Keep only one global median I² annotation in Extended Data Fig. 4g.",
    "",
    "PDF files are vector backups. PNG files are review previews only."
  ),
  illustrator_guide_path,
  useBytes = TRUE
)

# Verify that no display item lost all outputs.
expected_complete <- c(
  "Extended Data Fig. 1",
  "Extended Data Fig. 2",
  "Extended Data Fig. 3",
  "Extended Data Fig. 4"
)

complete_registry <- OUTPUT_REGISTRY %>%
  dplyr::filter(panel == "complete figure")

missing_complete <- setdiff(
  expected_complete,
  unique(complete_registry$display_item)
)

if (length(missing_complete) > 0) {
  stop(
    "Complete Extended Data outputs are missing for: ",
    paste(missing_complete, collapse = ", "),
    call. = FALSE
  )
}

if (!all(file.exists(OUTPUT_REGISTRY$path))) {
  stop(
    "The output registry contains missing files: ",
    paste(
      OUTPUT_REGISTRY$file[!file.exists(OUTPUT_REGISTRY$path)],
      collapse = ", "
    ),
    call. = FALSE
  )
}

manifest <- OUTPUT_REGISTRY %>%
  dplyr::mutate(
    exists = file.exists(path),
    size_bytes = file.info(path)$size,
    size_MB = round(size_bytes / 1024^2, 3)
  ) %>%
  dplyr::arrange(display_item, panel, format)

manifest_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_figure_manifest.csv"
)
readr::write_csv(manifest, manifest_path)

warnings_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_warnings.txt"
)
warnings_object <- warnings()
if (is.null(warnings_object)) {
  writeLines("No warnings were recorded.", warnings_path)
} else {
  writeLines(capture.output(warnings_object), warnings_path)
}

session_path <- file.path(
  diagnostics_root,
  "DEP_ExtendedData_sessionInfo.txt"
)
writeLines(
  capture.output(utils::sessionInfo()),
  session_path
)

pipeline_success <- (
  length(missing_complete) == 0 &&
    all(file.exists(OUTPUT_REGISTRY$path))
)

if (!isTRUE(pipeline_success)) {
  stop(
    "The final package audit did not pass. ",
    "Do not use the generated outputs.",
    call. = FALSE
  )
}

message("\n============================================================")
message(
  "Fixed Extended Data artboards confirmed: 180 mm width; heights ",
  paste(c(ED1_HEIGHT_MM, ED2_HEIGHT_MM, ED3_HEIGHT_MM, ED4_HEIGHT_MM), collapse = ", "),
  " mm."
)
message("DEP Extended Data package generated successfully.")
message("Output root: ", out_root)
message("Submission TIFFs: ", submission_root)
message("Illustrator master vectors: ", master_vector_root)
message("Illustrator panel vectors: ", panel_vector_root)
message("Legends: ", legend_path)
message("Source Data crosswalk: ", crosswalk_path)
message("Display/zero audit: ", display_label_audit_path)
message(
  "Balanced FDR-preservation summary: ",
  balanced_preservation_summary_path
)
message("Manifest: ", manifest_path)
message("Warnings: ", warnings_path)
message("============================================================\n")

###############################################################################
# END
###############################################################################

