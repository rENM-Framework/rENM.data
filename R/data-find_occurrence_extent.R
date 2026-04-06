#' Find the spatial extent of a species' eBird occurrence data
#'
#' Computes a centered percentile bounding box from occurrence files and
#' writes the result to an extent.txt file for use in downstream rENM
#' preprocessing and modeling workflows.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' Occurrence files are expected in:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/}
#'
#' with filenames of the form:
#'
#' \code{of-<year>.csv}
#'
#' Each file must contain coordinate columns equivalent to longitude and
#' latitude. Common aliases are automatically detected.
#'
#' \strong{Methods}:
#' A centered percentile bounding box is computed by selecting the central
#' \code{bbox_pct} percent of valid points along each axis and discarding
#' symmetric tails. For example:
#' \itemize{
#'   \item bbox_pct = 95 uses the (2.5, 97.5) quantiles
#'   \item bbox_pct = 100 uses the full min and max range
#' }
#'
#' All valid coordinates across all files are pooled prior to computing
#' quantiles.
#'
#' \strong{Coordinate validation}:
#' For each file, the function:
#' \itemize{
#'   \item Locates longitude and latitude columns using common aliases:
#'         longitude, lon, x and latitude, lat, y
#'   \item Converts values to numeric and discards invalid rows
#'   \item Retains only coordinates within:
#'         longitude in (-180, 180)
#'         latitude in (-90, 90)
#'   \item Removes zero-zero coordinates (0, 0)
#' }
#'
#' \strong{Outputs}:
#' The resulting bounding box is written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/extent.txt}
#'
#' with the format:
#' \itemize{
#'   \item Upper-left:  (<lon>, <lat>)
#'   \item Lower-right: (<lon>, <lat>)
#' }
#'
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' using the section title:
#'
#' \code{Processing summary (find_occurrence_extent)}
#'
#' \strong{Workflow context}:
#' This function is intended to be used after occurrence preprocessing
#' steps such as:
#' \itemize{
#'   \item get_ebird_occurrences()
#'   \item remove_duplicate_occurrences()
#'   \item thin_occurrences() or thin_occurrences2()
#'   \item tidy_occurrences()
#' }
#'
#' At this stage, the _occs directory should contain one or more
#' \code{of-<year>.csv} files with species and coordinate columns.
#'
#' @param alpha_code Character. Four-letter banding code for the species.
#'   Used to construct the run directory path.
#' @param bbox_pct Numeric. Percent of points to include in the centered
#'   bounding box with (0 < bbox_pct <= 100). Values less than 100 discard
#'   symmetric tails along each axis. Default = 99.
#'
#' @return Invisibly returns a List with the following components:
#' \itemize{
#'   \item alpha_code: Character species code
#'   \item bbox:
#'     \itemize{
#'       \item xmin: Minimum longitude
#'       \item xmax: Maximum longitude
#'       \item ymin: Minimum latitude
#'       \item ymax: Maximum latitude
#'     }
#'   \item paths:
#'     \itemize{
#'       \item extent: File path to extent.txt
#'       \item log: File path to the run log
#'     }
#'   \item counts:
#'     \itemize{
#'       \item files: Number of occurrence files scanned
#'       \item rows_read: Total rows read across all files
#'       \item points_used: Number of valid coordinate rows used
#'       \item bbox_pct: Percentile used for bounding box
#'     }
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes extent.txt to the _occs directory
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom utils read.csv
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#'   res_full <- find_occurrence_extent("CASP")
#'   res_full$bbox
#'   res_full$paths$extent
#'
#'   res_95 <- find_occurrence_extent("CASP", bbox_pct = 95)
#'   res_95$bbox
#' }
#'
#' @export
find_occurrence_extent <- function(alpha_code, bbox_pct = 99) {

  ## ---- validation ------------------------------------------------------------
  stopifnot(
    is.character(alpha_code),
    length(alpha_code) == 1,
    nzchar(alpha_code)
  )
  stopifnot(
    is.numeric(bbox_pct),
    length(bbox_pct) == 1,
    is.finite(bbox_pct)
  )
  if (bbox_pct <= 0 || bbox_pct > 100) {
    stop("bbox_pct must be in (0, 100].")
  }

  ## ---- helpers ---------------------------------------------------------------
  .expand   <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  .mkdir    <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  .now      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep_line <- function(n = 72L) paste(rep.int("-", n), collapse = "")
  .catln    <- function(...) { cat(paste0(..., "\n")); flush.console() }
  .norm     <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
  .find_col <- function(nms, candidates) {
    if (!length(nms)) return(NA_character_)
    nn <- .norm(nms)
    for (cand in candidates) {
      j <- match(.norm(cand), nn, nomatch = 0)
      if (j > 0) return(nms[j])
    }
    NA_character_
  }

  .write_log_section <- function(logfile, title, lines) {
    dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
    section <- c(
      "",
      .sep_line(),
      title,
      paste0("Timestamp: ", .now()),
      lines,
      ""
    )
    cat(paste0(section, collapse = "\n"), file = logfile, append = TRUE)
  }

  ## ---- paths -----------------------------------------------------------------
  project_dir <- rENM_project_dir()
  run_dir   <- .expand(file.path(project_dir, "runs", alpha_code))
  occs_dir  <- file.path(run_dir, "_occs")
  extent_fp <- file.path(occs_dir, "extent.txt")
  log_fp    <- file.path(run_dir, "_log.txt")
  .mkdir(occs_dir)

  ## ---- discover files --------------------------------------------------------
  patt  <- "^of-\\d+\\.csv$"
  files <- list.files(occs_dir, pattern = patt, full.names = TRUE)
  if (!length(files)) {
    stop("No of-<year>.csv files found in: ", occs_dir)
  }

  ## ---- console banner --------------------------------------------------------
  .catln(.sep_line())
  .catln(
    sprintf(
      "find_occurrence_extent: alpha_code=%s | bbox_pct=%.2f%%",
      alpha_code,
      bbox_pct
    )
  )
  .catln("Scanning directory: ", occs_dir)
  .catln(
    sprintf(
      "Found %d file(s): %s",
      length(files),
      paste(basename(files), collapse = ", ")
    )
  )
  .catln(.sep_line())

  ## ---- read and collect valid coordinates -----------------------------------
  all_lon <- numeric(0)
  all_lat <- numeric(0)
  total_rows <- 0L
  used_rows  <- 0L

  for (fp in files) {
    df <- try(
      utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE),
      silent = TRUE
    )
    if (inherits(df, "try-error")) {
      .catln("  ! Failed to read: ", fp)
      next
    }

    lon_col <- .find_col(names(df), c("longitude", "lon", "x"))
    lat_col <- .find_col(names(df), c("latitude", "lat", "y"))
    if (is.na(lon_col) || is.na(lat_col)) {
      .catln("  ! Skipping (missing lon/lat): ", basename(fp))
      next
    }

    lon <- suppressWarnings(as.numeric(df[[lon_col]]))
    lat <- suppressWarnings(as.numeric(df[[lat_col]]))
    total_rows <- total_rows + length(lon)

    ok <- is.finite(lon) & is.finite(lat) &
      lon >= -180 & lon <= 180 &
      lat >= -90  & lat <= 90 &
      !(lon == 0 & lat == 0)

    all_lon <- c(all_lon, lon[ok])
    all_lat <- c(all_lat, lat[ok])
    used_rows <- used_rows + sum(ok)
  }

  if (!length(all_lon)) {
    stop("No valid coordinates found across the of-*.csv files in: ", occs_dir)
  }

  .catln(
    sprintf(
      "Valid coordinate rows used: %d (of %d total rows read)",
      used_rows,
      total_rows
    )
  )

  ## ---- compute centered percentile bbox -------------------------------------
  p <- bbox_pct / 100
  tail_p <- (1 - p) / 2
  q_lo <- max(0, tail_p)
  q_hi <- min(1, 1 - tail_p)

  q_lon <- as.numeric(
    stats::quantile(
      all_lon,
      probs = c(q_lo, q_hi),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
  q_lat <- as.numeric(
    stats::quantile(
      all_lat,
      probs = c(q_lo, q_hi),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  lon_min <- min(q_lon)
  lon_max <- max(q_lon)
  lat_min <- min(q_lat)
  lat_max <- max(q_lat)

  ul_lon <- lon_min
  ul_lat <- lat_max
  lr_lon <- lon_max
  lr_lat <- lat_min

  ## ---- write extent.txt ------------------------------------------------------
  header <- c(
    "# extent.txt",
    "# Computed by find_occurrence_extent()",
    sprintf("# Timestamp: %s", .now()),
    sprintf(
      "# Source: %d file(s) in _occs (of-<year>.csv)",
      length(files)
    ),
    sprintf(
      "# Points used: %d of %d total rows",
      used_rows,
      total_rows
    ),
    sprintf(
      "# bbox_pct: %.2f%% (centered percentile box)",
      bbox_pct
    ),
    "# Coordinate order: (lon, lat)"
  )
  lines <- c(
    header,
    sprintf("Upper-left:  (%.6f, %.6f)", ul_lon, ul_lat),
    sprintf("Lower-right: (%.6f, %.6f)", lr_lon, lr_lat)
  )
  writeLines(lines, con = extent_fp)

  .catln("Wrote extent.txt -> ", extent_fp)
  .catln(sprintf("  Upper-left:  (%.6f, %.6f)", ul_lon, ul_lat))
  .catln(sprintf("  Lower-right: (%.6f, %.6f)", lr_lon, lr_lat))

  ## ---- append Processing summary to log -------------------------------------
  summary_lines <- c(
    sprintf("Alpha code:          %s", alpha_code),
    sprintf("Bounding box pct:    %.2f%%", bbox_pct),
    sprintf("Files scanned:       %d", length(files)),
    sprintf("Rows read (total):   %d", total_rows),
    sprintf("Points used:         %d", used_rows),
    "Extent (lon/lat):",
    sprintf("  Upper-left:        (%.6f, %.6f)", ul_lon, ul_lat),
    sprintf("  Lower-right:       (%.6f, %.6f)", lr_lon, lr_lat),
    sprintf("Output:              %s", extent_fp)
  )

  .write_log_section(
    logfile = log_fp,
    title   = "Processing summary (find_occurrence_extent)",
    lines   = summary_lines
  )

  .catln(.sep_line())
  .catln("Done. Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(list(
    alpha_code = alpha_code,
    bbox = list(
      xmin = lon_min,
      xmax = lon_max,
      ymin = lat_min,
      ymax = lat_max
    ),
    paths = list(
      extent = extent_fp,
      log    = log_fp
    ),
    counts = list(
      files       = length(files),
      rows_read   = total_rows,
      points_used = used_rows,
      bbox_pct    = bbox_pct
    )
  ))
}
