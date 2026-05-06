#' Derive a spatial extent from occurrence data
#'
#' Pools valid coordinates from all per-bin occurrence CSV files in
#' \code{_occs/}, computes a centered percentile bounding box, and writes the
#' result to \code{_occs/extent.txt} for use in downstream rENM workflows.
#'
#' @details
#' \strong{Percentile bounding box}
#'
#' The bounding box is computed by selecting the central \code{bbox_pct}
#' percent of valid points along each axis independently, discarding symmetric
#' tails. For example, \code{bbox_pct = 95} uses the 2.5th and 97.5th
#' quantiles; \code{bbox_pct = 100} uses the full coordinate range.
#'
#' \strong{Coordinate validation}
#'
#' For each file, longitude and latitude columns are detected automatically
#' from common aliases (\code{longitude}/\code{lon}/\code{x} and
#' \code{latitude}/\code{lat}/\code{y}). Values are coerced to numeric;
#' records outside the valid range (longitude -180 to 180, latitude -90 to 90)
#' or exactly at (0, 0)
#' are excluded. All surviving coordinates are pooled before quantile
#' computation.
#'
#' This function is typically called after \code{\link{tidy_occurrences}} has
#' moved finalized occurrence files into \code{_occs/}.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param bbox_pct Numeric scalar. Central percentage of points to include in
#'   the bounding box. Must be in (0, 100]. Default is 99.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{bbox} (list with \code{xmin}, \code{xmax}, \code{ymin},
#'   \code{ymax}), \code{paths} (list with \code{extent} and \code{log} file
#'   paths), and \code{counts} (list with \code{files}, \code{rows_read},
#'   \code{points_used}, \code{bbox_pct}).
#'
#' @seealso \code{\link{tidy_occurrences}}, \code{\link{find_range_extent}},
#'   \code{\link{set_extent}}, \code{\link{get_merra_variables}}
#'
#' @importFrom utils read.csv
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#' res <- find_occurrence_extent("CASP")
#' res$bbox
#'
#' res95 <- find_occurrence_extent("CASP", bbox_pct = 95)
#' find_occurrence_extent("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @export
find_occurrence_extent <- function(alpha_code, bbox_pct = 99, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(bbox_pct) || length(bbox_pct) != 1L || !is.finite(bbox_pct) ||
      bbox_pct <= 0 || bbox_pct > 100) {
    stop("`bbox_pct` must be a finite numeric scalar in (0, 100].", call. = FALSE)
  }

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir   <- .expand(file.path(project_root, "runs", alpha_code))
  occs_dir  <- file.path(run_dir, "_occs")
  extent_fp <- file.path(occs_dir, "extent.txt")
  log_fp    <- file.path(run_dir, "_log.txt")
  .mkdir(occs_dir)

  files <- list.files(occs_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No of-<year>.csv files found in: ", occs_dir, call. = FALSE)
  }

  ## ---- banner ---------------------------------------------------------------
  .catln(.sep_line())
  .catln(sprintf("find_occurrence_extent: alpha_code=%s | bbox_pct=%.2f%%",
                 alpha_code, bbox_pct))
  .catln("Scanning: ", occs_dir)
  .catln(sprintf("Found %d file(s): %s", length(files),
                 paste(basename(files), collapse = ", ")))
  .catln(.sep_line())

  ## ---- collect valid coordinates --------------------------------------------
  all_lon    <- numeric(0)
  all_lat    <- numeric(0)
  total_rows <- 0L
  used_rows  <- 0L

  for (fp in files) {
    df <- tryCatch(
      utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) { .catln("  ! Failed to read: ", fp); NULL }
    )
    if (is.null(df)) next

    lon_col <- .find_col(names(df), c("longitude", "lon", "x"))
    lat_col <- .find_col(names(df), c("latitude",  "lat", "y"))

    if (is.na(lon_col) || is.na(lat_col)) {
      .catln("  ! Skipping (no lon/lat columns): ", basename(fp))
      next
    }

    lon        <- suppressWarnings(as.numeric(df[[lon_col]]))
    lat        <- suppressWarnings(as.numeric(df[[lat_col]]))
    total_rows <- total_rows + length(lon)

    ok <- is.finite(lon) & is.finite(lat) &
      lon >= -180 & lon <= 180 &
      lat >=  -90 & lat <=  90 &
      !(lon == 0 & lat == 0)

    all_lon   <- c(all_lon, lon[ok])
    all_lat   <- c(all_lat, lat[ok])
    used_rows <- used_rows + sum(ok)
  }

  if (!length(all_lon)) {
    stop("No valid coordinates found across of-*.csv files in: ", occs_dir, call. = FALSE)
  }

  .catln(sprintf("Valid rows used: %d (of %d total)", used_rows, total_rows))

  ## ---- centered percentile bbox ---------------------------------------------
  p      <- bbox_pct / 100
  tail_p <- (1 - p) / 2
  q_lo   <- max(0, tail_p)
  q_hi   <- min(1, 1 - tail_p)

  q_lon <- as.numeric(stats::quantile(all_lon, probs = c(q_lo, q_hi),
                                      na.rm = TRUE, names = FALSE, type = 7))
  q_lat <- as.numeric(stats::quantile(all_lat, probs = c(q_lo, q_hi),
                                      na.rm = TRUE, names = FALSE, type = 7))

  lon_min <- min(q_lon); lon_max <- max(q_lon)
  lat_min <- min(q_lat); lat_max <- max(q_lat)

  ## ---- write extent.txt -----------------------------------------------------
  writeLines(c(
    "# extent.txt",
    "# Computed by find_occurrence_extent()",
    sprintf("# Timestamp: %s", .now()),
    sprintf("# Source: %d file(s) in _occs/", length(files)),
    sprintf("# Points used: %d of %d total rows", used_rows, total_rows),
    sprintf("# bbox_pct: %.2f%% (centered percentile box)", bbox_pct),
    "# Coordinate order: (lon, lat)",
    sprintf("Upper-left:  (%.6f, %.6f)", lon_min, lat_max),
    sprintf("Lower-right: (%.6f, %.6f)", lon_max, lat_min)
  ), con = extent_fp)

  .catln("Wrote extent.txt: ", extent_fp)
  .catln(sprintf("  Upper-left:  (%.6f, %.6f)", lon_min, lat_max))
  .catln(sprintf("  Lower-right: (%.6f, %.6f)", lon_max, lat_min))

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (find_occurrence_extent)", c(
    sprintf("Alpha code:        %s", alpha_code),
    sprintf("Bounding box pct:  %.2f%%", bbox_pct),
    sprintf("Files scanned:     %d", length(files)),
    sprintf("Rows read (total): %d", total_rows),
    sprintf("Points used:       %d", used_rows),
    "Extent (lon/lat):",
    sprintf("  Upper-left:      (%.6f, %.6f)", lon_min, lat_max),
    sprintf("  Lower-right:     (%.6f, %.6f)", lon_max, lat_min),
    sprintf("Output:            %s", extent_fp)
  ))

  .catln("Done. Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(list(
    alpha_code = alpha_code,
    bbox       = list(xmin = lon_min, xmax = lon_max,
                      ymin = lat_min, ymax = lat_max),
    paths      = list(extent = extent_fp, log = log_fp),
    counts     = list(files       = length(files),
                      rows_read   = total_rows,
                      points_used = used_rows,
                      bbox_pct    = bbox_pct)
  ))
}
