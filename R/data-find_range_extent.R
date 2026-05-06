#' Derive a spatial extent from a species USGS GAP range polygon
#'
#' Reads a USGS GAP range shapefile for a species, computes a padded
#' geographic bounding box in WGS84 coordinates, and writes the result
#' to \code{_occs/extent.txt} within the species run directory.
#'
#' Existing \code{extent.txt} files are automatically backed up prior
#' to writing a new extent definition.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   matching the \code{ALPHA.CODE} field in the project species table.
#'
#' @param pad_pct Numeric scalar. Symmetric percent padding applied to
#'   the bounding box. Positive values expand the extent and negative
#'   values contract it. Must be greater than -100. Default is
#'   \code{2}.
#'
#' @param project_dir Character scalar or \code{NULL}. Path to the
#'   rENM project root. If \code{NULL}, the project directory is
#'   resolved using \code{rENM_project_dir()}.
#'
#' @details
#' The species lookup table is expected at \code{data/_species.csv}
#' under the project root directory.
#'
#' The function identifies the appropriate GAP range identifier using
#' the \code{GAP.RANGE} field associated with the supplied
#' \code{alpha_code}. The corresponding shapefile is expected at
#' \code{data/shapefiles/GAP_RANGE/GAP_RANGE.shp} under the project
#' root, where \code{GAP_RANGE} is the value of that field.
#'
#' If the shapefile coordinate reference system is not longitude/latitude
#' WGS84 (EPSG:4326), it is transformed before extent calculation.
#' The resulting bounding box is symmetrically padded relative to its
#' center using \code{pad_pct}. The result is written to
#' \code{runs/ALPHA_CODE/_occs/extent.txt}.
#'
#' @return Invisibly returns a named list with numeric elements
#'   \code{ul_lon}, \code{ul_lat}, \code{lr_lon}, and \code{lr_lat}
#'   (upper-left and lower-right corners in decimal degrees).
#'
#' @seealso
#'   \code{\link{find_occurrence_extent}},
#'   \code{\link{set_extent}},
#'   \code{\link{get_merra_variables}}
#'
#' @importFrom sf st_bbox st_is_longlat st_make_valid st_read st_transform
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' res <- find_range_extent("CASP")
#'
#' find_range_extent(
#'   alpha_code = "CASP",
#'   pad_pct = 5
#' )
#' }
#'
#' @export
find_range_extent <- function(alpha_code, pad_pct = 2, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(pad_pct) || length(pad_pct) != 1L || !is.finite(pad_pct)) {
    stop("`pad_pct` must be a finite numeric scalar.", call. = FALSE)
  }
  if (pad_pct <= -100) {
    stop("`pad_pct` must be greater than -100.", call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required but not installed.", call. = FALSE)
  }

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir   <- file.path(project_root, "runs", alpha_code)
  occs_dir  <- file.path(run_dir, "_occs")
  log_fp    <- file.path(run_dir, "_log.txt")
  extent_fp <- file.path(occs_dir, "extent.txt")
  .mkdir(occs_dir)

  ## ---- read species table and resolve GAP range -----------------------------
  species_csv <- .expand(file.path(project_root, "data", "_species.csv"))
  if (!file.exists(species_csv)) {
    stop("Species table not found at: ", species_csv, call. = FALSE)
  }

  sp <- utils::read.csv(species_csv, stringsAsFactors = FALSE, check.names = FALSE)

  # Normalize column names for lookup (matches rENM.core convention)
  norm      <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
  col_norm  <- norm(names(sp))
  alpha_idx <- match("ALPHACODE", col_norm, nomatch = 0L)
  gap_idx   <- match("GAPRANGE",  col_norm, nomatch = 0L)

  if (alpha_idx == 0L || gap_idx == 0L) {
    stop(
      "Species table must contain ALPHA.CODE and GAP.RANGE columns.\n",
      "Available: ", paste(names(sp), collapse = ", "),
      call. = FALSE
    )
  }

  alpha_col <- names(sp)[alpha_idx]
  gap_col   <- names(sp)[gap_idx]

  row_idx <- which(toupper(trimws(sp[[alpha_col]])) == toupper(trimws(alpha_code)))
  if (!length(row_idx)) {
    stop("No row found with alpha code '", alpha_code, "' in: ", species_csv, call. = FALSE)
  }

  gap_range <- trimws(sp[[gap_col]][row_idx[1L]])
  .catln("[find_range_extent] GAP.RANGE='", gap_range, "'")

  ## ---- locate shapefile -----------------------------------------------------
  shp_path <- .expand(file.path(project_root, "data", "shapefiles",
                                gap_range, paste0(gap_range, ".shp")))
  if (!file.exists(shp_path)) {
    stop("Shapefile not found: ", shp_path, call. = FALSE)
  }
  .catln("[find_range_extent] Reading: ", shp_path)

  ## ---- read, validate CRS, compute bbox -------------------------------------
  sfobj <- sf::st_read(shp_path, quiet = TRUE)

  if (sf::st_is_longlat(sfobj)) {
    geom_ll <- sfobj
    .catln("[find_range_extent] CRS is lon/lat; using as-is.")
  } else {
    .catln("[find_range_extent] Projecting to EPSG:4326 (WGS84).")
    geom_ll <- sf::st_transform(sfobj, 4326L)
  }

  bb <- sf::st_bbox(sf::st_make_valid(geom_ll))

  xmin <- as.numeric(bb[["xmin"]]); xmax <- as.numeric(bb[["xmax"]])
  ymin <- as.numeric(bb[["ymin"]]); ymax <- as.numeric(bb[["ymax"]])

  if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) {
    stop("Non-finite bounding box coordinates from shapefile.", call. = FALSE)
  }

  width  <- xmax - xmin
  height <- ymax - ymin

  if (width <= 0 || height <= 0) {
    stop("Invalid bounding box: non-positive width or height.", call. = FALSE)
  }

  ## ---- apply symmetric padding ----------------------------------------------
  factor     <- 1 + pad_pct / 100
  cx         <- (xmin + xmax) / 2
  cy         <- (ymin + ymax) / 2
  ul_lon     <- cx - (width  * factor) / 2
  lr_lon     <- cx + (width  * factor) / 2
  ul_lat     <- cy + (height * factor) / 2
  lr_lat     <- cy - (height * factor) / 2

  .catln(sprintf("[find_range_extent] Padded bbox: UL=(%.6f, %.6f) LR=(%.6f, %.6f)",
                 ul_lon, ul_lat, lr_lon, lr_lat))

  ## ---- backup and write extent.txt ------------------------------------------
  if (file.exists(extent_fp)) {
    ts         <- gsub("[^0-9]", "", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    backup_path <- file.path(occs_dir, sprintf("extent_backup_%s.txt", ts))
    file.rename(extent_fp, backup_path)
    .catln("[find_range_extent] Backed up existing extent.txt to: ", basename(backup_path))
  }

  writeLines(c(
    "# extent.txt",
    "# Computed by find_range_extent()",
    paste0("# Timestamp: ", .now()),
    paste0("# Source: ", basename(shp_path), " via find_range_extent()"),
    "# Points used: n/a",
    paste0("# pad_pct: ", formatC(pad_pct, format = "f", digits = 1),
           "% (symmetric percent padding)"),
    "# Coordinate order: (lon, lat)",
    sprintf("Upper-left:  (%.6f, %.6f)", ul_lon, ul_lat),
    sprintf("Lower-right: (%.6f, %.6f)", lr_lon, lr_lat)
  ), con = extent_fp)

  .catln("[find_range_extent] Wrote extent.txt: ", extent_fp)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (find_range_extent)", c(
    sprintf("Alpha code:   %s", alpha_code),
    sprintf("Species CSV:  %s", species_csv),
    sprintf("GAP.RANGE:    %s", gap_range),
    sprintf("Shapefile:    %s", shp_path),
    sprintf("Output file:  %s", extent_fp),
    sprintf("Pad pct:      %.1f%%", pad_pct),
    sprintf("UL (lon,lat): (%.6f, %.6f)", ul_lon, ul_lat),
    sprintf("LR (lon,lat): (%.6f, %.6f)", lr_lon, lr_lat)
  ))

  invisible(list(
    ul_lon = ul_lon,
    ul_lat = ul_lat,
    lr_lon = lr_lon,
    lr_lat = lr_lat
  ))
}
