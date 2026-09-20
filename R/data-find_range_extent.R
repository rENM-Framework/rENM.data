#' Derive a spatial extent from a species USGS GAP range polygon
#'
#' Reads a USGS GAP range shapefile for a species, buffers the range
#' polygon outward by a fixed real-world margin, and writes the buffered
#' polygon's geographic bounding box to \code{_occs/extent.txt} within the
#' species run directory. The buffered polygon itself is saved alongside
#' \code{extent.txt} for downstream use.
#'
#' Existing \code{extent.txt} files are automatically backed up prior
#' to writing a new extent definition.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   matching the \code{ALPHA.CODE} field in the project species table.
#'
#' @param buffer_km Numeric scalar. Outward buffer applied to the GAP
#'   range polygon, in kilometers, before its bounding box is taken.
#'   Must be non-negative. Default is \code{250}.
#'
#'   The buffered polygon is species-specific. The distance is not. The same
#'   250 km is used for every species. Sizing it from a species' own computed
#'   centroid velocity would be circular, since that velocity is an output of
#'   the analysis rather than a known input to it. The 250 km
#'   figure derives from Huang, Sauer & Dubayah (2017), who tracked the
#'   abundance-weighted geographic centroid of 57 permanent resident North
#'   American bird species over 44 years using Breeding Bird Survey data
#'   and found an average, multidirectional centroid velocity of
#'   5.89 km/yr; projected across rENM's 45-year study window that gives
#'   approximately 265 km, and 250 km is used as a slightly conservative
#'   round figure. That estimate is drawn from permanent residents
#'   specifically and may not be the right margin for strongly migratory
#'   species. It is a provisional, literature-informed default that should
#'   be revisited as more species run through the pipeline.
#'
#'   Huang, Q., J. R. Sauer, and R. O. Dubayah. 2017. Multidirectional
#'   abundance shifts among North American birds and the relative
#'   influence of multifaceted climate factors. Global Change Biology
#'   23:3610-3622. \doi{10.1111/gcb.13683}
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
#' Three distinct geometric objects are involved, and they are not
#' interchangeable:
#' \itemize{
#'   \item The raw GAP polygon, the true irregular shape read from the
#'     shapefile. Range-based summary statistics are computed against
#'     this polygon, loaded directly from the shapefile by the functions
#'     that need it; this function does not alter it.
#'   \item The buffered GAP polygon, the raw polygon buffered outward by
#'     \code{buffer_km} in EPSG:5070 (USA Contiguous Albers Equal Area
#'     Conic), the equal-area CRS already used for CONUS vector work
#'     elsewhere in the framework. It is written to
#'     \code{runs/ALPHA_CODE/_occs/range_buffered.gpkg} in that same CRS,
#'     so that boundary/buffer-ring statistics can use the polygon itself
#'     rather than only its bounding box.
#'   \item The extent bounding box, the rectangle enclosing the buffered
#'     polygon, written to \code{runs/ALPHA_CODE/_occs/extent.txt} in
#'     WGS84 (EPSG:4326). This is the practical computational domain for
#'     rasters, environmental-variable extraction, and model
#'     construction.
#' }
#'
#' Because the rectangle only has to contain the buffered polygon, its
#' margin exceeds \code{buffer_km} at points away from the polygon's
#' widest extent. That is expected.
#'
#' @return Invisibly returns a named list with numeric elements
#'   \code{ul_lon}, \code{ul_lat}, \code{lr_lon}, and \code{lr_lat}
#'   (upper-left and lower-right corners in decimal degrees), plus
#'   \code{buffered_polygon}, the path to the saved buffered polygon.
#'
#' @seealso
#'   \code{\link{find_occurrence_extent}},
#'   \code{\link{set_extent}},
#'   \code{\link{get_merra_variables}}
#'
#' @importFrom sf st_bbox st_buffer st_make_valid st_read st_sf st_transform
#' @importFrom sf st_union st_write
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' res <- find_range_extent("CASP")
#'
#' find_range_extent(
#'   alpha_code = "CASP",
#'   buffer_km = 300
#' )
#' }
#'
#' @export
find_range_extent <- function(alpha_code, buffer_km = 250, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(buffer_km) || length(buffer_km) != 1L || !is.finite(buffer_km)) {
    stop("`buffer_km` must be a finite numeric scalar.", call. = FALSE)
  }
  if (buffer_km < 0) {
    stop("`buffer_km` must be non-negative.", call. = FALSE)
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
  buffer_fp <- file.path(occs_dir, "range_buffered.gpkg")
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

  ## ---- read and buffer the range polygon ------------------------------------
  # Buffering happens in an equal-area projected CRS, not in lon/lat degrees,
  # so `buffer_km` is a true ground distance in every direction. A degree pad
  # varies with latitude and corresponds to no fixed real-world distance.
  buffer_crs <- "EPSG:5070"

  sfobj    <- sf::st_read(shp_path, quiet = TRUE)
  raw_geom <- sf::st_union(sf::st_make_valid(sfobj))

  .catln("[find_range_extent] Buffering range polygon by ",
         format(buffer_km, trim = TRUE), " km in ", buffer_crs, ".")

  raw_eq      <- sf::st_transform(raw_geom, buffer_crs)
  buffered_eq <- sf::st_buffer(raw_eq, dist = buffer_km * 1000)

  ## ---- persist the buffered polygon -----------------------------------------
  # Saved rather than discarded once its bounding box is known: the
  # buffer-ring boundary statistics need this polygon itself, since the
  # ring is the difference between it and the raw GAP polygon.
  if (file.exists(buffer_fp)) unlink(buffer_fp)
  sf::st_write(sf::st_sf(geometry = buffered_eq), buffer_fp, quiet = TRUE)
  .catln("[find_range_extent] Wrote buffered polygon: ", buffer_fp)

  ## ---- bounding box of the buffered polygon ---------------------------------
  bb <- sf::st_bbox(sf::st_transform(buffered_eq, 4326L))

  ul_lon <- as.numeric(bb[["xmin"]]); lr_lon <- as.numeric(bb[["xmax"]])
  lr_lat <- as.numeric(bb[["ymin"]]); ul_lat <- as.numeric(bb[["ymax"]])

  if (!all(is.finite(c(ul_lon, lr_lon, lr_lat, ul_lat)))) {
    stop("Non-finite bounding box coordinates from buffered polygon.", call. = FALSE)
  }
  if (lr_lon <= ul_lon || ul_lat <= lr_lat) {
    stop("Invalid bounding box: non-positive width or height.", call. = FALSE)
  }

  .catln(sprintf("[find_range_extent] Buffered bbox: UL=(%.6f, %.6f) LR=(%.6f, %.6f)",
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
    paste0("# buffer_km: ", format(buffer_km, trim = TRUE),
           " (fixed real-world buffer applied to source polygon;",
           " box below is the buffered polygon's bounding box)"),
    paste0("# Buffer CRS: ", buffer_crs, " (equal-area; see methods)"),
    paste0("# Buffered polygon: ", basename(buffer_fp)),
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
    sprintf("Buffer file:  %s", buffer_fp),
    sprintf("Buffer:       %s km in %s", format(buffer_km, trim = TRUE), buffer_crs),
    sprintf("UL (lon,lat): (%.6f, %.6f)", ul_lon, ul_lat),
    sprintf("LR (lon,lat): (%.6f, %.6f)", lr_lon, lr_lat)
  ))

  invisible(list(
    ul_lon           = ul_lon,
    ul_lat           = ul_lat,
    lr_lon           = lr_lon,
    lr_lat           = lr_lat,
    buffered_polygon = buffer_fp
  ))
}
