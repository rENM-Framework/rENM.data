#' Find the spatial extent of a species' USGS GAP range
#'
#' Computes a bounding box from a species GAP range polygon and writes
#' the result to extent.txt for use in rENM preprocessing workflows.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' The species table is read from:
#'
#' \code{<project_dir>/data/}
#' \code{_species.csv}
#'
#' with fallback to:
#'
#' \code{<project_dir>/data/}
#' \code{species.csv}
#'
#' The row where \code{ALPHA.CODE == alpha_code} is located and the
#' corresponding \code{GAP.RANGE} value is extracted.
#'
#' The GAP range shapefile is then loaded from:
#'
#' \code{<project_dir>/data/shapefiles/}
#' \code{<gap_range>/<gap_range>.shp}
#'
#' \strong{Methods}:
#' The shapefile geometry is validated and its coordinate reference
#' system (CRS) is checked. If necessary, the geometry is transformed to
#' WGS84 (EPSG:4326).
#'
#' A bounding box is computed in longitude and latitude degrees and then
#' adjusted using a symmetric padding factor defined by \code{pad_pct}.
#'
#' The padding is applied relative to the bounding box center. For
#' example:
#' \itemize{
#'   \item pad_pct = 2 expands width and height by 2 percent
#'   \item pad_pct = -2 contracts width and height by 2 percent
#' }
#'
#' Values less than or equal to -100 are invalid.
#'
#' \strong{Outputs}:
#' The resulting bounding box is written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/extent.txt}
#'
#' Any existing extent.txt file is backed up prior to writing.
#'
#' The output file contains:
#' \itemize{
#'   \item Upper-left:  (<lon>, <lat>)
#'   \item Lower-right: (<lon>, <lat>)
#' }
#'
#' Coordinates are written in decimal degrees with six decimal places of
#' precision.
#'
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' using the section title:
#'
#' \code{Processing summary (find_range_extent)}
#'
#' \strong{Workflow context}:
#' This function is intended for workflows that rely on GAP range
#' polygons to define spatial extents for downstream analyses. It can be
#' used alongside occurrence-based preprocessing steps such as:
#' \itemize{
#'   \item get_ebird_occurrences()
#'   \item remove_duplicate_occurrences()
#'   \item thin_occurrences() or thin_occurrences2()
#'   \item tidy_occurrences()
#' }
#'
#' \strong{Data requirements}:
#' The species CSV must contain:
#' \itemize{
#'   \item ALPHA.CODE: species identifier
#'   \item GAP.RANGE: shapefile directory name
#' }
#'
#' The shapefile must exist at the expected path and contain valid
#' geometry.
#'
#' @param alpha_code Character. Species alpha code matching ALPHA.CODE in
#'   the species CSV file.
#' @param pad_pct Numeric. Percent padding applied to the bounding box.
#'   Positive values expand and negative values contract the box. Must be
#'   greater than -100. Default = 2.
#'
#' @return List (invisibly returned) with elements:
#' \itemize{
#'   \item ul_lon: Numeric upper-left longitude
#'   \item ul_lat: Numeric upper-left latitude
#'   \item lr_lon: Numeric lower-right longitude
#'   \item lr_lat: Numeric lower-right latitude
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes extent.txt to the _occs directory
#'   \item Backs up any existing extent.txt file
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom sf st_read st_crs st_is_longlat st_transform st_make_valid st_bbox
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#'   res <- find_range_extent("CASP", pad_pct = 2)
#'   res
#'
#'   res_tight <- find_range_extent("CASP", pad_pct = -1)
#' }
#'
#' @export
find_range_extent <- function(alpha_code, pad_pct = 2) {
  # ----------------------------- helpers --------------------------------------
  .ts  <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .f6  <- function(x) sprintf("%.6f", as.numeric(x))
  .sep <- function(n = 72L) paste(rep.int("-", n), collapse = "")
  .step <- function(...) {
    cat(paste0("[find_range_extent] ", paste0(..., collapse = ""), "\n"))
    flush.console()
  }

  .expand <- function(p) {
    normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  }

  .append_log <- function(path, lines) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    con <- file(path, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con = con)
  }

  .backup_file <- function(path) {
    if (!file.exists(path)) return(NULL)
    ts <- gsub("[^0-9]", "", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    bname <- sprintf("extent_backup_%s.txt", ts)
    backup_path <- file.path(dirname(path), bname)
    file.rename(path, backup_path)
    backup_path
  }

  # ----------------------------- validation -----------------------------------
  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nzchar(alpha_code))
  stopifnot(is.numeric(pad_pct), length(pad_pct) == 1, is.finite(pad_pct))

  if (pad_pct <= -100) {
    stop("pad_pct must be greater than -100 (percent).")
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required but not installed.", call. = FALSE)
  }

  # ----------------------------- project directory ----------------------------
  project_dir <- rENM_project_dir()

  # ----------------------------- paths ----------------------------------------
  species_csv_candidates <- c(
    file.path(project_dir, "data", "_species.csv"),
    file.path(project_dir, "data", "species.csv")
  )

  species_csv_candidates <- vapply(
    species_csv_candidates,
    .expand,
    character(1)
  )

  species_csv <- species_csv_candidates[file.exists(species_csv_candidates)]

  if (!length(species_csv)) {
    stop(
      "Could not find species CSV. Tried: ",
      paste(species_csv_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  species_csv <- species_csv[1]
  .step("Using species CSV: ", species_csv)

  # Run and output paths
  run_dir  <- file.path(project_dir, "runs", alpha_code)
  occs_dir <- file.path(run_dir, "_occs")
  log_file <- file.path(run_dir, "_log.txt")
  out_file <- file.path(occs_dir, "extent.txt")
  dir.create(occs_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------- read species table ---------------------------
  sp <- utils::read.csv(species_csv, stringsAsFactors = FALSE, check.names = FALSE)

  if (!all(c("ALPHA.CODE", "GAP.RANGE") %in% names(sp))) {
    stop(
      "Species CSV must contain 'ALPHA.CODE' and 'GAP.RANGE' columns.",
      call. = FALSE
    )
  }

  row_idx <- which(sp$ALPHA.CODE == alpha_code)

  if (!length(row_idx)) {
    stop(
      "No row found in species CSV with ALPHA.CODE == '", alpha_code, "'.",
      call. = FALSE
    )
  }

  gap_range <- sp$GAP.RANGE[row_idx[1]]
  .step(sprintf("Resolved GAP.RANGE='%s'", gap_range))

  # ----------------------------- shapefile path -------------------------------
  shp_dir  <- file.path(project_dir, "data", "shapefiles", gap_range)
  shp_path <- file.path(shp_dir, paste0(gap_range, ".shp"))

  if (!file.exists(shp_path)) {
    stop("Shapefile not found: ", shp_path, call. = FALSE)
  }

  .step("Shapefile path: ", shp_path)

  # ----------------------------- rest of function unchanged --------------------
  .step("Reading shapefile via sf::st_read()...")
  sfobj <- sf::st_read(shp_path, quiet = TRUE)

  crs_obj <- sf::st_crs(sfobj)
  if (is.null(crs_obj) || is.na(crs_obj$epsg)) {
    .step("CRS is missing or has no EPSG code; proceeding with provided geometry.")
  }

  if (sf::st_is_longlat(sfobj)) {
    geom_ll <- sfobj
    .step("CRS appears to be longitude/latitude. Using as-is.")
  } else {
    .step("CRS is projected or non-longlat. Transforming to EPSG:4326 (WGS84).")
    geom_ll <- sf::st_transform(sfobj, 4326)
  }

  rng <- sf::st_make_valid(geom_ll)

  bb <- sf::st_bbox(rng)

  xmin <- as.numeric(bb[["xmin"]])
  xmax <- as.numeric(bb[["xmax"]])
  ymin <- as.numeric(bb[["ymin"]])
  ymax <- as.numeric(bb[["ymax"]])

  if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) {
    stop("Non-finite bounding box coordinates from shapefile.", call. = FALSE)
  }

  .step(sprintf(
    "Original bbox (lon/lat): xmin=%.6f, xmax=%.6f, ymin=%.6f, ymax=%.6f",
    xmin, xmax, ymin, ymax
  ))

  width  <- xmax - xmin
  height <- ymax - ymin

  if (width <= 0 || height <= 0) {
    stop("Invalid bbox: non-positive width or height.", call. = FALSE)
  }

  factor <- 1 + pad_pct / 100
  new_width  <- width  * factor
  new_height <- height * factor

  cx <- (xmin + xmax) / 2
  cy <- (ymin + ymax) / 2

  xmin2 <- cx - new_width  / 2
  xmax2 <- cx + new_width  / 2
  ymin2 <- cy - new_height / 2
  ymax2 <- cy + new_height / 2

  ul_lon <- xmin2
  ul_lat <- ymax2
  lr_lon <- xmax2
  lr_lat <- ymin2

  .step(sprintf(
    "Padded bbox (lon/lat): xmin=%.6f, xmax=%.6f, ymin=%.6f, ymax=%.6f",
    xmin2, xmax2, ymin2, ymax2
  ))

  backup_path <- .backup_file(out_file)
  if (!is.null(backup_path)) {
    .step("Backed up existing extent.txt to: ", backup_path)
  }

  tz_stamp <- .ts()

  header_lines <- c(
    "# extent.txt",
    "# Computed by find_range_extent()",
    paste0("# Timestamp: ", tz_stamp),
    paste0("# Source: ", basename(shp_path), " via find_range_extent()"),
    "# Points used: n/a",
    paste0("# pad_pct: ", formatC(pad_pct, format = "f", digits = 1), "% (symmetric percent padding)"),
    "# Coordinate order: (lon, lat)"
  )

  extent_lines <- c(
    sprintf("Upper-left:  (%s, %s)", .f6(ul_lon), .f6(ul_lat)),
    sprintf("Lower-right: (%s, %s)", .f6(lr_lon), .f6(lr_lat))
  )

  writeLines(c(header_lines, extent_lines), con = out_file)
  .step("Wrote extent.txt.")

  ul_fmt <- sprintf("(%s, %s)", .f6(ul_lon), .f6(ul_lat))
  lr_fmt <- sprintf("(%s, %s)", .f6(lr_lon), .f6(lr_lat))

  summary_block <- c(
    "",
    .sep(),
    "Processing summary (find_range_extent)",
    sprintf("%-18s %s", "Timestamp:", tz_stamp),
    sprintf("%-18s %s", "Alpha code:", alpha_code),
    sprintf("%-18s %s", "Species CSV:", species_csv),
    sprintf("%-18s %s", "GAP.RANGE:", gap_range),
    sprintf("%-18s %s", "Shapefile:", shp_path),
    sprintf("%-18s %s", "Output file:", out_file),
    sprintf("%-18s %s", "Pad pct:", sprintf("%.1f%%", pad_pct)),
    sprintf("%-18s %s", "UL (lon,lat):", ul_fmt),
    sprintf("%-18s %s", "LR (lon,lat):", lr_fmt)
  )

  .append_log(log_file, summary_block)
  cat(paste0(paste(summary_block, collapse = "\n"), "\n"))

  invisible(list(
    ul_lon = ul_lon,
    ul_lat = ul_lat,
    lr_lon = lr_lon,
    lr_lat = lr_lat
  ))
}
