#' Crop and export MERRA-2 predictor rasters for a species run
#'
#' For each 5-year bin from 1980 to 2020, locates MERRA-2 raster files from
#' the \code{m2} and \code{mc} source directories, crops them to the species'
#' spatial extent, and writes the results to the run's \code{_vars/} directory.
#'
#' @details
#' \strong{Input locations}
#'
#' Rasters (\code{.tif}, \code{.tiff}, or \code{.asc}) are read from:
#' \itemize{
#'   \item \code{PROJECT_DIR/data/merra/m2/YEAR/}
#'   \item \code{PROJECT_DIR/data/merra/mc/YEAR/}
#' }
#' When multiple files share the same basename, GeoTIFF is preferred over ASC.
#' Use \code{m2_vars} and \code{mc_vars} to restrict processing to specific
#' variable names (matched by basename, case-insensitive).
#'
#' \strong{Spatial extent}
#'
#' The bounding box is read from \code{RUN_DIR/_occs/extent.txt}, which must
#' be created first by \code{\link{find_occurrence_extent}},
#' \code{\link{find_range_extent}}, or \code{\link{set_extent}}.
#'
#' \strong{CRS handling}
#'
#' Geographic rasters (lon/lat CRS or EPSG:4326) with 0-360 longitude are
#' rotated to -180-180 before cropping. Projected rasters are cropped by
#' projecting the WGS84 bounding box polygon into the raster CRS. Rasters that
#' produce no valid cells after cropping are skipped.
#'
#' \strong{Output}
#'
#' Cropped rasters are written to \code{RUN_DIR/_vars/YEAR/}. For
#' \code{.asc} output, GDAL sidecar files (\code{.aux.xml}, \code{.xml},
#' \code{.prj}) are removed. Missing CRS on output is assigned WGS84.
#'
#' @param alpha_code Character scalar. Four-letter species code.
#' @param m2_vars Character vector or \code{NULL}. Variable basenames to
#'   include from the \code{m2} directory. \code{NULL} includes all.
#' @param mc_vars Character vector or \code{NULL}. Variable basenames to
#'   include from the \code{mc} directory. \code{NULL} includes all.
#' @param file_type Character. Output format(s): \code{".tif"}, \code{".asc"},
#'   or both. Default is both.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a data frame with one row per candidate raster
#'   and columns \code{year}, \code{family}, \code{basename}, \code{in_file},
#'   \code{tif_written}, \code{asc_written}, and \code{notes} (reason for
#'   skipping, or empty string on success).
#'
#' @seealso \code{\link{find_occurrence_extent}}, \code{\link{find_range_extent}},
#'   \code{\link{set_extent}}
#'
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' get_merra_variables("CASP")
#' get_merra_variables("CASP", m2_vars = c("bio1", "bio12"),
#'                     project_dir = "/projects/rENM")
#' }
#'
#' @export
get_merra_variables <- function(alpha_code,
                                m2_vars    = NULL,
                                mc_vars    = NULL,
                                file_type  = c(".tif", ".asc"),
                                project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  file_type <- unique(match.arg(file_type, choices = c(".tif", ".asc"), several.ok = TRUE))

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. Install with: install.packages('terra')", call. = FALSE)
  }

  .t              <- getNamespace("terra")
  .basename_noext <- function(x) tools::file_path_sans_ext(basename(x))

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir      <- .expand(file.path(project_root, "runs", alpha_code))
  vars_dir     <- file.path(run_dir, "_vars")
  m2_root      <- .expand(file.path(project_root, "data", "merra", "m2"))
  mc_root      <- .expand(file.path(project_root, "data", "merra", "mc"))
  log_fp       <- file.path(run_dir, "_log.txt")
  .mkdir(vars_dir)

  bins <- seq(1980L, 2020L, by = 5L)

  ## ---- parse extent.txt -----------------------------------------------------
  extent_fp <- file.path(run_dir, "_occs", "extent.txt")
  if (!file.exists(extent_fp)) {
    stop(
      "extent.txt not found at: ", extent_fp, "\n",
      "Run find_occurrence_extent(), find_range_extent(), or set_extent() first.",
      call. = FALSE
    )
  }

  .parse_pair <- function(s) {
    m <- regexec("\\(([-0-9.]+)\\s*,\\s*([-0-9.]+)\\)", s, perl = TRUE)
    v <- regmatches(s, m)
    if (length(v) == 1L && length(v[[1L]]) == 3L) as.numeric(v[[1L]][2:3])
    else c(NA_real_, NA_real_)
  }

  txt    <- readLines(extent_fp, warn = FALSE)
  ul     <- .parse_pair(grep("Upper-left",  txt, ignore.case = TRUE, value = TRUE)[1L])
  lr     <- .parse_pair(grep("Lower-right", txt, ignore.case = TRUE, value = TRUE)[1L])
  if (any(is.na(c(ul, lr)))) stop("Failed to parse coordinates from extent.txt.", call. = FALSE)

  xmin  <- min(ul[1L], lr[1L]); xmax <- max(ul[1L], lr[1L])
  ymin  <- min(ul[2L], lr[2L]); ymax <- max(ul[2L], lr[2L])
  wgs84 <- "EPSG:4326"
  bb    <- list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

  ## ---- build WGS84 bbox polygon ---------------------------------------------
  xs         <- c(bb$xmin, bb$xmax, bb$xmax, bb$xmin, bb$xmin)
  ys         <- c(bb$ymin, bb$ymin, bb$ymax, bb$ymax, bb$ymin)
  bbox_wgs84 <- .t$vect(cbind(xs, ys), type = "polygons", atts = NULL, crs = wgs84)

  ## ---- candidate discovery --------------------------------------------------
  .list_family <- function(root, family, year, vars = NULL) {
    year_dir <- file.path(root, as.character(year))
    if (!dir.exists(year_dir)) return(data.frame())
    files <- list.files(year_dir, pattern = "\\.(tif|tiff|asc)$",
                        full.names = TRUE, ignore.case = TRUE)
    if (!length(files)) return(data.frame())
    base <- .basename_noext(files)
    df   <- data.frame(family = family, year = year, in_file = files,
                       basename = base, stringsAsFactors = FALSE)
    if (!is.null(vars)) df <- df[.norm_names(df$basename) %in% .norm_names(vars), , drop = FALSE]
    if (!nrow(df)) return(df)
    # prefer tif over asc when basenames collide
    ord <- order(df$basename, grepl("\\.asc$", df$in_file, ignore.case = TRUE))
    df  <- df[ord, , drop = FALSE]
    df[!duplicated(df$basename), , drop = FALSE]
  }

  all_rows <- do.call(rbind, lapply(bins, function(y) {
    rbind(.list_family(m2_root, "m2", y, m2_vars),
          .list_family(mc_root, "mc", y, mc_vars))
  }))

  if (!nrow(all_rows)) {
    stop("No candidate MERRA-2 rasters found under m2/mc for the requested bins.", call. = FALSE)
  }

  ## ---- raster helpers -------------------------------------------------------
  .read_rast <- function(fp) {
    r <- try(.t$rast(fp), silent = TRUE)
    if (inherits(r, "try-error")) { .catln("  ! Failed to read: ", basename(fp)); return(NULL) }
    r
  }

  .is_geographic <- function(r) {
    info <- .t$crs(r, proj = TRUE)
    if (is.na(info) || !nzchar(trimws(info))) return(FALSE)
    grepl("longlat", tolower(info), fixed = TRUE) || grepl("4326", info)
  }

  .crop_geographic <- function(r, poly) {
    ex <- .t$ext(r)
    if (ex$xmin >= 0 && ex$xmax <= 360) {
      .catln("  - Rotating 0-360 longitude to -180/180.")
      r <- .t$rotate(r)
    }
    .t$crop(r, poly, snap = "near")
  }

  .crop_projected <- function(r, poly) {
    r_crs <- terra::crs(r, proj = TRUE)
    if (is.na(r_crs) || !nzchar(trimws(r_crs))) {
      .catln("  - Skipping (no CRS)."); return(NULL)
    }
    bb_proj <- try(.t$project(poly, r_crs), silent = TRUE)
    if (inherits(bb_proj, "try-error")) {
      .catln("  - Skipping (bbox projection failed)."); return(NULL)
    }
    .t$crop(r, bb_proj, snap = "near")
  }

  .write_tif <- function(x, fp) {
    r_crs <- terra::crs(x, proj = TRUE)
    if (is.na(r_crs) || !nzchar(trimws(r_crs))) terra::crs(x) <- wgs84
    ok <- try(.t$writeRaster(x, fp, overwrite = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) { .catln("  ! Write error (tif): ", fp); return(FALSE) }
    .catln("  * Wrote: ", fp); TRUE
  }

  .write_asc <- function(x, fp) {
    ok1 <- try(terra::writeRaster(x, fp, filetype = "AAIGrid",
                                  gdal = c("DECIMAL_PRECISION=6", "FORCE_CELLSIZE=TRUE"),
                                  NAflag = -9999, overwrite = TRUE), silent = TRUE)
    if (inherits(ok1, "try-error")) {
      ok2 <- try(terra::writeRaster(x, fp, NAflag = -9999, overwrite = TRUE), silent = TRUE)
      if (inherits(ok2, "try-error")) {
        msg <- paste(stats::na.omit(c(as.character(ok1), as.character(ok2))), collapse = " | ")
        .catln("  ! Write error (asc): ", fp, " - ", msg)
        return(FALSE)
      }
    }
    .catln("  * Wrote: ", fp)
    asc_root <- tools::file_path_sans_ext(fp)
    for (sc in c(paste0(fp, ".aux.xml"), paste0(asc_root, ".xml"), paste0(asc_root, ".prj"))) {
      if (file.exists(sc)) {
        ok_rm <- try(file.remove(sc), silent = TRUE)
        if (!inherits(ok_rm, "try-error") && ok_rm) .catln("    - Removed sidecar: ", basename(sc))
      }
    }
    TRUE
  }

  ## ---- main loop ------------------------------------------------------------
  .catln(.sep_line())
  .catln("get_merra_variables: alpha_code=", alpha_code)
  .catln(sprintf("Extent (WGS84): xmin=%.6f xmax=%.6f ymin=%.6f ymax=%.6f",
                 bb$xmin, bb$xmax, bb$ymin, bb$ymax))
  .catln("M2 root:     ", m2_root)
  .catln("MC root:     ", mc_root)
  .catln("Output root: ", vars_dir)
  .catln(.sep_line())

  out_rows <- list()

  for (y in bins) {
    this <- all_rows[all_rows$year == y, , drop = FALSE]
    if (!nrow(this)) next

    out_dir <- file.path(vars_dir, as.character(y))
    .mkdir(out_dir)
    .catln(.sep_line())
    .catln("Year bin: ", y, " (", nrow(this), " candidate raster(s))")

    wrote_tif <- 0L; wrote_asc <- 0L

    for (i in seq_len(nrow(this))) {
      fam   <- this$family[i]
      in_fp <- this$in_file[i]
      base  <- this$basename[i]

      .catln(" - ", fam, "/", base, " (", basename(in_fp), ")")
      r <- .read_rast(in_fp)

      if (is.null(r)) {
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year = y, family = fam, basename = base, in_file = in_fp,
          tif_written = FALSE, asc_written = FALSE, notes = "read_failed",
          stringsAsFactors = FALSE)
        next
      }

      r_crop <- if (.is_geographic(r)) {
        try(.crop_geographic(r, bbox_wgs84), silent = TRUE)
      } else {
        try(.crop_projected(r, bbox_wgs84), silent = TRUE)
      }

      if (inherits(r_crop, "try-error") || is.null(r_crop)) {
        .catln("  - Skipping (crop failed).")
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year = y, family = fam, basename = base, in_file = in_fp,
          tif_written = FALSE, asc_written = FALSE, notes = "crop_failed",
          stringsAsFactors = FALSE)
        next
      }

      vals <- try(.t$values(r_crop), silent = TRUE)
      if (inherits(vals, "try-error") || all(is.na(vals))) {
        .catln("  - Skipping (no valid cells after crop).")
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year = y, family = fam, basename = base, in_file = in_fp,
          tif_written = FALSE, asc_written = FALSE, notes = "no_cells_after_crop",
          stringsAsFactors = FALSE)
        next
      }

      base_out      <- file.path(out_dir, base)
      wrote_this_tif <- FALSE
      wrote_this_asc <- FALSE

      if (".tif" %in% file_type && .write_tif(r_crop, paste0(base_out, ".tif"))) {
        wrote_tif <- wrote_tif + 1L; wrote_this_tif <- TRUE
      }
      if (".asc" %in% file_type && .write_asc(r_crop, paste0(base_out, ".asc"))) {
        wrote_asc <- wrote_asc + 1L; wrote_this_asc <- TRUE
      }

      out_rows[[length(out_rows) + 1L]] <- data.frame(
        year = y, family = fam, basename = base, in_file = in_fp,
        tif_written = wrote_this_tif, asc_written = wrote_this_asc,
        notes = if (!wrote_this_tif && !wrote_this_asc) "write_failed" else "",
        stringsAsFactors = FALSE)
    }

    .catln(sprintf("  - Wrote %d .tif and %d .asc file(s) to %s",
                   wrote_tif, wrote_asc, out_dir))
  }

  out_df <- if (length(out_rows)) do.call(rbind, out_rows) else data.frame(
    year = integer(0), family = character(0), basename = character(0),
    in_file = character(0), tif_written = logical(0), asc_written = logical(0),
    notes = character(0), stringsAsFactors = FALSE)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (get_merra_variables)", c(
    sprintf("Alpha code:     %s", alpha_code),
    sprintf("Extent (WGS84): xmin=%.6f xmax=%.6f ymin=%.6f ymax=%.6f",
            bb$xmin, bb$xmax, bb$ymin, bb$ymax),
    sprintf("Input m2 root:  %s", m2_root),
    sprintf("Input mc root:  %s", mc_root),
    sprintf("Output root:    %s", vars_dir),
    sprintf("Rasters processed: %d", nrow(out_df)),
    sprintf("Rasters written:   %d", sum(out_df$tif_written | out_df$asc_written))
  ))

  .catln(.sep_line())
  .catln("Done. Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(out_df)
}
