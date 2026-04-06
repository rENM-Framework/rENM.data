#' Assemble a run-tailored set of MERRA-2 variables for further processing
#'
#' Crops and exports MERRA-2 variables into 5-year bins for a species,
#' writing run-specific predictor rasters for downstream rENM modeling.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' MERRA-2 rasters are read from:
#'
#' \code{<project_dir>/data/merra/m2/<year>/}
#' \code{<project_dir>/data/merra/mc/<year>/}
#'
#' where <year> corresponds to 5-year bins (1980, 1985, ..., 2020).
#'
#' Files may be in .tif, .tiff, or .asc format. If multiple files share
#' the same basename, duplicates are resolved with preference for .tif.
#'
#' The spatial extent is read from:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/extent.txt}
#'
#' which must contain "Upper-left" and "Lower-right" coordinates in
#' (lon, lat) order.
#'
#' \strong{Methods}:
#' For each 5-year bin, the function:
#' \itemize{
#'   \item Identifies candidate rasters from m2 and mc directories
#'   \item Filters by basename using m2_vars and mc_vars if provided
#'   \item De-duplicates files by basename with preference for GeoTIFF
#'   \item Crops rasters to the specified bounding box
#' }
#'
#' Geographic rasters:
#' \itemize{
#'   \item Detects longitude range (0-360) and rotates to (-180, 180)
#'   \item Handles antimeridian crossing during cropping
#' }
#'
#' Projected rasters:
#' \itemize{
#'   \item Builds a WGS84 bounding box polygon
#'   \item Projects the polygon into the raster CRS
#'   \item Crops using the projected geometry
#' }
#'
#' Only rasters with valid cells after cropping are written.
#'
#' \strong{Outputs}:
#' Cropped rasters are written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_vars/<year>/}
#'
#' Output formats are controlled by \code{file_type}:
#' \itemize{
#'   \item .tif: GeoTIFF with CRS preserved or assigned (WGS84 if missing)
#'   \item .asc: ASCII grid; sidecar files (.aux.xml, .xml, .prj) are removed
#' }
#'
#' \strong{Log file}:
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' including:
#' \itemize{
#'   \item alpha code
#'   \item extent used (xmin, xmax, ymin, ymax)
#'   \item input m2 and mc roots
#'   \item output directory
#'   \item per-year processing results
#'   \item skipped files and reasons
#' }
#'
#' \strong{Data requirements}:
#' The extent.txt file must exist and be correctly formatted. Input
#' raster directories must contain at least one valid raster file.
#'
#' @param alpha_code Character. Four-letter species code used to locate
#'   the run directory and extent file.
#' @param m2_vars Character, NULL. Optional vector of variable basenames
#'   to include from the m2 directory.
#' @param mc_vars Character, NULL. Optional vector of variable basenames
#'   to include from the mc directory.
#' @param file_type Character. Output formats to write. Allowed values
#'   are ".tif" and ".asc".
#'
#' @return Invisibly returns a Data frame summarizing processing results:
#' \itemize{
#'   \item year: Integer bin year
#'   \item family: Character indicating m2 or mc source
#'   \item basename: Character variable name
#'   \item in_file: Character input file path
#'   \item tif_written: Logical indicating GeoTIFF output success
#'   \item asc_written: Logical indicating ASCII output success
#'   \item notes: Character processing notes or failure reason
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes cropped raster files to disk
#'   \item Creates output directories as needed
#'   \item Removes ASCII sidecar files when applicable
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#'   get_merra_variables("CASP")
#' }
#'
#' @export
get_merra_variables <- function(alpha_code,
                                m2_vars   = NULL,
                                mc_vars   = NULL,
                                file_type = c(".tif", ".asc")) {

  ## ---- basic validation ------------------------------------------------------
  stopifnot(is.character(alpha_code), length(alpha_code) == 1L, nzchar(alpha_code))
  file_type <- unique(match.arg(file_type, choices = c(".tif", ".asc"), several.ok = TRUE))

  ## ---- namespace checks ------------------------------------------------------
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Please install.packages('terra').")
  }
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("Please install.packages('tools').")
  }
  .t <- getNamespace("terra")

  ## ---- helpers ---------------------------------------------------------------
  .expand     <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  .mkdir      <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  .now        <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep        <- function(n = 72L) paste(rep.int("-", n), collapse = "")
  .catln      <- function(...) { cat(paste0(..., "\n")); flush.console() }
  .norm       <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
  .basename_noext <- function(x) tools::file_path_sans_ext(basename(x))

  .err_txt <- function(e) {
    if (inherits(e, "try-error")) as.character(e) else NA_character_
  }

  ## ---- paths (UPDATED ONLY) --------------------------------------------------
  bins <- seq(1980L, 2020L, by = 5L)

  project_dir <- rENM_project_dir()

  run_dir  <- .expand(file.path(project_dir, "runs", alpha_code))
  vars_dir <- file.path(run_dir, "_vars")
  .mkdir(vars_dir)

  m2_root <- .expand(file.path(project_dir, "data", "merra", "m2"))
  mc_root <- .expand(file.path(project_dir, "data", "merra", "mc"))

  ## ---- logging helper --------------------------------------------------------
  log_file <- file.path(run_dir, "_log.txt")
  .append_log <- function(title, lines) {
    hdr <- c("", .sep(), title, paste0("Timestamp: ", .now()))
    cat(paste0(c(hdr, lines, ""), collapse = "\n"), file = log_file, append = TRUE)
    invisible(log_file)
  }

  ## ---- parse extent.txt ------------------------------------------------------
  .parse_pair <- function(s) {
    m <- regexec("\\(([-0-9.]+)\\s*,\\s*([-0-9.]+)\\)", s, perl = TRUE)
    v <- regmatches(s, m)
    if (length(v) == 1L && length(v[[1]]) == 3L) {
      as.numeric(v[[1]][2:3])
    } else {
      c(NA_real_, NA_real_)
    }
  }
  .parse_extent <- function(fp) {
    txt <- readLines(fp, warn = FALSE)
    ul_line <- grep("Upper-left", txt, ignore.case = TRUE, value = TRUE)
    lr_line <- grep("Lower-right", txt, ignore.case = TRUE, value = TRUE)
    if (!length(ul_line) || !length(lr_line)) {
      stop("extent.txt must contain 'Upper-left' and 'Lower-right' lines.")
    }
    ul <- .parse_pair(ul_line[1L])
    lr <- .parse_pair(lr_line[1L])
    if (any(is.na(c(ul, lr)))) stop("Failed to parse coordinates from extent.txt.")
    list(ul = ul, lr = lr)
  }

  extent_fp <- file.path(run_dir, "_occs", "extent.txt")
  if (!file.exists(extent_fp)) {
    stop("extent.txt not found at: ", extent_fp, "\n",
         "Run find_occurrence_extent() before get_merra_variables().")
  }
  ex <- .parse_extent(extent_fp)

  xmin <- min(ex$ul[1L], ex$lr[1L])
  xmax <- max(ex$ul[1L], ex$lr[1L])
  ymin <- min(ex$lr[2L], ex$ul[2L])
  ymax <- max(ex$lr[2L], ex$ul[2L])

  wgs84 <- "EPSG:4326"
  bb <- list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

  ## ---- build WGS84 bbox polygon ----------------------------------------------
  .bbox_polygon_wgs84 <- function(bb) {
    xs <- c(bb$xmin, bb$xmax, bb$xmax, bb$xmin, bb$xmin)
    ys <- c(bb$ymin, bb$ymin, bb$ymax, bb$ymax, bb$ymin)
    .t$vect(cbind(xs, ys), type = "polygons", atts = NULL, crs = wgs84)
  }
  bbox_wgs84 <- .bbox_polygon_wgs84(bb)

  ## ---- candidate discovery ---------------------------------------------------
  .list_family <- function(root, family, year, vars = NULL) {
    year_dir <- file.path(root, as.character(year))
    if (!dir.exists(year_dir)) return(data.frame())
    files <- list.files(year_dir, pattern = "\\.(tif|tiff|asc)$",
                        full.names = TRUE, ignore.case = TRUE)
    if (!length(files)) return(data.frame())

    base <- .basename_noext(files)
    df <- data.frame(
      family   = family,
      year     = year,
      in_file  = files,
      basename = base,
      stringsAsFactors = FALSE
    )
    if (!is.null(vars)) {
      keep <- .norm(df$basename) %in% .norm(vars)
      df <- df[keep, , drop = FALSE]
    }
    if (!nrow(df)) return(df)

    # De-duplicate by basename, preferring tif/tiff over asc
    ord <- order(df$basename,
                 grepl("\\.asc$", df$in_file, ignore.case = TRUE),
                 decreasing = FALSE)
    df <- df[ord, , drop = FALSE]
    df <- df[!duplicated(df$basename), , drop = FALSE]
    df
  }

  all_rows <- do.call(rbind, lapply(bins, function(y) {
    rbind(
      .list_family(m2_root, "m2", y, m2_vars),
      .list_family(mc_root, "mc", y, mc_vars)
    )
  }))
  if (!nrow(all_rows)) {
    stop("No candidate MERRA-2 rasters found under m2/mc roots for the requested bins.")
  }

  ## ---- raster helpers --------------------------------------------------------
  .read_rast <- function(fp) {
    ok <- try(.t$rast(fp), silent = TRUE)
    if (inherits(ok, "try-error")) {
      .catln("  ! Failed to read raster: ", basename(fp), " - ", .err_txt(ok))
      return(NULL)
    }
    ok
  }

  .is_geographic <- function(r) {
    info <- .t$crs(r, proj = TRUE)
    if (is.na(info) || !nzchar(trimws(info))) return(FALSE)
    grepl("longlat", tolower(info), fixed = TRUE) || grepl("4326", info)
  }

  .maybe_rotate_0_360 <- function(r) {
    ex <- .t$ext(r)
    if (ex$xmin >= 0 && ex$xmax <= 360) {
      .catln("  - Rotating 0-360 longitude raster into -180 to 180.")
      r <- .t$rotate(r)
    }
    r
  }

  .crop_geographic <- function(r, bb_poly) {
    r <- .maybe_rotate_0_360(r)
    .t$crop(r, bb_poly, snap = "near")
  }

  .crop_projected <- function(r, bb_poly) {
    r_crs <- terra::crs(r, proj = TRUE)
    if (is.na(r_crs) || !nzchar(trimws(r_crs))) {
      .catln("  - Skipping (no CRS defined for projected raster).")
      return(NULL)
    }
    bb_proj <- try(.t$project(bb_poly, r_crs), silent = TRUE)
    if (inherits(bb_proj, "try-error")) {
      .catln("  - Skipping (failed to project bbox polygon to raster CRS).")
      return(NULL)
    }
    .t$crop(r, bb_proj, snap = "near")
  }

  .write_tif <- function(x, out_fp) {
    x_crs <- terra::crs(x, proj = TRUE)
    if (is.na(x_crs) || !nzchar(trimws(x_crs))) terra::crs(x) <- wgs84
    ok <- try(.t$writeRaster(x, out_fp, overwrite = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) {
      .catln("  ! Write error (tif): ", out_fp, " - ", .err_txt(ok))
      FALSE
    } else {
      .catln("  * Wrote: ", out_fp)
      TRUE
    }
  }

  # .write_asc <- function(x, out_fp) {
  #   ok1 <- try(.t$writeRaster(
  #     x, out_fp,
  #     filetype = "AAIGrid",
  #     gdal     = c("DECIMAL_PRECISION=6", "FORCE_CELLSIZE=TRUE"),
  #     NAflag   = -9999,
  #     overwrite = TRUE
  #   ), silent = TRUE)
  #
  #   if (!inherits(ok1, "try-error")) {
  #     .catln("  * Wrote: ", out_fp)
  #     return(TRUE)
  #   }
  #
  #   ok2 <- try(.t$writeRaster(x, out_fp, NAflag = -9999, overwrite = TRUE), silent = TRUE)
  #   if (!inherits(ok2, "try-error")) {
  #     .catln("  * Wrote: ", out_fp)
  #     return(TRUE)
  #   }
  #
  #   msg <- paste(stats::na.omit(c(.err_txt(ok1), .err_txt(ok2))), collapse = " | ")
  #   .catln("  ! Write error (asc): ", out_fp, " - ", msg)
  #   FALSE
  # }

  .write_asc <- function(x, out_fp) {
    # Attempt 1: explicit AAIGrid
    ok1 <- try(
      terra::writeRaster(
        x, out_fp,
        filetype = "AAIGrid",
        gdal     = c("DECIMAL_PRECISION=6", "FORCE_CELLSIZE=TRUE"),
        NAflag   = -9999,
        overwrite = TRUE
      ),
      silent = TRUE
    )

    # Attempt 2 (fallback): generic writeRaster()
    if (inherits(ok1, "try-error")) {
      ok2 <- try(
        terra::writeRaster(
          x, out_fp,
          NAflag = -9999,
          overwrite = TRUE
        ),
        silent = TRUE
      )
    } else {
      ok2 <- ok1
    }

    # If both attempts failed, issue error message
    if (inherits(ok2, "try-error")) {
      msg <- paste(stats::na.omit(c(as.character(ok1), as.character(ok2))), collapse = " | ")
      .catln("  ! Write error (asc): ", out_fp, " - ", msg)
      return(FALSE)
    }

    .catln("  * Wrote: ", out_fp)

    # ---- Remove ASCII-grid sidecars -----------------------------------------
    # Terra/GDAL often generate the following beside .asc:
    #   <file>.asc.aux.xml
    #   <file>.xml
    #   <file>.prj
    asc_root <- tools::file_path_sans_ext(out_fp)

    sidecars <- c(
      paste0(out_fp, ".aux.xml"),        # bio1.asc.aux.xml
      paste0(asc_root, ".xml"),          # bio1.xml
      paste0(asc_root, ".prj")           # bio1.prj
    )

    for (sc in sidecars) {
      if (file.exists(sc)) {
        ok_rm <- try(file.remove(sc), silent = TRUE)
        if (!inherits(ok_rm, "try-error") && ok_rm) {
          .catln("    - Removed sidecar: ", basename(sc))
        }
      }
    }
    # -------------------------------------------------------------------------

    TRUE
  }

  ## ---- main loop -------------------------------------------------------------
  out_rows <- list()

  .catln(.sep())
  .catln("get_merra_variables: alpha_code=", alpha_code)
  .catln(sprintf("Extent (WGS84): xmin=%.6f xmax=%.6f ymin=%.6f ymax=%.6f",
                 bb$xmin, bb$xmax, bb$ymin, bb$ymax))
  .catln("M2 root: ", m2_root)
  .catln("MC root: ", mc_root)
  .catln("Output root: ", vars_dir)
  .catln(.sep())

  for (y in bins) {
    this <- all_rows[all_rows$year == y, , drop = FALSE]
    if (!nrow(this)) next

    out_dir <- file.path(vars_dir, as.character(y))
    .mkdir(out_dir)

    .catln(.sep())
    .catln("Year bin: ", y, " (", nrow(this), " candidate rasters)")

    wrote_tif <- 0L
    wrote_asc <- 0L

    for (i in seq_len(nrow(this))) {
      fam   <- this$family[i]
      in_fp <- this$in_file[i]
      base  <- this$basename[i]

      .catln(" - ", fam, "/", base, " (", basename(in_fp), ")")
      r <- .read_rast(in_fp)
      if (is.null(r)) {
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year        = y,
          family      = fam,
          basename    = base,
          in_file     = in_fp,
          tif_written = FALSE,
          asc_written = FALSE,
          notes       = "read_failed",
          stringsAsFactors = FALSE
        )
        next
      }

      # Decide cropping mode
      is_geo <- .is_geographic(r)
      if (is_geo) {
        r_crop <- try(.crop_geographic(r, bbox_wgs84), silent = TRUE)
      } else {
        r_crop <- try(.crop_projected(r, bbox_wgs84), silent = TRUE)
      }
      if (inherits(r_crop, "try-error") || is.null(r_crop)) {
        .catln("  - Skipping (crop failed).")
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year        = y,
          family      = fam,
          basename    = base,
          in_file     = in_fp,
          tif_written = FALSE,
          asc_written = FALSE,
          notes       = "crop_failed",
          stringsAsFactors = FALSE
        )
        next
      }

      # Ensure something remains after crop
      vals <- try(.t$values(r_crop), silent = TRUE)
      if (inherits(vals, "try-error") || all(is.na(vals))) {
        .catln("  - Skipping (no valid cells after crop).")
        out_rows[[length(out_rows) + 1L]] <- data.frame(
          year        = y,
          family      = fam,
          basename    = base,
          in_file     = in_fp,
          tif_written = FALSE,
          asc_written = FALSE,
          notes       = "no_cells_after_crop",
          stringsAsFactors = FALSE
        )
        next
      }

      # Write outputs
      base_out <- file.path(out_dir, base)
      wrote_this_tif <- FALSE
      wrote_this_asc <- FALSE

      if (".tif" %in% file_type) {
        out_fp_tif <- paste0(base_out, ".tif")
        if (.write_tif(r_crop, out_fp_tif)) {
          wrote_tif <- wrote_tif + 1L
          wrote_this_tif <- TRUE
        }
      }
      if (".asc" %in% file_type) {
        out_fp_asc <- paste0(base_out, ".asc")
        if (.write_asc(r_crop, out_fp_asc)) {
          wrote_asc <- wrote_asc + 1L
          wrote_this_asc <- TRUE
        }
      }

      out_rows[[length(out_rows) + 1L]] <- data.frame(
        year        = y,
        family      = fam,
        basename    = base,
        in_file     = in_fp,
        tif_written = wrote_this_tif,
        asc_written = wrote_this_asc,
        notes       = if (!wrote_this_tif && !wrote_this_asc) "write_failed" else "",
        stringsAsFactors = FALSE
      )
    }

    .catln(sprintf("  - Wrote %d .tif and %d .asc file(s) to %s",
                   wrote_tif, wrote_asc, out_dir))
  }

  out_df <- if (length(out_rows)) {
    do.call(rbind, out_rows)
  } else {
    data.frame(
      year        = integer(0),
      family      = character(0),
      basename    = character(0),
      in_file     = character(0),
      tif_written = logical(0),
      asc_written = logical(0),
      notes       = character(0),
      stringsAsFactors = FALSE
    )
  }

  ## ---- log summary -----------------------------------------------------------
  .append_log(
    "Processing summary (get_merra_variables)",
    c(
      sprintf("Alpha code:          %s", alpha_code),
      sprintf("Extent (WGS84):     xmin=%.6f xmax=%.6f ymin=%.6f ymax=%.6f",
              bb$xmin, bb$xmax, bb$ymin, bb$ymax),
      sprintf("Input m2 root:      %s", m2_root),
      sprintf("Input mc root:      %s", mc_root),
      sprintf("Output root:        %s", vars_dir),
      "Notes:",
      " - UL/LR order auto-detected; longitudes normalized to [-180,180].",
      " - Geographic rasters: rotated 0-360 to -180 to 180 if needed; antimeridian handled.",
      " - Projected rasters: bbox polygon projected to raster CRS for cropping.",
      " - .tif: preserved/assigned CRS; .asc: sidecars removed where possible.",
      " - Wrote any crop with cells (even if stats reported NA)."
    )
  )

  .catln(.sep())
  .catln("Done. Log updated: ", log_file)
  .catln(.sep())

  invisible(out_df)
}
