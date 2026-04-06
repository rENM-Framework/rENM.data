#' Set the spatial extent to be used in further model processing
#'
#' Create or update an \code{extent.txt} file for a species run,
#' preserving the format expected by downstream functions.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' The extent file is written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/extent.txt}
#'
#' Coordinates are provided as:
#' \itemize{
#'   \item \code{ul}: upper-left (lon, lat)
#'   \item \code{lr}: lower-right (lon, lat)
#' }
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Validates coordinate structure and numeric ranges
#'   \item Checks for canonical orientation of bounding box corners
#'   \item Renames any existing \code{extent.txt} to a backup file
#'   \item Writes a new \code{extent.txt} using standardized formatting
#'   \item Appends a processing summary to the run log
#' }
#'
#' Coordinates must satisfy:
#' \itemize{
#'   \item longitude in (-180, 180)
#'   \item latitude in (-90, 90)
#' }
#'
#' Orientation checks:
#' \itemize{
#'   \item Upper-left latitude should exceed lower-right latitude
#'   \item Upper-left longitude should be less than lower-right longitude
#' }
#'
#' Non-canonical orientation is reported but not corrected.
#'
#' \strong{Outputs}:
#' The \code{extent.txt} file contains:
#'
#' \itemize{
#'   \item Header metadata (timestamp, source, coordinate order)
#'   \item Upper-left and lower-right coordinates in (lon, lat) format
#' }
#'
#' Existing files are backed up as:
#'
#' \code{extent.txt.original}
#'
#' or a timestamped variant if a backup already exists.
#'
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' \strong{Data requirements}:
#' Coordinates must be numeric length-2 vectors with finite values.
#'
#' @param alpha_code Character(1). Species alpha code used to locate the
#'   run directory under \code{<project_dir>/runs/<alpha_code>/}.
#' @param ul Numeric(2). Upper-left corner (lon, lat) of the bounding box.
#'   Default \code{c(-128.0, 49.0)}.
#' @param lr Numeric(2). Lower-right corner (lon, lat) of the bounding box.
#'   Default \code{c(-66.5, 24.0)}.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{alpha_code}: Species code
#'   \item \code{extent}: List with \code{ul} and \code{lr} coordinate pairs
#'   \item \code{paths}: List with paths to \code{extent.txt} and log file
#'   \item \code{notes}: Character vector of orientation warnings or NULL
#' }
#'
#' @examples
#' \dontrun{
#'   set_extent("CASP")
#'   set_extent("AMRO", ul = c(-130, 55), lr = c(-60, 20))
#' }
#'
#' @export
set_extent <- function(alpha_code,
                       ul = c(-128.0, 49.0),
                       lr = c(-66.5, 24.0)) {

  ## ---- helpers ----------------------------------------------------------------
  .expand   <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  .mkdir    <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  .now      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep_line <- function(n = 72L) paste(rep.int("-", n), collapse = "")
  .catln    <- function(...) { cat(paste0(..., "\n")); flush.console() }

  .append_log <- function(log_fp, title, lines) {
    dir.create(dirname(log_fp), recursive = TRUE, showWarnings = FALSE)
    section <- c(
      "",
      .sep_line(),
      title,
      paste0("Timestamp: ", .now()),
      lines,
      ""
    )
    cat(paste0(section, collapse = "\n"), file = log_fp, append = TRUE)
  }

  ## ---- validation -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  .is_num_pair <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  if (!.is_num_pair(ul)) stop("`ul` must be numeric c(lon, lat).", call. = FALSE)
  if (!.is_num_pair(lr)) stop("`lr` must be numeric c(lon, lat).", call. = FALSE)

  ul <- as.numeric(ul); names(ul) <- c("lon", "lat")
  lr <- as.numeric(lr); names(lr) <- c("lon", "lat")

  in_range <- function(lon, lat) isTRUE(lon >= -180 && lon <= 180 && lat >= -90 && lat <= 90)
  if (!in_range(ul["lon"], ul["lat"]) || !in_range(lr["lon"], lr["lat"])) {
    stop("Coordinates out of bounds. lon must be in [-180, 180], lat in [-90, 90].", call. = FALSE)
  }

  orientation_notes <- character(0)
  if (ul["lat"] < lr["lat"]) {
    orientation_notes <- c(orientation_notes,
                           "Upper-left latitude is less than lower-right latitude (non-standard orientation).")
  }
  if (ul["lon"] > lr["lon"]) {
    orientation_notes <- c(orientation_notes,
                           "Upper-left longitude is greater than lower-right longitude (non-standard orientation).")
  }

  ## ---- paths ------------------------------------------------------------------
  project_dir <- rENM_project_dir()
  run_dir  <- .expand(file.path(project_dir, "runs", alpha_code))
  occs_dir <- file.path(run_dir, "_occs")
  log_fp   <- file.path(run_dir, "_log.txt")
  extent_fp <- file.path(occs_dir, "extent.txt")

  .mkdir(occs_dir)

  ## ---- console banner ---------------------------------------------------------
  .catln(.sep_line())
  .catln(sprintf("set_extent: alpha_code=%s", alpha_code))
  .catln(sprintf("UL (lon,lat): (%.6f, %.6f)", ul["lon"], ul["lat"]))
  .catln(sprintf("LR (lon,lat): (%.6f, %.6f)", lr["lon"], lr["lat"]))
  if (length(orientation_notes)) {
    .catln("Note: ", paste(orientation_notes, collapse = "; "))
  }
  .catln("Writing directory: ", occs_dir)
  .catln(.sep_line())

  ## ---- rotate existing extent.txt --------------------------------------------
  backup_fp <- file.path(occs_dir, "extent.txt.original")
  rotated_to <- NULL
  did_rotate <- FALSE

  if (file.exists(extent_fp)) {
    if (file.exists(backup_fp)) {
      timestamp <- gsub("[- :]", "", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
      backup_fp <- file.path(occs_dir, paste0("extent.txt.original_", timestamp))
    }
    ok <- file.rename(extent_fp, backup_fp)
    if (!ok) {
      warning("Could not rename existing extent.txt.", call. = FALSE)
    } else {
      rotated_to <- backup_fp
      did_rotate <- TRUE
      .catln("Existing extent.txt found -> renamed to: ", basename(backup_fp))
    }
  }

  ## ---- write extent.txt -------------------------------------------------------
  header <- c(
    "# extent.txt",
    "# Computed by set_extent()",
    paste0("# Timestamp: ", .now()),
    "# Source: manual entry via set_extent()",
    "# Points used: n/a",
    "# bbox_pct: n/a",
    "# Coordinate order: (lon, lat)"
  )
  body <- c(
    sprintf("Upper-left:  (%.6f, %.6f)", ul["lon"], ul["lat"]),
    sprintf("Lower-right: (%.6f, %.6f)", lr["lon"], lr["lat"])
  )

  writeLines(c(header, body), con = extent_fp)
  .catln("Wrote new extent.txt to: ", extent_fp)

  ## ---- log --------------------------------------------------------------------
  details <- c(
    sprintf("Alpha code:          %s", alpha_code),
    sprintf("UL (lon,lat):        (%.6f, %.6f)", ul["lon"], ul["lat"]),
    sprintf("LR (lon,lat):        (%.6f, %.6f)", lr["lon"], lr["lat"]),
    sprintf("Wrote:               %s", extent_fp),
    sprintf("Rotated prior:       %s", if (did_rotate) rotated_to else "(none)")
  )
  .append_log(log_fp, "Processing summary (set_extent)", details)
  .catln("Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(list(
    alpha_code = alpha_code,
    extent = list(ul = ul, lr = lr),
    paths = list(extent = extent_fp, log = log_fp),
    notes = if (length(orientation_notes)) orientation_notes else NULL
  ))
}
