#' Set the spatial extent for a species run
#'
#' Writes an \code{extent.txt} file for a species run from explicit bounding
#' box coordinates. Any existing \code{extent.txt} is renamed to a backup
#' before writing. Use \code{\link{find_occurrence_extent}} or
#' \code{\link{find_range_extent}} to derive the extent automatically from
#' data rather than supplying coordinates manually.
#'
#' @details
#' Coordinates are provided as upper-left and lower-right corners in
#' (longitude, latitude) order. Non-canonical orientation (e.g., upper-left
#' latitude below lower-right latitude) is reported as a warning but not
#' corrected. Valid ranges are longitude -180 to 180 and latitude -90 to 90.
#'
#' An existing \code{extent.txt} is backed up as \code{extent.txt.original},
#' or as a timestamped variant if a backup already exists.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param ul Numeric vector of length 2. Upper-left corner \code{c(lon, lat)}.
#'   Default \code{c(-128.0, 49.0)}.
#' @param lr Numeric vector of length 2. Lower-right corner \code{c(lon, lat)}.
#'   Default \code{c(-66.5, 24.0)}.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{extent} (list with \code{ul} and \code{lr} named numeric vectors),
#'   \code{paths} (list with \code{extent} and \code{log} file paths), and
#'   \code{notes} (character vector of orientation warnings, or \code{NULL}).
#'
#' @seealso \code{\link{find_occurrence_extent}}, \code{\link{find_range_extent}}
#'
#' @examples
#' \dontrun{
#' set_extent("CASP")
#' set_extent("AMRO", ul = c(-130, 55), lr = c(-60, 20))
#' set_extent("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @export
set_extent <- function(alpha_code,
                       ul = c(-128.0, 49.0),
                       lr = c(-66.5,  24.0),
                       project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  .is_num_pair <- function(x) is.numeric(x) && length(x) == 2L && all(is.finite(x))
  if (!.is_num_pair(ul)) stop("`ul` must be a finite numeric vector of length 2: c(lon, lat).", call. = FALSE)
  if (!.is_num_pair(lr)) stop("`lr` must be a finite numeric vector of length 2: c(lon, lat).", call. = FALSE)

  ul <- as.numeric(ul); names(ul) <- c("lon", "lat")
  lr <- as.numeric(lr); names(lr) <- c("lon", "lat")

  in_range <- function(lon, lat) isTRUE(lon >= -180 && lon <= 180 && lat >= -90 && lat <= 90)
  if (!in_range(ul["lon"], ul["lat"]) || !in_range(lr["lon"], lr["lat"])) {
    stop("Coordinates out of bounds. Longitude must be in [-180, 180], latitude in [-90, 90].",
         call. = FALSE)
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

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir   <- .expand(file.path(project_root, "runs", alpha_code))
  occs_dir  <- file.path(run_dir, "_occs")
  log_fp    <- file.path(run_dir, "_log.txt")
  extent_fp <- file.path(occs_dir, "extent.txt")

  .mkdir(occs_dir)

  ## ---- banner ---------------------------------------------------------------
  .catln(.sep_line())
  .catln(sprintf("set_extent: alpha_code=%s", alpha_code))
  .catln(sprintf("UL (lon,lat): (%.6f, %.6f)", ul["lon"], ul["lat"]))
  .catln(sprintf("LR (lon,lat): (%.6f, %.6f)", lr["lon"], lr["lat"]))
  if (length(orientation_notes)) .catln("Note: ", paste(orientation_notes, collapse = "; "))
  .catln(.sep_line())

  ## ---- backup existing extent.txt -------------------------------------------
  backup_fp  <- file.path(occs_dir, "extent.txt.original")
  rotated_to <- NULL

  if (file.exists(extent_fp)) {
    if (file.exists(backup_fp)) {
      ts        <- gsub("[- :]", "", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
      backup_fp <- file.path(occs_dir, paste0("extent.txt.original_", ts))
    }
    ok <- file.rename(extent_fp, backup_fp)
    if (!ok) {
      warning("Could not rename existing extent.txt.", call. = FALSE)
    } else {
      rotated_to <- backup_fp
      .catln("Existing extent.txt backed up to: ", basename(backup_fp))
    }
  }

  ## ---- write extent.txt -----------------------------------------------------
  writeLines(c(
    "# extent.txt",
    "# Computed by set_extent()",
    paste0("# Timestamp: ", .now()),
    "# Source: manual entry via set_extent()",
    "# Points used: n/a",
    "# bbox_pct: n/a",
    "# Coordinate order: (lon, lat)",
    sprintf("Upper-left:  (%.6f, %.6f)", ul["lon"], ul["lat"]),
    sprintf("Lower-right: (%.6f, %.6f)", lr["lon"], lr["lat"])
  ), con = extent_fp)

  .catln("Wrote extent.txt: ", extent_fp)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (set_extent)", c(
    sprintf("Alpha code:    %s", alpha_code),
    sprintf("UL (lon,lat):  (%.6f, %.6f)", ul["lon"], ul["lat"]),
    sprintf("LR (lon,lat):  (%.6f, %.6f)", lr["lon"], lr["lat"]),
    sprintf("Wrote:         %s", extent_fp),
    sprintf("Backed up:     %s", if (!is.null(rotated_to)) rotated_to else "(none)")
  ))

  .catln("Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(list(
    alpha_code = alpha_code,
    extent     = list(ul = ul, lr = lr),
    paths      = list(extent = extent_fp, log = log_fp),
    notes      = if (length(orientation_notes)) orientation_notes else NULL
  ))
}
