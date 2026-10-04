#' Clip occurrence records to the modeled extent
#'
#' Removes, from each per-bin occurrence CSV in \code{_occs/tmp/}, every
#' record that falls outside the extent box written to
#' \code{_occs/extent.txt}, and overwrites each file in place.
#'
#' @details
#' The predictor rasters cover only the extent box, so a record outside it
#' has no predictor values and \code{sdm} drops it when a model is fitted.
#' Removing such records before thinning and capping keeps them from using
#' up the per-bin record cap. Without this step, a species whose eBird
#' records extend well beyond its GAP range (into Mexico or Canada, for
#' example) can reach its models with far fewer records per bin than the cap
#' suggests.
#'
#' The clip uses the box, not the buffered range polygon, because the box is
#' the area the predictors cover. Records on cells that are \code{NA} inside
#' the box, such as open water, are still dropped later by \code{sdm}.
#'
#' Run this after \code{find_range_extent()} (or \code{set_extent()}) and
#' before \code{thin_occurrences2()} and \code{limit_record_count()}.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{files_processed} (character vector of processed file paths),
#'   \code{counts_before} and \code{counts_after} (named integer vectors),
#'   \code{extent} (named numeric vector: \code{xmin}, \code{xmax},
#'   \code{ymin}, \code{ymax}), \code{output_dir}, and \code{log_file}.
#'
#' @seealso \code{\link{find_range_extent}}, \code{\link{thin_occurrences2}},
#'   \code{\link{limit_record_count}}, \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#' clip_occurrences("CASP")
#' clip_occurrences("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @family occurrence processing
#' @export
clip_occurrences <- function(alpha_code, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir   <- file.path(project_root, "runs", alpha_code)
  tmp_dir   <- file.path(run_dir, "_occs", "tmp")
  extent_fp <- file.path(run_dir, "_occs", "extent.txt")
  log_fp    <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Occurrence tmp directory not found for '", alpha_code, "'.\n",
      "Expected: ", tmp_dir, "\n",
      "Run get_ebird_occurrences() first.",
      call. = FALSE
    )
  }
  if (!file.exists(extent_fp)) {
    stop(
      "extent.txt not found at: ", extent_fp, "\n",
      "Run find_range_extent() or set_extent() first.",
      call. = FALSE
    )
  }

  ## ---- parse extent.txt -----------------------------------------------------
  # Same format get_merra_variables() reads: "Upper-left:  (lon, lat)" and
  # "Lower-right: (lon, lat)".
  .parse_pair <- function(s) {
    m <- regexec("\\(([-0-9.]+)\\s*,\\s*([-0-9.]+)\\)", s, perl = TRUE)
    v <- regmatches(s, m)
    if (length(v) == 1L && length(v[[1L]]) == 3L) as.numeric(v[[1L]][2:3])
    else c(NA_real_, NA_real_)
  }

  txt <- readLines(extent_fp, warn = FALSE)
  ul  <- .parse_pair(grep("Upper-left",  txt, ignore.case = TRUE, value = TRUE)[1L])
  lr  <- .parse_pair(grep("Lower-right", txt, ignore.case = TRUE, value = TRUE)[1L])
  if (any(is.na(c(ul, lr)))) {
    stop("Failed to parse coordinates from extent.txt.", call. = FALSE)
  }

  ext <- c(xmin = min(ul[1L], lr[1L]), xmax = max(ul[1L], lr[1L]),
           ymin = min(ul[2L], lr[2L]), ymax = max(ul[2L], lr[2L]))

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  ## ---- process each file ----------------------------------------------------
  counts_before <- integer(0)
  counts_after  <- integer(0)

  for (f in files) {
    dat <- tryCatch(
      utils::read.csv(f, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Failed to read: ", f, " | ", conditionMessage(e), call. = FALSE)
      }
    )

    if (!all(c("longitude", "latitude") %in% names(dat))) {
      stop("Missing required columns (longitude, latitude) in: ", f, call. = FALSE)
    }

    # Inclusive bounds, matching how a point on the box edge falls inside the
    # cropped predictor rasters.
    keep <- !is.na(dat$longitude) & !is.na(dat$latitude) &
      dat$longitude >= ext[["xmin"]] & dat$longitude <= ext[["xmax"]] &
      dat$latitude  >= ext[["ymin"]] & dat$latitude  <= ext[["ymax"]]

    utils::write.csv(dat[keep, , drop = FALSE], f, row.names = FALSE)

    counts_before <- c(counts_before, nrow(dat))
    counts_after  <- c(counts_after, sum(keep))

    message(sprintf("Processed %s | kept %d of %d inside the extent",
                    basename(f), sum(keep), nrow(dat)))
    if (!any(keep)) {
      warning("No records inside the extent in ", basename(f), call. = FALSE)
    }
  }

  names(counts_before) <- basename(files)
  names(counts_after)  <- basename(files)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (clip_occurrences)", c(
    sprintf("Alpha code:      %s", alpha_code),
    sprintf("Extent:          lon %.6f to %.6f, lat %.6f to %.6f",
            ext[["xmin"]], ext[["xmax"]], ext[["ymin"]], ext[["ymax"]]),
    sprintf("Files processed: %d", length(files)),
    "Records inside the extent (kept / read):",
    paste0("  - ", names(counts_after), ": ", counts_after, " / ", counts_before)
  ))

  invisible(list(
    alpha_code      = alpha_code,
    files_processed = files,
    counts_before   = counts_before,
    counts_after    = counts_after,
    extent          = ext,
    output_dir      = tmp_dir,
    log_file        = log_fp
  ))
}
