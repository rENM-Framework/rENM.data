#' Apply spatial thinning to occurrence records (sequential)
#'
#' Applies a greedy sequential thinning algorithm to per-bin occurrence CSV
#' files in \code{_occs/tmp/}, enforcing a minimum separation distance between
#' retained points. Files are overwritten in place.
#'
#' @details
#' The thinning algorithm processes records in input order. The first point is
#' always retained. Each subsequent point is retained only if its Haversine
#' distance to every previously retained point is at least
#' \code{thin_distance} kilometers. This is a greedy nearest-neighbor filter
#' and does not guarantee an optimal solution, but is fast and deterministic
#' for a given input order.
#'
#' Use \code{\link{thin_occurrences2}} for parallel execution on large datasets.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param thin_distance Numeric scalar. Minimum separation distance in
#'   kilometers between retained points. Must be positive.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{files_processed}, \code{points_before} (named integer vector),
#'   \code{points_after} (named integer vector), \code{output_dir}, and
#'   \code{log_file}.
#'
#' @seealso \code{\link{thin_occurrences2}},
#'   \code{\link{remove_duplicate_occurrences}},
#'   \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#' thin_occurrences("CASP", thin_distance = 10)
#' thin_occurrences("CASP", thin_distance = 10, project_dir = "/projects/rENM")
#' }
#'
#' @family occurrence processing
#' @export
thin_occurrences <- function(alpha_code, thin_distance, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(thin_distance) || length(thin_distance) != 1L ||
      !is.finite(thin_distance) || thin_distance <= 0) {
    stop("`thin_distance` must be a positive numeric scalar (km).", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir  <- file.path(project_root, "runs", alpha_code)
  tmp_dir  <- file.path(run_dir, "_occs", "tmp")
  log_fp   <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Occurrence tmp directory not found for '", alpha_code, "'.\n",
      "Expected: ", tmp_dir, "\n",
      "Run prior preprocessing steps first.",
      call. = FALSE
    )
  }

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  ## ---- Haversine distance (km) ----------------------------------------------
  haversine_km <- function(lon1, lat1, lon2, lat2) {
    R      <- 6371
    to_rad <- pi / 180
    dlon   <- (lon2 - lon1) * to_rad
    dlat   <- (lat2 - lat1) * to_rad
    a <- sin(dlat / 2)^2 +
      cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
    2 * R * atan2(sqrt(a), sqrt(1 - a))
  }

  ## ---- process each file ----------------------------------------------------
  processed_files <- character(0)
  before_counts   <- integer(0)
  after_counts    <- integer(0)

  for (f in files) {
    dat <- utils::read.csv(f, stringsAsFactors = FALSE)

    if (!all(c("longitude", "latitude") %in% names(dat))) {
      stop("Missing required columns (longitude, latitude) in: ", f, call. = FALSE)
    }

    n_before <- nrow(dat)
    keep     <- logical(n_before)

    for (i in seq_len(n_before)) {
      if (i == 1L) { keep[i] <- TRUE; next }
      dists  <- sapply(which(keep), function(j) {
        haversine_km(dat$longitude[i], dat$latitude[i],
                     dat$longitude[j], dat$latitude[j])
      })
      keep[i] <- all(dists >= thin_distance)
    }

    dat_thinned <- dat[keep, , drop = FALSE]
    n_after     <- nrow(dat_thinned)

    utils::write.csv(dat_thinned, f, row.names = FALSE)

    processed_files <- c(processed_files, f)
    before_counts   <- c(before_counts, n_before)
    after_counts    <- c(after_counts, n_after)

    message(sprintf("Processed %s | before=%d after=%d", basename(f), n_before, n_after))
  }

  names(before_counts) <- basename(processed_files)
  names(after_counts)  <- basename(processed_files)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (thin_occurrences)", c(
    sprintf("Alpha code:         %s", alpha_code),
    sprintf("Thin distance (km): %s", thin_distance),
    sprintf("Files processed:    %d", length(processed_files)),
    "Counts by file:",
    paste0("  - ", names(before_counts), ": ", before_counts, " -> ", after_counts)
  ))

  invisible(list(
    alpha_code      = alpha_code,
    files_processed = processed_files,
    points_before   = before_counts,
    points_after    = after_counts,
    output_dir      = tmp_dir,
    log_file        = log_fp
  ))
}
