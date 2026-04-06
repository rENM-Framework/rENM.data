#' Apply spatial thinning to occurrence records (sequential version)
#'
#' Applies spatial thinning to occurrence records to reduce spatial
#' clustering by enforcing a minimum distance between retained points.
#' This helps mitigate spatial sampling bias in ecological niche models.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' Spatial thinning is typically applied after duplicate removal and
#' prior to model fitting. The thinning algorithm enforces a minimum
#' separation distance between retained occurrence points, reducing
#' over-representation of heavily sampled locations.
#'
#' The rENM project directory is resolved using
#' \code{\link{rENM_project_dir}}, ensuring that all filesystem operations
#' are portable and do not rely on hard-coded paths such as ~/rENM.
#'
#' \strong{Inputs}
#' This function operates on per-bin occurrence CSV files located in:
#'
#' \code{<rENM_project_dir>/runs/<alpha_code>/_occs/tmp/}
#'
#' Files must follow the naming pattern:
#'
#' \code{of-<year>.csv}
#'
#' Each file must contain columns:
#'
#' \itemize{
#'   \item longitude
#'   \item latitude
#' }
#'
#' \strong{Outputs}
#' Thinned datasets are written back to the same directory:
#'
#' \code{<rENM_project_dir>/runs/<alpha_code>/_occs/tmp/}
#'
#' Files are overwritten in place.
#'
#' \strong{Methods}
#' A sequential thinning algorithm is applied:
#'
#' \itemize{
#'   \item Records are processed in input order
#'   \item The first point is always retained
#'   \item Each subsequent point is compared to all retained points
#'   \item Distances are computed using the Haversine formula
#'   \item Points closer than the specified threshold are removed
#' }
#'
#' The Haversine distance is computed in kilometers using a spherical
#' Earth approximation.
#'
#' \strong{Data requirements}
#' The temporary occurrence directory must exist and contain one or more
#' CSV files matching the expected naming pattern and column structure.
#'
#' Centralizing project-directory resolution via
#' \code{\link{rENM_project_dir}} ensures that all functions in the rENM
#' package operate relative to a single configurable project root,
#' improving reproducibility and enabling use in HPC or multi-project
#' environments.
#'
#' @param alpha_code Character. Four-letter species alpha code (e.g.,
#' "CASP").
#'
#' @param thin_distance Numeric. Minimum allowed distance (in kilometers)
#' between retained occurrence points.
#'
#' @param project_dir Character. Optional rENM project root directory. If
#' NULL, the directory is resolved using
#' \code{\link{rENM_project_dir}} via:
#'
#' \enumerate{
#'   \item Explicit project_dir argument
#'   \item options(rENM.project_dir)
#'   \item RENM_PROJECT_DIR environment variable
#' }
#'
#' Providing project_dir explicitly is recommended for reproducible
#' scripts and package testing.
#'
#' @return
#' List. Invisibly returns a structured list containing:
#'
#' \describe{
#'   \item{alpha_code}{Character. Species alpha code.}
#'   \item{files_processed}{Character vector. Full paths to processed
#'   occurrence files.}
#'   \item{points_before}{Named integer vector. Number of occurrence
#'   points before thinning for each file.}
#'   \item{points_after}{Named integer vector. Number of occurrence
#'   points after thinning for each file.}
#'   \item{output_dir}{Character. Directory containing thinned occurrence
#'   files.}
#'   \item{log_file}{Character. Path to the _log.txt file updated with
#'   processing details.}
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Overwrites occurrence CSV files with thinned datasets
#'   \item Reduces point counts according to the thinning distance
#'   \item Appends a processing summary to the log file
#' }
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#'
#' # Standard workflow
#' thin_occurrences("CASP", thin_distance = 10)
#'
#' # Reproducible script
#' thin_occurrences("CASP", thin_distance = 10,
#'                  project_dir = "/projects/rENM")
#'
#' }
#'
#' @seealso
#' \link{rENM_project_dir}
#'
#' @export
thin_occurrences <- function(alpha_code,
                             thin_distance,
                             project_dir = NULL) {

  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  if (!is.numeric(thin_distance) || length(thin_distance) != 1L || thin_distance <= 0) {
    stop("`thin_distance` must be a positive numeric scalar (km).", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  # ---------------------------------------------------------------------------
  # Resolve project directory (CRAN-compliant)
  # ---------------------------------------------------------------------------
  # All filesystem paths are derived from this root. This replaces any
  # previous hard-coded paths such as ~/rENM and ensures portability.
  project_root <- rENM_project_dir(project_dir)

  # ---------------------------------------------------------------------------
  # Define directory structure
  # ---------------------------------------------------------------------------
  run_dir  <- file.path(project_root, "runs", alpha_code)
  occ_dir  <- file.path(run_dir, "_occs")
  tmp_dir  <- file.path(occ_dir, "tmp")
  log_file <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Occurrence tmp directory not found for alpha_code '", alpha_code, "'.\n",
      "Expected: ", tmp_dir, "\n",
      "Run prior preprocessing steps first.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Identify occurrence files
  # ---------------------------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helper: Haversine distance (km)
  # ---------------------------------------------------------------------------
  haversine_km <- function(lon1, lat1, lon2, lat2) {
    R <- 6371
    to_rad <- pi / 180
    dlon <- (lon2 - lon1) * to_rad
    dlat <- (lat2 - lat1) * to_rad

    a <- sin(dlat/2)^2 +
      cos(lat1 * to_rad) * cos(lat2 * to_rad) *
      sin(dlon/2)^2

    2 * R * atan2(sqrt(a), sqrt(1 - a))
  }

  processed_files <- character(0)
  before_counts   <- integer(0)
  after_counts    <- integer(0)

  # ---------------------------------------------------------------------------
  # Process each file
  # ---------------------------------------------------------------------------
  for (f in files) {

    dat <- utils::read.csv(f, stringsAsFactors = FALSE)

    if (!all(c("longitude", "latitude") %in% names(dat))) {
      stop("Missing required columns in file: ", f, call. = FALSE)
    }

    n_before <- nrow(dat)

    keep <- logical(n_before)

    for (i in seq_len(n_before)) {

      if (i == 1) {
        keep[i] <- TRUE
        next
      }

      dists <- sapply(which(keep), function(j) {
        haversine_km(
          dat$longitude[i], dat$latitude[i],
          dat$longitude[j], dat$latitude[j]
        )
      })

      keep[i] <- all(dists >= thin_distance)
    }

    dat_thinned <- dat[keep, , drop = FALSE]

    n_after <- nrow(dat_thinned)

    utils::write.csv(dat_thinned, f, row.names = FALSE)

    processed_files <- c(processed_files, f)
    before_counts   <- c(before_counts, n_before)
    after_counts    <- c(after_counts, n_after)

    message(sprintf(
      "Processed %s | before=%d after=%d",
      basename(f), n_before, n_after
    ))
  }

  names(before_counts) <- basename(processed_files)
  names(after_counts)  <- basename(processed_files)

  # ---------------------------------------------------------------------------
  # Logging
  # ---------------------------------------------------------------------------
  .now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep <- function(n = 72L) paste(rep("-", n), collapse = "")

  log_lines <- c(
    "",
    .sep(),
    "Processing summary (thin_occurrences)",
    paste0("Timestamp: ", .now()),
    paste0("Alpha code: ", alpha_code),
    paste0("Thin distance (km): ", thin_distance),
    paste0("Files processed: ", length(processed_files)),
    "Counts by file:",
    paste0("  - ", names(before_counts),
           ": ", before_counts, " -> ", after_counts)
  )

  cat(paste(log_lines, collapse = "\n"), file = log_file, append = TRUE)

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  invisible(list(
    alpha_code = alpha_code,
    files_processed = processed_files,
    points_before = before_counts,
    points_after = after_counts,
    output_dir = tmp_dir,
    log_file = log_file
  ))
}
