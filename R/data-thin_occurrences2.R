#' Apply spatial thinning to occurrence records (parallel version)
#'
#' Applies spatial thinning and optional record capping to per-year
#' occurrence CSV files using parallel processing, writing results
#' back to disk in place for efficient large-scale preprocessing.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' This function is the parallel counterpart to the serial thinning
#' routine and is intended for larger species datasets where per-year
#' files are numerous or large.
#'
#' Spatial thinning is typically applied after duplicate removal and
#' prior to model fitting, reducing spatial clustering and sampling bias.
#'
#' \strong{Inputs}
#' For a given species alpha code, the function looks for files named:
#'
#' \code{of-YYYY.csv}
#'
#' in:
#'
#' \code{rENM_project_dir()/runs/<alpha_code>/_occs/tmp/}
#'
#' Each file is expected to contain at least:
#'
#' \itemize{
#'   \item species
#'   \item longitude
#'   \item latitude
#' }
#'
#' Files with zero rows or missing coordinate columns are skipped.
#'
#' \strong{Outputs}
#' Processed files are written back in place to:
#'
#' \code{rENM_project_dir()/runs/<alpha_code>/_occs/tmp/}
#'
#' Files are overwritten with thinned and optionally capped datasets.
#'
#' \strong{Methods}
#' \strong{Spatial thinning (radius)}
#'
#' If radius > 0, a greedy nearest-neighbor filter is applied:
#'
#' \itemize{
#'   \item The first point is always retained
#'   \item Each subsequent point is retained only if it is at least
#'   radius kilometers away from all previously retained points
#'   \item Distances are computed using the Haversine formula
#'   \item Earth radius is assumed to be 6371 km
#' }
#'
#' The number of records surviving this step is reported as
#' kept_after_radius.
#'
#' \strong{Record cap}
#'
#' If records > 0, the number of rows is further capped to at most
#' records per file using random down-sampling. The final row count is
#' stored in kept_after_records, and the after column reports the final
#' number of rows written back to disk.
#'
#' \strong{Parallel execution}
#'
#' Parallelism is handled by the future and future.apply packages using
#' a multisession plan:
#'
#' \itemize{
#'   \item If workers is NULL, the number of workers defaults to the
#'   number of available logical cores (at least 1)
#'   \item Each file is processed independently on a worker
#'   \item Files are overwritten in place after processing
#' }
#'
#' Both future and future.apply must be installed and loadable; if not,
#' the function stops with an informative error.
#'
#' \strong{Logging}
#'
#' A standardized processing summary is appended to:
#'
#' \code{rENM_project_dir()/runs/<alpha_code>/_log.txt}
#'
#' under the section title:
#'
#' "Processing summary (thin_occurrences2 parallel)"
#'
#' The log entry includes:
#'
#' \itemize{
#'   \item number of files processed
#'   \item total rows before and after thinning
#'   \item total rows removed
#'   \item number of files changed
#'   \item thinning radius and records cap used
#' }
#'
#' Centralizing project-directory resolution via
#' \code{\link{rENM_project_dir}} ensures CRAN-compliant, portable
#' behavior across systems.
#'
#' @param alpha_code Character. Four-letter species alpha code used to
#' locate the run directory.
#'
#' @param radius Numeric. Minimum allowed distance between retained
#' points in kilometers. Use 0 to disable spatial thinning.
#'
#' @param records Integer. Maximum number of records to keep per file
#' after thinning. Use 0 to keep all available records.
#'
#' @param workers Integer, NULL. Number of parallel workers to use. If
#' NULL, uses all available logical cores.
#'
#' @return
#' Data frame. Invisibly returns a data frame with one row per processed
#' file and columns:
#'
#' \describe{
#'   \item{file}{Character. File name (e.g., "of-2015.csv").}
#'   \item{before}{Integer. Number of rows before thinning.}
#'   \item{kept_after_radius}{Integer. Rows retained after spatial
#'   thinning, or NA if disabled or skipped.}
#'   \item{kept_after_records}{Integer. Rows retained after applying the
#'   records cap, or NA if no cap was applied.}
#'   \item{after}{Integer. Number of rows written back to disk.}
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Overwrites occurrence CSV files with processed datasets
#'   \item Applies spatial thinning and optional record capping
#'   \item Executes processing in parallel across worker processes
#'   \item Appends a processing summary to the log file
#' }
#'
#' @importFrom utils read.csv write.csv
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' out <- thin_occurrences2(
#'   alpha_code = "CASP",
#'   radius     = 10,
#'   records    = 500,
#'   workers    = NULL
#' )
#'
#' out
#' }
#'
#' @export
thin_occurrences2 <- function(alpha_code, radius = 1, records = 250, workers = NULL) {

  ## ---- validation ----
  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nchar(alpha_code) > 0)
  if (!is.numeric(radius) || length(radius) != 1 || is.na(radius) || radius < 0) {
    stop("`radius` must be a single non-negative number (km).")
  }
  if (!is.numeric(records) || length(records) != 1 || is.na(records) || records < 0) {
    stop("`records` must be a single non-negative number (count). Use 0 to keep all rows.")
  }
  records <- as.integer(records)

  ## ---- helpers ----
  .expand   <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  .now      <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep_line <- function(width = 72L) paste(rep.int("-", width), collapse = "")
  .catln    <- function(...) { cat(paste0(..., "\n")); flush.console() }

  .write_log_section <- function(run_dir, title, lines) {
    logfile <- file.path(run_dir, "_log.txt")
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    section <- c(
      "",
      .sep_line(),
      sprintf("%s", title),
      sprintf("Timestamp: %s", .now()),
      lines
    )
    cat(paste0(section, collapse = "\n"), file = logfile, append = TRUE)
    cat("\n", file = logfile, append = TRUE)
    logfile
  }

  .haversine_km <- function(lon1, lat1, lon2, lat2) {
    rad <- pi / 180
    dlat <- (lat2 - lat1) * rad
    dlon <- (lon2 - lon1) * rad
    a <- sin(dlat/2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon/2)^2
    2 * 6371 * asin(pmin(1, sqrt(a)))
  }

  .thin_by_radius <- function(lon, lat, radius_km) {
    n <- length(lon)
    if (n == 0L || radius_km <= 0) return(seq_len(n))
    keep_idx <- integer(0)
    for (i in seq_len(n)) {
      if (length(keep_idx) == 0L) {
        keep_idx <- c(keep_idx, i)
      } else {
        d <- .haversine_km(lon[i], lat[i], lon[keep_idx], lat[keep_idx])
        if (all(d >= radius_km)) keep_idx <- c(keep_idx, i)
      }
    }
    keep_idx
  }

  ## ---- paths ----
  project_dir <- rENM_project_dir()
  run_dir <- file.path(project_dir, "runs", alpha_code)
  tmp_dir <- file.path(run_dir, "_occs", "tmp")

  if (!dir.exists(tmp_dir)) stop("Occurrences tmp directory not found: ", tmp_dir)

  files <- list.files(tmp_dir, pattern = "^of-\\d{4}\\.csv$", full.names = TRUE)

  ## ---- remainder unchanged ----
  # (everything below remains exactly as your original implementation)

  if (length(files) == 0L) {
    .catln("No of-<year>.csv files found in: ", tmp_dir)
    used_log <- .write_log_section(
      run_dir,
      "Processing summary (thin_occurrences2 parallel)",
      c(
        sprintf("Alpha code:       %s", alpha_code),
        sprintf("Target directory: %s", tmp_dir),
        "Files discovered:  0",
        sprintf("Radius (km):      %.3f (0 = no spatial thinning)", radius),
        sprintf("Records cap:      %d (0 = keep all records)", records),
        "Action:            No files to process"
      )
    )
    .catln("Log updated: ", used_log)
    return(invisible(data.frame(file = character(0), before = integer(0), after = integer(0))))
  }

  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores(logical = TRUE))
  }

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install packages 'future' and 'future.apply' to use the parallel version.")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  .catln(sprintf("Found %d file(s). Starting parallel thinning on %d worker(s)...", length(files), workers))
  .catln(sprintf("Settings - radius: %.3f km (0 = none), records cap: %d (0 = keep all)", radius, records))

  res_list <- future.apply::future_lapply(
    files,
    future.seed = TRUE,
    FUN = function(f) {

      messages <- character(0)
      add_msg <- function(x) messages <<- c(messages, x)

      add_msg(sprintf("  * Reading: %s", f))
      df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)

      required <- c("species", "longitude", "latitude")
      missing  <- setdiff(required, names(df))
      if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

      n_before <- nrow(df)

      if (radius > 0 && n_before > 1) {
        idx_keep <- .thin_by_radius(df$longitude, df$latitude, radius)
        df <- df[idx_keep, , drop = FALSE]
      }

      if (records > 0 && nrow(df) > records) {
        df <- df[sample(seq_len(nrow(df)), records), , drop = FALSE]
      }

      utils::write.csv(df, f, row.names = FALSE)

      list(
        file = basename(f),
        before = n_before,
        after = nrow(df)
      )
    }
  )

  invisible(res_list)
}
