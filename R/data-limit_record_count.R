#' Limit occurrence records per bin by random sampling
#'
#' Randomly downsamples occurrence records within each per-bin CSV file
#' to a maximum of \code{record_count} rows while preserving all records
#' when counts fall below the specified threshold.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' This function operates on occurrence files located in:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/tmp/}
#'
#' Each file follows the naming pattern:
#'
#' \code{of-<year>.csv}
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Reads each \code{of-<year>.csv} file
#'   \item Counts the number of rows
#'   \item Randomly samples up to \code{record_count} rows without
#'         replacement
#'   \item Writes the sampled dataset back to the same file
#' }
#'
#' Sampling uses base R \code{sample()}, allowing reproducibility when a
#' random seed is set externally (e.g., \code{set.seed()}).
#'
#' \strong{Outputs}:
#' Each input file is overwritten with the sampled dataset in:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/tmp/}
#'
#' \strong{Log file}:
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' including total counts before and after sampling and per-file
#' summaries.
#'
#' \strong{Typical use cases}:
#' \itemize{
#'   \item Balancing sample sizes across temporal bins
#'   \item Sensitivity analysis of model performance vs. sample size
#'   \item Reducing computational load for downstream modeling
#' }
#'
#' @param alpha_code Character. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param record_count Integer. Maximum number of records to retain per
#'   bin. Must be positive.
#' @param project_dir Character, NULL. Optional project root directory.
#'   If NULL, resolved via \code{rENM_project_dir()}.
#'
#' @return List (invisibly returned) with elements:
#' \itemize{
#'   \item \code{alpha_code}: Character species alpha code
#'   \item \code{files_processed}: Character vector of processed file paths
#'   \item \code{counts_before}: Named integer vector of original counts
#'   \item \code{counts_after}: Named integer vector of retained counts
#'   \item \code{output_dir}: Character directory containing processed files
#'   \item \code{log_file}: Character path to the run log file
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Overwrites occurrence CSV files in the tmp directory
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#'   limit_record_count("CASP", record_count = 1000)
#'
#'   limit_record_count(
#'     "CASP",
#'     record_count = 500,
#'     project_dir = "/projects/rENM"
#'   )
#' }
#'
#' @seealso \code{get_ebird_occurrences()},
#'   \code{remove_duplicate_occurrences()},
#'   \code{thin_occurrences()},
#'   \code{tidy_occurrences()}
#'
#' @family occurrence processing
#' @export
limit_record_count <- function(alpha_code,
                               record_count,
                               project_dir = NULL) {

  ## ---- validation -----------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  if (!is.numeric(record_count) || length(record_count) != 1L ||
      !is.finite(record_count) || record_count <= 0) {
    stop("`record_count` must be a positive numeric(1).", call. = FALSE)
  }

  alpha_code  <- toupper(trimws(alpha_code))
  record_count <- as.integer(record_count)

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)

  run_dir  <- file.path(project_root, "runs", alpha_code)
  tmp_dir  <- file.path(run_dir, "_occs", "tmp")
  log_file <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Occurrence tmp directory not found for alpha_code '", alpha_code, "'.\n",
      "Expected: ", tmp_dir,
      call. = FALSE
    )
  }

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  ## ---- helpers --------------------------------------------------------------
  .now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep <- function(n = 72L) paste(rep("-", n), collapse = "")
  .catln <- function(...) { cat(paste0(..., "\n")); flush.console() }

  ## ---- processing -----------------------------------------------------------
  processed_files <- character(0)
  before_counts   <- integer(0)
  after_counts    <- integer(0)

  .catln(.sep())
  .catln("limit_record_count: alpha_code=", alpha_code)
  .catln("Record count limit: ", record_count)
  .catln("Directory: ", tmp_dir)
  .catln(.sep())

  for (f in files) {

    dat <- tryCatch(
      utils::read.csv(f, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Failed to read file: ", f, " | ", conditionMessage(e), call. = FALSE)
      }
    )

    n_before <- nrow(dat)

    if (n_before > record_count) {
      idx <- sample.int(n_before, size = record_count, replace = FALSE)
      dat <- dat[idx, , drop = FALSE]
    }

    n_after <- nrow(dat)

    utils::write.csv(dat, f, row.names = FALSE)

    processed_files <- c(processed_files, f)
    before_counts   <- c(before_counts, n_before)
    after_counts    <- c(after_counts, n_after)

    .catln("  ", basename(f), " | before=", n_before, " after=", n_after)
  }

  names(before_counts) <- basename(processed_files)
  names(after_counts)  <- basename(processed_files)

  total_before <- sum(before_counts)
  total_after  <- sum(after_counts)

  ## ---- console summary ------------------------------------------------------
  .catln(.sep())
  .catln("Done.")
  .catln("Total before: ", total_before)
  .catln("Total after:  ", total_after)
  .catln("Log file: ", log_file)
  .catln(.sep())

  ## ---- logging --------------------------------------------------------------
  log_lines <- c(
    "",
    .sep(),
    "Processing summary (limit_record_count)",
    paste0("Timestamp:          ", .now()),
    paste0("Alpha code:         ", alpha_code),
    paste0("Record count limit: ", record_count),
    paste0("Files processed:    ", length(processed_files)),
    paste0("Total before:       ", total_before),
    paste0("Total after:        ", total_after),
    "Per-file counts:",
    paste0("  - ", names(before_counts),
           ": ", before_counts, " -> ", after_counts)
  )

  cat(paste(log_lines, collapse = "\n"), file = log_file, append = TRUE)

  ## ---- return ---------------------------------------------------------------
  invisible(list(
    alpha_code      = alpha_code,
    files_processed = processed_files,
    counts_before   = before_counts,
    counts_after    = after_counts,
    output_dir      = tmp_dir,
    log_file        = log_file
  ))
}
