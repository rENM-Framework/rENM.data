#' Remove duplicate occurrence records from per-bin CSV files
#'
#' Reads each per-bin occurrence CSV in \code{_occs/tmp/}, removes rows with
#' identical longitude and latitude values (retaining the first occurrence of
#' each unique coordinate pair), and overwrites the file in place.
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
#'   \code{duplicates_removed} (named integer vector of counts per file),
#'   \code{output_dir}, and \code{log_file}.
#'
#' @seealso \code{\link{get_ebird_occurrences}}, \code{\link{thin_occurrences}},
#'   \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#' remove_duplicate_occurrences("CASP")
#' remove_duplicate_occurrences("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @family occurrence processing
#' @export
remove_duplicate_occurrences <- function(alpha_code, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
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
      "Run get_ebird_occurrences() first.",
      call. = FALSE
    )
  }

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  ## ---- process each file ----------------------------------------------------
  processed_files <- character(0)
  removed_counts  <- integer(0)

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

    n_before   <- nrow(dat)
    dat_unique <- dat[!duplicated(dat[, c("longitude", "latitude")]), ]
    n_removed  <- n_before - nrow(dat_unique)

    utils::write.csv(dat_unique, f, row.names = FALSE)

    processed_files <- c(processed_files, f)
    removed_counts  <- c(removed_counts, n_removed)

    message(sprintf("Processed %s | removed %d duplicate(s)", basename(f), n_removed))
  }

  names(removed_counts) <- basename(processed_files)

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (remove_duplicate_occurrences)", c(
    sprintf("Alpha code:      %s", alpha_code),
    sprintf("Files processed: %d", length(processed_files)),
    "Duplicates removed by file:",
    paste0("  - ", names(removed_counts), ": ", removed_counts)
  ))

  invisible(list(
    alpha_code         = alpha_code,
    files_processed    = processed_files,
    duplicates_removed = removed_counts,
    output_dir         = tmp_dir,
    log_file           = log_fp
  ))
}
