#' Remove duplicate occurrence records from EBD occurrence datasets
#'
#' Removes duplicate spatial records (longitude/latitude pairs) from
#' per-bin occurrence CSV files generated during earlier rENM
#' preprocessing steps.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' Occurrence files are read from:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/tmp/}
#'
#' Each file is expected to follow the naming pattern:
#'
#' \itemize{
#'   \item \code{of-<year>.csv}
#' }
#'
#' and contain columns named \code{longitude} and \code{latitude}.
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Reads each per-bin CSV file
#'   \item Identifies duplicate rows based on identical longitude and latitude
#'   \item Retains the first occurrence of each unique coordinate pair
#'   \item Overwrites the original file with the de-duplicated dataset
#' }
#'
#' Duplicate records are defined strictly as rows sharing identical
#' longitude and latitude values.
#'
#' \strong{Outputs}:
#' Cleaned occurrence files are written back to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/tmp/}
#'
#' replacing the original files.
#'
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' under the section:
#'
#' \code{Processing summary (remove_duplicate_occurrences)}
#'
#' \strong{Workflow context}:
#' This function is typically used after temporal binning and prior to
#' spatial thinning or modeling steps. It assumes that occurrence files
#' have been created by \code{get_ebird_occurrences()} or equivalent
#' preprocessing steps.
#'
#' @param alpha_code Character(1). Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character(1) or NULL. Optional project root directory.
#'   If NULL, resolved via \code{rENM_project_dir()}.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{alpha_code}: Species alpha code
#'   \item \code{files_processed}: Character vector of processed file paths
#'   \item \code{duplicates_removed}: Named integer vector of counts removed per file
#'   \item \code{output_dir}: Directory containing cleaned occurrence files
#'   \item \code{log_file}: Path to the run log file
#' }
#'
#' @examples
#' \dontrun{
#'   remove_duplicate_occurrences("CASP")
#'
#'   remove_duplicate_occurrences(
#'     "CASP",
#'     project_dir = "/projects/rENM"
#'   )
#' }
#'
#' @family occurrence processing
#' @export
remove_duplicate_occurrences <- function(alpha_code, project_dir = NULL) {

  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  # ---------------------------------------------------------------------------
  # Resolve project directory
  # ---------------------------------------------------------------------------
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
      "Run get_ebird_occurrences() first.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Identify occurrence files
  # ---------------------------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)

  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep <- function(n = 72L) paste(rep("-", n), collapse = "")

  processed_files <- character(0)
  removed_counts  <- integer(0)

  # ---------------------------------------------------------------------------
  # Process each file
  # ---------------------------------------------------------------------------
  for (f in files) {

    dat <- tryCatch(
      utils::read.csv(f, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Failed to read file: ", f, " | ", conditionMessage(e), call. = FALSE)
      }
    )

    if (!all(c("longitude", "latitude") %in% names(dat))) {
      stop("Missing required columns in file: ", f, call. = FALSE)
    }

    n_before <- nrow(dat)

    dat_unique <- dat[!duplicated(dat[, c("longitude", "latitude")]), ]

    n_after   <- nrow(dat_unique)
    n_removed <- n_before - n_after

    utils::write.csv(dat_unique, f, row.names = FALSE)

    processed_files <- c(processed_files, f)
    removed_counts  <- c(removed_counts, n_removed)

    message(sprintf("Processed %s | removed %d duplicates", basename(f), n_removed))
  }

  names(removed_counts) <- basename(processed_files)

  # ---------------------------------------------------------------------------
  # Logging
  # ---------------------------------------------------------------------------
  log_lines <- c(
    .sep(),
    "Processing summary (remove_duplicate_occurrences)",
    paste0("Timestamp: ", .now()),
    paste0("Alpha code: ", alpha_code),
    paste0("Files processed: ", length(processed_files)),
    "Duplicate removal by file:",
    paste0("  - ", names(removed_counts), ": ", removed_counts)
  )

  cat(paste(log_lines, collapse = "\n"), file = log_file, append = TRUE)

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  invisible(list(
    alpha_code = alpha_code,
    files_processed = processed_files,
    duplicates_removed = removed_counts,
    output_dir = tmp_dir,
    log_file = log_file
  ))
}
