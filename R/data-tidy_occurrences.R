#' Finalize occurrence files by moving from temporary to main directory
#'
#' Moves per-bin occurrence CSV files from a temporary staging directory
#' into the main occurrence directory for a species run, overwriting
#' existing files and removing the temporary directory upon completion.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' This function is typically the final cleanup step after:
#'
#' \itemize{
#'   \item occurrence retrieval
#'   \item duplicate removal
#'   \item spatial thinning
#' }
#'
#' The rENM project directory is resolved using
#' \code{\link{rENM_project_dir}}, ensuring that all filesystem operations
#' are portable and do not rely on hard-coded paths such as ~/rENM.
#'
#' \strong{Inputs}
#' Occurrence files are expected in:
#'
#' \code{<project_root>/runs/<alpha_code>/_occs/tmp/}
#'
#' Files must be CSV format with names corresponding to time bins.
#'
#' \strong{Outputs}
#' Files are copied into:
#'
#' \code{<project_root>/runs/<alpha_code>/_occs/}
#'
#' Existing files with the same name are overwritten.
#'
#' \strong{Methods}
#' This function performs a destructive finalization step:
#'
#' \itemize{
#'   \item Files in _occs/tmp/ are copied into _occs/
#'   \item Existing files in _occs/ with the same name are overwritten
#'   \item The _occs/tmp/ directory is deleted after successful transfer
#' }
#'
#' This step signals that occurrence preprocessing is complete and that
#' the resulting files are ready for downstream modeling.
#'
#' \strong{Data requirements}
#' The temporary directory must exist and contain at least one CSV file.
#'
#' Centralizing project-directory resolution via
#' \code{\link{rENM_project_dir}} ensures consistent behavior across all
#' rENM functions and supports reproducibility, portability, and HPC use.
#'
#' @param alpha_code Character. Four-letter species alpha code (e.g.,
#' "CASP").
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
#' @return
#' List. Invisibly returns a structured list containing:
#'
#' \describe{
#'   \item{alpha_code}{Character. Species alpha code.}
#'   \item{files_moved}{Character vector. Full paths to files moved into
#'   the destination directory.}
#'   \item{destination_dir}{Character. Path to the final _occs/
#'   directory.}
#'   \item{log_file}{Character. Path to the _log.txt file updated with
#'   processing details.}
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Copies CSV files from the temporary directory into the main
#'   occurrence directory
#'   \item Overwrites existing files with matching names
#'   \item Deletes the temporary directory after transfer
#'   \item Appends a processing summary to the log file
#' }
#'
#' @examples
#' \dontrun{
#'
#' tidy_occurrences("CASP")
#'
#' tidy_occurrences("CASP", project_dir = "/projects/rENM")
#'
#' }
#'
#' @seealso
#' \link{rENM_project_dir}
#'
#' @export
tidy_occurrences <- function(alpha_code, project_dir = NULL) {

  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  # ---------------------------------------------------------------------------
  # Resolve project directory (CRAN-compliant)
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
      "Temporary occurrence directory not found for alpha_code '", alpha_code, "'.\n",
      "Expected: ", tmp_dir,
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Identify files to move
  # ---------------------------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No files found in temporary directory: ", tmp_dir, call. = FALSE)
  }

  moved_files <- character(0)

  # ---------------------------------------------------------------------------
  # Move files (overwrite existing)
  # ---------------------------------------------------------------------------
  for (f in files) {

    dest <- file.path(occ_dir, basename(f))

    success <- file.copy(f, dest, overwrite = TRUE)

    if (!success) {
      stop("Failed to copy file: ", f, " -> ", dest, call. = FALSE)
    }

    moved_files <- c(moved_files, dest)

    message("Moved: ", basename(f))
  }

  # ---------------------------------------------------------------------------
  # Remove temporary directory
  # ---------------------------------------------------------------------------
  unlink(tmp_dir, recursive = TRUE, force = TRUE)

  message("Removed temporary directory: ", tmp_dir)

  # ---------------------------------------------------------------------------
  # Logging
  # ---------------------------------------------------------------------------
  .now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep <- function(n = 72L) paste(rep("-", n), collapse = "")

  log_lines <- c(
    "",
    .sep(),
    "Processing summary (tidy_occurrences)",
    paste0("Timestamp: ", .now()),
    paste0("Alpha code: ", alpha_code),
    paste0("Files moved: ", length(moved_files)),
    "Moved files:",
    paste0("  - ", basename(moved_files)),
    paste0("Temporary directory removed: ", tmp_dir)
  )

  cat(paste(log_lines, collapse = "\n"), file = log_file, append = TRUE)

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  invisible(list(
    alpha_code = alpha_code,
    files_moved = moved_files,
    destination_dir = occ_dir,
    log_file = log_file
  ))
}
