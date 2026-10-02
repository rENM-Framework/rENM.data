#' Finalize occurrence files by moving them from staging to the main directory
#'
#' Copies all CSV files from \code{_occs/tmp/} into \code{_occs/}, overwrites
#' any existing files with matching names, and removes the temporary directory.
#' This is the final cleanup step in the occurrence preprocessing pipeline,
#' signaling that the files are ready for downstream modeling.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{files_moved} (character vector of destination paths),
#'   \code{destination_dir}, and \code{log_file}.
#'
#' @seealso \code{\link{get_ebird_occurrences}},
#'   \code{\link{remove_duplicate_occurrences}},
#'   \code{\link{thin_occurrences}}
#'
#' @importFrom utils flush.console
#'
#' @examples
#' \dontrun{
#' tidy_occurrences("CASP")
#' tidy_occurrences("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @export
tidy_occurrences <- function(alpha_code, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  alpha_code <- toupper(trimws(alpha_code))

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir  <- file.path(project_root, "runs", alpha_code)
  occ_dir  <- file.path(run_dir, "_occs")
  tmp_dir  <- file.path(occ_dir, "tmp")
  log_fp   <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Temporary occurrence directory not found for '", alpha_code, "'.\n",
      "Expected: ", tmp_dir,
      call. = FALSE
    )
  }

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No CSV files found in temporary directory: ", tmp_dir, call. = FALSE)
  }

  ## ---- copy files -----------------------------------------------------------
  moved_files <- character(0)

  for (f in files) {
    dest    <- file.path(occ_dir, basename(f))
    success <- file.copy(f, dest, overwrite = TRUE)

    if (!success) {
      stop("Failed to copy: ", f, " -> ", dest, call. = FALSE)
    }

    moved_files <- c(moved_files, dest)
    message("Moved: ", basename(f))
  }

  ## ---- remove tmp directory -------------------------------------------------
  # unlink() reports failure through its return value and the directory's
  # continued existence, not an error. Inside a Dropbox folder the sync client
  # can briefly hold files written moments earlier, so a single attempt
  # sometimes fails; the directory was then left behind while the log said
  # it had been removed. Retry briefly, then say what actually happened.
  for (attempt in 1:5) {
    unlink(tmp_dir, recursive = TRUE, force = TRUE)
    if (!dir.exists(tmp_dir)) break
    Sys.sleep(1)
  }
  tmp_removed <- !dir.exists(tmp_dir)
  if (tmp_removed) {
    message("Removed temporary directory: ", tmp_dir)
  } else {
    warning("Could not remove temporary directory after 5 attempts: ", tmp_dir,
            ". Nothing reads it; it can be deleted by hand.", call. = FALSE)
  }

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (tidy_occurrences)", c(
    sprintf("Alpha code:               %s", alpha_code),
    sprintf("Files moved:              %d", length(moved_files)),
    "Moved files:",
    paste0("  - ", basename(moved_files)),
    sprintf("Temporary directory %s: %s",
            if (tmp_removed) "removed" else "NOT removed", tmp_dir)
  ))

  invisible(list(
    alpha_code      = alpha_code,
    files_moved     = moved_files,
    destination_dir = occ_dir,
    log_file        = log_fp
  ))
}
