#' Downsample occurrence records per bin to a maximum count
#'
#' Randomly samples up to \code{record_count} rows from each per-bin
#' occurrence CSV in \code{_occs/tmp/}, overwriting each file in place.
#' Files with fewer rows than \code{record_count} are left unchanged.
#' Reproducible sampling requires setting a seed externally with
#' \code{set.seed()} before calling this function.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param record_count Integer scalar. Maximum records to retain per bin.
#'   Must be positive. Default is 250.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{files_processed}, \code{counts_before} (named integer vector),
#'   \code{counts_after} (named integer vector), \code{output_dir}, and
#'   \code{log_file}.
#'
#' @seealso \code{\link{get_ebird_occurrences}},
#'   \code{\link{remove_duplicate_occurrences}},
#'   \code{\link{thin_occurrences}}, \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' \dontrun{
#' limit_record_count("CASP", record_count = 250)
#' limit_record_count("CASP", record_count = 500, project_dir = "/projects/rENM")
#' }
#'
#' @family occurrence processing
#' @export
limit_record_count <- function(alpha_code, record_count = 250, project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(record_count) || length(record_count) != 1L ||
      !is.finite(record_count) || record_count <= 0) {
    stop("`record_count` must be a positive numeric scalar.", call. = FALSE)
  }

  alpha_code   <- toupper(trimws(alpha_code))
  record_count <- as.integer(record_count)

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir  <- file.path(project_root, "runs", alpha_code)
  tmp_dir  <- file.path(run_dir, "_occs", "tmp")
  log_fp   <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop(
      "Occurrence tmp directory not found for '", alpha_code, "'.\n",
      "Expected: ", tmp_dir,
      call. = FALSE
    )
  }

  ## ---- discover files -------------------------------------------------------
  files <- list.files(tmp_dir, pattern = "^of-\\d+\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No occurrence files found in: ", tmp_dir, call. = FALSE)
  }

  ## ---- process each file ----------------------------------------------------
  .catln(.sep_line())
  .catln("limit_record_count: alpha_code=", alpha_code)
  .catln("Record count limit: ", record_count)
  .catln("Directory:          ", tmp_dir)
  .catln(.sep_line())

  processed_files <- character(0)
  before_counts   <- integer(0)
  after_counts    <- integer(0)

  for (f in files) {
    dat <- tryCatch(
      utils::read.csv(f, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Failed to read: ", f, " | ", conditionMessage(e), call. = FALSE)
      }
    )

    n_before <- nrow(dat)

    if (n_before > record_count) {
      dat <- dat[sample.int(n_before, size = record_count, replace = FALSE), , drop = FALSE]
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

  .catln(.sep_line())
  .catln("Done. Total before: ", total_before, " | Total after: ", total_after)
  .catln(.sep_line())

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (limit_record_count)", c(
    sprintf("Alpha code:         %s", alpha_code),
    sprintf("Record count limit: %d", record_count),
    sprintf("Files processed:    %d", length(processed_files)),
    sprintf("Total before:       %d", total_before),
    sprintf("Total after:        %d", total_after),
    "Per-file counts:",
    paste0("  - ", names(before_counts), ": ", before_counts, " -> ", after_counts)
  ))

  invisible(list(
    alpha_code      = alpha_code,
    files_processed = processed_files,
    counts_before   = before_counts,
    counts_after    = after_counts,
    output_dir      = tmp_dir,
    log_file        = log_fp
  ))
}
