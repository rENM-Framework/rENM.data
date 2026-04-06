#' Initialize an rENM run directory structure
#'
#' Create the directory layout required for a single rENM modeling run
#' identified by a four-letter species alpha code. Missing directories
#' are created while existing directories are preserved.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' The base project directory is resolved using:
#'
#' \code{rENM_project_dir()}
#'
#' The run directory is constructed as:
#'
#' \code{<project_dir>/runs/<alpha_code>}
#'
#' \strong{Directory structure}:
#' The following directories are created if they do not already exist:
#'
#' \itemize{
#'   \item \code{_occs/} occurrence input/output files
#'   \item \code{_vars/} cropped predictor rasters
#'   \item \code{Trends/} trend analysis outputs
#'   \item \code{Summaries/} plots, tables, and reports
#'   \item \code{TimeSeries/<year>/model/} model outputs
#'   \item \code{TimeSeries/<year>/occs/} processed occurrences
#'   \item \code{TimeSeries/<year>/vars/} environmental predictors
#' }
#'
#' Additional subdirectories under \code{Trends/} include:
#'
#' \itemize{
#'   \item \code{centroids/}
#'   \item \code{suitability/}
#'   \item \code{variables/}
#' }
#'
#' Time series directories are created for 5-year bins spanning
#' 1980 to 2020.
#'
#' \strong{Outputs}:
#' A standardized log file is created or appended at:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' recording created and existing directories.
#'
#' \strong{Methods}:
#' \itemize{
#'   \item Resolves project directory
#'   \item Constructs run-specific directory paths
#'   \item Creates missing directories recursively
#'   \item Records created and pre-existing paths
#'   \item Appends a processing summary to the run log
#' }
#'
#' \strong{Data requirements}:
#' No input data are required. This function initializes structure only
#' and does not validate or process data files.
#'
#' @param alpha_code Character(1). Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character(1) or NULL. Optional project root directory.
#'   If NULL, resolved via \code{rENM_project_dir()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{alpha_code}: Species alpha code
#'   \item \code{run_dir}: Full path to the run directory
#'   \item \code{created}: Character vector of newly created directories
#'   \item \code{existing}: Character vector of directories that already existed
#'   \item \code{log_file}: Path to the run log file
#' }
#'
#' @examples
#' \dontrun{
#'   Sys.setenv(RENM_PROJECT_DIR = "~/rENM")
#'   run <- set_up_run("CASP")
#'
#'   run <- set_up_run("CASP", project_dir = "/projects/rENM")
#'
#'   run$run_dir
#'   run$created
#' }
#'
#' @family run preparation
#' @export
set_up_run <- function(alpha_code, project_dir = NULL) {

  ## ---- validation ------------------------------------------------------------
  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nzchar(alpha_code))

  ## ---- small helpers ---------------------------------------------------------
  .mkdir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  .now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep_line <- function(n = 72L) paste(rep.int("-", n), collapse = "")
  .catln <- function(...) { cat(paste0(..., "\n")); flush.console() }

  .append_log <- function(title, lines, logfile) {
    section <- c("", .sep_line(), title, paste0("Timestamp: ", .now()), lines, "")
    cat(paste0(section, collapse = "\n"), file = logfile, append = TRUE)
  }

  ## ---- resolve project directory --------------------------------------------
  project_root <- rENM_project_dir(project_dir)

  ## ---- construct run paths ---------------------------------------------------
  runs_root <- file.path(project_root, "runs")
  run_dir <- file.path(runs_root, alpha_code)

  occs_dir <- file.path(run_dir, "_occs")
  vars_dir <- file.path(run_dir, "_vars")
  trends <- file.path(run_dir, "Trends")
  summaries <- file.path(run_dir, "Summaries")
  ts_root <- file.path(run_dir, "TimeSeries")

  years <- seq(1980, 2020, by = 5)

  log_fp <- file.path(run_dir, "_log.txt")

  ## ---- banner ----------------------------------------------------------------
  .catln(.sep_line())
  .catln(sprintf("set_up_run: alpha_code=%s", alpha_code))
  .catln("Project root: ", project_root)
  .catln("Run directory: ", run_dir)
  .catln(.sep_line())

  ## ---- create directory structure -------------------------------------------
  created <- character(0)
  ensured <- character(0)

  for (d in c(runs_root, run_dir, occs_dir, vars_dir, trends, summaries, ts_root)) {
    if (!dir.exists(d)) {
      .mkdir(d)
      created <- c(created, d)
      .catln("  + created: ", d)
    } else {
      ensured <- c(ensured, d)
    }
  }

  ## ---- Trends subdirectories -------------------------------------------------
  for (leaf in c("centroids", "suitability", "variables")) {
    d <- file.path(trends, leaf)

    if (!dir.exists(d)) {
      .mkdir(d)
      created <- c(created, d)
      .catln("  + created: ", d)
    } else {
      ensured <- c(ensured, d)
    }
  }

  ## ---- TimeSeries structure --------------------------------------------------
  for (y in years) {

    for (leaf in c("model", "vars", "occs")) {

      d <- file.path(ts_root, as.character(y), leaf)

      if (!dir.exists(d)) {
        .mkdir(d)
        created <- c(created, d)
        .catln("  + created: ", d)
      } else {
        ensured <- c(ensured, d)
      }
    }
  }

  ## ---- summary ---------------------------------------------------------------
  n_created <- length(created)
  n_existing <- length(ensured)

  .catln(.sep_line())
  .catln(sprintf("Done. Created %d new path(s); %d already existed.", n_created, n_existing))
  .catln("Log file: ", log_fp)
  .catln(.sep_line())

  ## ---- append log ------------------------------------------------------------
  .append_log("Processing summary (set_up_run)", c(
    sprintf("Alpha code:           %s", alpha_code),
    sprintf("Project root:         %s", project_root),
    sprintf("Run root:             %s", run_dir),
    sprintf("Created (count):      %d", n_created),
    if (n_created) "Created paths:" else NULL,
    if (n_created) paste0("  - ", created) else NULL,
    sprintf("Pre-existing (count): %d", n_existing),
    if (n_existing) "Pre-existing paths:" else NULL,
    if (n_existing) paste0("  - ", ensured) else NULL
  ), log_fp)

  invisible(list(
    alpha_code = alpha_code,
    run_dir = run_dir,
    created = created,
    existing = ensured,
    log_file = log_fp
  ))
}
