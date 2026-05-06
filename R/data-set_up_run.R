#' Initialize an rENM run directory structure
#'
#' Creates the directory layout required for a single rENM modeling run
#' identified by a four-letter species alpha code. Missing directories are
#' created; existing directories are preserved.
#'
#' @details
#' The run directory is constructed as \code{PROJECT_DIR/runs/ALPHA_CODE}
#' and populated with the following structure:
#'
#' \preformatted{
#' runs/<alpha_code>/
#'   _occs/                      occurrence input/output files
#'   _vars/                      cropped predictor rasters
#'   Trends/
#'     centroids/
#'     suitability/
#'     variables/
#'   Summaries/                  plots, tables, and reports
#'   TimeSeries/
#'     <year>/
#'       model/                  model outputs
#'       occs/                   processed occurrences
#'       vars/                   environmental predictors
#' }
#'
#' Time series directories are created for 5-year bins spanning 1980 to 2020.
#' A run log is created or appended at \code{RUN_DIR/_log.txt}.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{alpha_code},
#'   \code{run_dir}, \code{created} (character vector of newly created
#'   directory paths), \code{existing} (character vector of pre-existing
#'   paths), and \code{log_file}.
#'
#' @seealso \code{\link[rENM.core]{rENM_project_dir}}
#'
#' @examples
#' \dontrun{
#' run <- set_up_run("CASP")
#' run <- set_up_run("CASP", project_dir = "/projects/rENM")
#' run$run_dir
#' run$created
#' }
#'
#' @family run preparation
#' @export
set_up_run <- function(alpha_code, project_dir = NULL) {

  ## ---- validation ------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  ## ---- resolve project directory --------------------------------------------
  project_root <- rENM_project_dir(project_dir)

  ## ---- construct run paths --------------------------------------------------
  runs_root <- file.path(project_root, "runs")
  run_dir   <- file.path(runs_root, alpha_code)
  occs_dir  <- file.path(run_dir, "_occs")
  vars_dir  <- file.path(run_dir, "_vars")
  trends    <- file.path(run_dir, "Trends")
  summaries <- file.path(run_dir, "Summaries")
  ts_root   <- file.path(run_dir, "TimeSeries")
  years     <- seq(1980L, 2020L, by = 5L)
  log_fp    <- file.path(run_dir, "_log.txt")

  ## ---- banner ---------------------------------------------------------------
  .catln(.sep_line())
  .catln(sprintf("set_up_run: alpha_code=%s", alpha_code))
  .catln("Project root:  ", project_root)
  .catln("Run directory: ", run_dir)
  .catln(.sep_line())

  ## ---- create directories ---------------------------------------------------
  created <- character(0)
  ensured <- character(0)

  ensure <- function(d) {
    if (!dir.exists(d)) {
      .mkdir(d)
      created <<- c(created, d)
      .catln("  + created: ", d)
    } else {
      ensured <<- c(ensured, d)
    }
  }

  for (d in c(runs_root, run_dir, occs_dir, vars_dir, trends, summaries, ts_root)) {
    ensure(d)
  }

  for (leaf in c("centroids", "suitability", "variables")) {
    ensure(file.path(trends, leaf))
  }

  for (y in years) {
    for (leaf in c("model", "vars", "occs")) {
      ensure(file.path(ts_root, as.character(y), leaf))
    }
  }

  ## ---- summary --------------------------------------------------------------
  n_created  <- length(created)
  n_existing <- length(ensured)

  .catln(.sep_line())
  .catln(sprintf("Done. Created %d new path(s); %d already existed.", n_created, n_existing))
  .catln("Log file: ", log_fp)
  .catln(.sep_line())

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (set_up_run)", c(
    sprintf("Alpha code:           %s", alpha_code),
    sprintf("Project root:         %s", project_root),
    sprintf("Run root:             %s", run_dir),
    sprintf("Created (count):      %d", n_created),
    if (n_created)  c("Created paths:",    paste0("  - ", created))  else NULL,
    sprintf("Pre-existing (count): %d", n_existing),
    if (n_existing) c("Pre-existing paths:", paste0("  - ", ensured)) else NULL
  ))

  invisible(list(
    alpha_code = alpha_code,
    run_dir    = run_dir,
    created    = created,
    existing   = ensured,
    log_file   = log_fp
  ))
}
