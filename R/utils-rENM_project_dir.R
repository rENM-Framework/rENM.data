#' Resolve and validate an rENM project directory on the local filesystem
#'
#' Determines the rENM project root directory using a defined precedence
#' order and verifies that it exists, returning a normalized absolute
#' path suitable for use by package functions.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' This helper determines the rENM project root directory and verifies
#' that it exists. It is intended to centralize project-directory
#' configuration for all functions in the package, ensuring
#' CRAN-friendly behavior with no hard-coded home paths, no assumptions
#' about working directory, and clear error messages.
#'
#' The project directory is resolved using the following precedence:
#'
#' \enumerate{
#'   \item The explicit project_dir argument (recommended for scripts/tests)
#'   \item options(rENM.project_dir = "...")
#'   \item The environment variable RENM_PROJECT_DIR
#' }
#'
#' If a candidate value is found, it is expanded (e.g., "~"),
#' normalized to an absolute path, and verified to exist using
#' normalizePath(..., mustWork = TRUE). If no valid directory can be
#' determined, the function errors with actionable setup instructions.
#'
#' \strong{Inputs}
#' The function accepts an optional project directory path or attempts
#' to resolve it from configuration sources including options and
#' environment variables.
#'
#' \strong{Outputs}
#' Returns a normalized absolute path to an existing directory.
#'
#' \strong{Methods}
#' Resolution is performed using a helper that:
#'
#' \itemize{
#'   \item Validates that input is a non-empty scalar character string
#'   \item Expands "~" using path.expand
#'   \item Normalizes and verifies existence using normalizePath with
#'   mustWork = TRUE
#' }
#'
#' If a candidate fails validation or does not exist, it is ignored and
#' the next configuration source is attempted.
#'
#' \strong{Data requirements}
#' The resolved directory must exist on the local filesystem. This
#' function does not validate directory contents (e.g., presence of
#' data/_species.csv); content validation is the responsibility of
#' calling functions.
#'
#' \strong{Configuration}
#' Users can configure the project directory in several ways:
#'
#' \itemize{
#'   \item Session-only: \code{Sys.setenv(RENM_PROJECT_DIR = "~/rENM")}
#'   \item Session-only: \code{options(rENM.project_dir = "~/rENM")}
#'   \item Persistent: add \code{RENM_PROJECT_DIR=~/rENM} to
#'   \code{~/.Renviron}
#' }
#'
#' In package examples and tests, it is best practice to pass
#' project_dir explicitly (e.g., a temporary directory or
#' \code{system.file("extdata", package = "rENM")}).
#'
#' @param project_dir Character. Optional path to the rENM project root.
#' If provided, it takes precedence over options and environment
#' variables.
#'
#' @return
#' Character. A length-1 string representing the normalized absolute
#' path to the existing project directory.
#'
#' Side effects:
#' \itemize{
#'   \item Expands and normalizes filesystem paths
#'   \item Validates existence of the resolved directory
#'   \item Throws an error if no valid directory can be resolved
#' }
#'
#' @examples
#' \dontrun{
#' # Preferred for reproducible scripts:
#' rENM_project_dir("/projects/rENM")
#'
#' # Convenience for interactive use:
#' Sys.setenv(RENM_PROJECT_DIR = "~/rENM")
#' rENM_project_dir()
#'
#' # Or using an R option:
#' options(rENM.project_dir = "~/rENM")
#' rENM_project_dir()
#' }
#'
#' @seealso
#' \link{normalizePath}
#'
#' @export
rENM_project_dir <- function(project_dir = NULL) {

  # Internal helper: accept a candidate path, ensure it is a usable scalar string,
  # expand "~", and convert to a normalized absolute path that must exist.
  #
  # Returning NULL means "no usable path found" and allows the caller to try the
  # next configuration source (argument -> option -> environment variable).
  resolve <- function(p) {
    if (is.null(p)) return(NULL)                       # unset source
    if (!is.character(p) || length(p) != 1) return(NULL)  # must be a scalar string
    if (is.na(p) || !nzchar(p)) return(NULL)           # reject NA/empty strings

    # path.expand() handles "~" reliably across platforms.
    # normalizePath(..., mustWork=TRUE) ensures it exists and returns an absolute path.
    normalizePath(path.expand(p), mustWork = TRUE)
  }

  # 1) Explicit argument (most reproducible; best for scripts/tests)
  p <- resolve(project_dir)
  if (!is.null(p)) return(p)

  # 2) Package option (session-level convenience)
  p <- resolve(getOption("rENM.project_dir"))
  if (!is.null(p)) return(p)

  # 3) Environment variable (convenient and can be made persistent via .Renviron)
  p <- resolve(Sys.getenv("RENM_PROJECT_DIR", unset = NA_character_))
  if (!is.null(p)) return(p)

  # If nothing worked, fail with a helpful, actionable message rather than
  # letting downstream file reads fail with confusing errors.
  stop(
    "rENM project directory not found.\n",
    "Provide `project_dir =` or set one of:\n",
    "  options(rENM.project_dir = '/path/to/rENM')\n",
    "  Sys.setenv(RENM_PROJECT_DIR = '/path/to/rENM')",
    call. = FALSE
  )
}
