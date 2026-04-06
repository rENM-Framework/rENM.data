#' Look up species metadata by four-letter alpha code
#'
#' Retrieves a single species record from the rENM species metadata
#' table stored within the project directory using a standardized
#' alpha code lookup.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' The project directory is resolved using
#' \code{\link{rENM_project_dir}}, which centralizes configuration for
#' all filesystem access in the rENM package and avoids hard-coded
#' paths such as ~/rENM. This design improves portability,
#' reproducibility, and CRAN compliance.
#'
#' \strong{Inputs}
#' The species table is expected at:
#'
#' \code{<rENM_project_dir>/data/_species.csv}
#'
#' The lookup uses a four-letter alpha code (e.g., eBird or USGS banding
#' codes) and is case-insensitive. The function tolerates variation in
#' column names (e.g., COMMON.NAME vs Common Name, SCIENTIFIC.NAME vs
#' SciName, ALPHA.CODE vs SPECIESCODE).
#'
#' The species table must contain a column representing the four-letter
#' alpha code. The function attempts to match any of the following
#' (case-insensitive, punctuation ignored):
#'
#' \itemize{
#'   \item ALPHA.CODE
#'   \item SPECIESCODE
#'   \item BANDINGCODE
#'   \item ALPHACODE
#'   \item CODE
#' }
#'
#' \strong{Outputs}
#' Returned column names are standardized regardless of the source file.
#'
#' The returned data frame always contains these columns:
#'
#' \itemize{
#'   \item COMMON.NAME
#'   \item SCIENTIFIC.NAME
#'   \item EBD.RECORDS
#'   \item EBD.RANGE
#'   \item GAP.RANGE
#' }
#'
#' \strong{Methods}
#' Column-name normalization is performed by converting names to
#' uppercase and removing punctuation so that variations such as
#' COMMON.NAME, Common Name, and common_name resolve to the same
#' comparison key.
#'
#' \strong{Data requirements}
#' The species table must exist, be readable, and contain at least one
#' row with the required columns and a valid alpha code field.
#'
#' @param alpha_code Character. Four-letter species alpha code (e.g.,
#' "CASP"). Matching is case-insensitive.
#'
#' @param project_dir Character. Optional rENM project root directory.
#' If NULL, the directory is resolved using
#' \code{\link{rENM_project_dir}} using the following precedence:
#'
#' \enumerate{
#'   \item Explicit project_dir argument
#'   \item options(rENM.project_dir)
#'   \item RENM_PROJECT_DIR environment variable
#' }
#'
#' Supplying project_dir explicitly is recommended for reproducible
#' scripts and for package tests.
#'
#' @return
#' Data frame. Returns a one-row data frame containing standardized
#' species metadata:
#'
#' \describe{
#'   \item{COMMON.NAME}{Character. Common (English) species name.}
#'   \item{SCIENTIFIC.NAME}{Character. Binomial scientific name.}
#'   \item{EBD.RECORDS}{Numeric. Approximate number of eBird occurrence
#'   records.}
#'   \item{EBD.RANGE}{Character. Range description from eBird metadata.}
#'   \item{GAP.RANGE}{Character. Distribution metadata from USGS GAP
#'   analysis.}
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Reads the species metadata table from disk
#'   \item Validates column structure and content
#'   \item Standardizes output column names
#' }
#'
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#'
#' # Standard lookup
#' info <- get_species_info("CASP")
#'
#' # Reproducible script with explicit project directory
#' info <- get_species_info("CASP", project_dir = "/projects/rENM")
#'
#' }
#'
#' @export
get_species_info <- function(alpha_code, project_dir = NULL) {

  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character(1).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Resolve rENM project directory
  # ---------------------------------------------------------------------------
  # Centralized configuration ensures no hard-coded filesystem paths
  project_root <- rENM_project_dir(project_dir)

  csv_path <- file.path(project_root, "data", "_species.csv")

  if (!file.exists(csv_path)) {
    stop("Species CSV not found at: ", csv_path, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Read species table
  # ---------------------------------------------------------------------------
  df <- tryCatch(
    utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop(
        "Failed to read CSV at ", csv_path, ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!nrow(df)) {
    stop("The species table is present but empty: ", csv_path, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Column-name normalization helper
  # ---------------------------------------------------------------------------
  # Converts names to uppercase and strips punctuation so that
  # "COMMON.NAME", "Common Name", and "common_name" all normalize
  # to the same comparison string.
  norm <- function(x) gsub("[^A-Z0-9]", "", toupper(x))

  col_norm <- norm(names(df))

  # ---------------------------------------------------------------------------
  # Identify alpha-code column
  # ---------------------------------------------------------------------------
  alpha_candidates <- c("ALPHACODE", "ALPHA", "SPECIESCODE", "SPSC", "BANDINGCODE")

  alpha_col_idx <- match(TRUE, col_norm %in% alpha_candidates, nomatch = 0L)

  if (alpha_col_idx == 0L) {
    stop(
      "Could not locate the alpha-code column. Looked for one of: ",
      paste(alpha_candidates, collapse = ", "),
      ". Available columns are: ", paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }

  alpha_col <- names(df)[alpha_col_idx]

  # ---------------------------------------------------------------------------
  # Required output columns
  # ---------------------------------------------------------------------------
  wanted      <- c("COMMON.NAME", "SCIENTIFIC.NAME", "EBD.RECORDS", "EBD.RANGE", "GAP.RANGE")
  wanted_norm <- norm(wanted)

  match_idx <- match(wanted_norm, col_norm, nomatch = 0L)

  if (any(match_idx == 0L)) {
    missing <- wanted[match_idx == 0L]
    stop(
      "Missing expected column(s) in species table: ",
      paste(missing, collapse = ", "),
      ". Available columns are: ", paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }

  actual_wanted <- names(df)[match_idx]

  # ---------------------------------------------------------------------------
  # Filter by alpha code
  # ---------------------------------------------------------------------------
  key <- toupper(trimws(alpha_code))

  matches <- df[toupper(trimws(df[[alpha_col]])) == key, , drop = FALSE]

  if (nrow(matches) == 0L) {
    stop("No species found for alpha code '", alpha_code, "'.", call. = FALSE)
  }

  if (nrow(matches) > 1L) {
    stop(
      "Multiple species matched alpha code '", alpha_code, "'. ",
      "Please ensure codes are unique in the species table.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Return standardized metadata
  # ---------------------------------------------------------------------------
  out <- matches[, actual_wanted, drop = FALSE]
  names(out) <- wanted
  rownames(out) <- NULL

  out
}
