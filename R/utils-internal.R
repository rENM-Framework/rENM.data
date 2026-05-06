# utils-internal.R
# Internal helper functions for rENM.data
#
# Follows the pattern established in rENM.core/R/utils-internal.R.
# Not exported. Functions are available package-wide without qualification.


# Path helpers ----------------------------------------------------------------

#' @noRd
.expand <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)

#' @noRd
.mkdir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)


# Console output --------------------------------------------------------------

#' @noRd
.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

#' @noRd
.sep_line <- function(n = 72L) paste(rep.int("-", n), collapse = "")

#' @noRd
.catln <- function(...) { cat(paste0(..., "\n")); flush.console() }


# Logging ---------------------------------------------------------------------

# Append a titled, timestamped section to a run log file.
# log_fp  -- full path to the log file (created if absent)
# title   -- section heading string
# lines   -- character vector of body lines
#' @noRd
.append_log <- function(log_fp, title, lines) {
  dir.create(dirname(log_fp), recursive = TRUE, showWarnings = FALSE)
  section <- c(
    "",
    .sep_line(),
    title,
    paste0("Timestamp: ", .now()),
    lines,
    ""
  )
  cat(paste0(section, collapse = "\n"), file = log_fp, append = TRUE)
  invisible(log_fp)
}


# Column detection ------------------------------------------------------------

# Normalize a character vector to lowercase alphanumeric only.
#' @noRd
.norm_names <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Find the first column in col_names that matches any of candidates
# after normalization. Returns NA_character_ if none found.
#' @noRd
.find_col <- function(col_names, candidates) {
  if (!length(col_names)) return(NA_character_)
  normed <- .norm_names(col_names)
  for (cand in candidates) {
    idx <- match(.norm_names(cand), normed, nomatch = 0L)
    if (idx > 0L) return(col_names[idx])
  }
  NA_character_
}
