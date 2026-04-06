#' Read and process an eBird EBD occurrence record file
#'
#' Extracts occurrence records for a single species from an eBird
#' Basic Dataset (EBD) text file, filters valid coordinates, and
#' organizes the data into temporal bins for rENM workflows.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' The EBD text file is located at:
#'
#' \code{<project_dir>/data/ebird/<ebd_dataset>/}
#' \code{<ebd_dataset>.txt}
#'
#' where <ebd_dataset> is obtained from \code{get_species_info()} via
#' the EBD.RECORDS field.
#'
#' \strong{Methods}:
#' The function performs the following steps:
#' \itemize{
#'   \item Reads the EBD text file using streaming input
#'   \item Filters records to valid decimal-degree coordinates
#'   \item Writes a cleaned CSV copy to the _occs directory
#'   \item Parses observation dates into years
#'   \item Assigns records to 5-year bins spanning 1980 to 2024
#'   \item Writes one CSV file per bin to the _occs/tmp directory
#' }
#'
#' Longitude and latitude columns are automatically detected using
#' common aliases. Coordinates are converted to numeric and filtered
#' to valid ranges:
#' \itemize{
#'   \item longitude in (-180, 180)
#'   \item latitude in (-90, 90)
#'   \item excludes (0, 0) coordinates
#' }
#'
#' Duplicate records are retained.
#'
#' \strong{Outputs}:
#' A cleaned CSV copy is written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/<ebd_dataset>.csv}
#'
#' Per-bin occurrence files are written to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_occs/tmp/}
#'
#' with filenames:
#'
#' \itemize{
#'   \item of-<year>.csv
#' }
#'
#' A processing summary is appended to:
#'
#' \code{<project_dir>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' using the section title:
#'
#' \code{Processing summary (get_ebird_occurrences)}
#'
#' \strong{Data requirements}:
#' The EBD file must contain longitude, latitude, and date columns or
#' recognizable aliases. At least one valid coordinate record must be
#' present.
#'
#' @param alpha_code Character. Four-letter species code.
#' @param project_dir Character, NULL. Optional project root directory.
#'   If NULL, resolved using \code{rENM_project_dir()} via options or
#'   environment variables.
#'
#' @return Invisibly returns a List with the following components:
#' \itemize{
#'   \item csv_copy: File path to the cleaned CSV
#'   \item bin_dir: Directory containing per-bin CSV files
#'   \item bin_files: Character vector of per-bin file paths
#'   \item log_file: File path to the run log
#'   \item n_read: Integer number of rows read from the EBD file
#'   \item n_valid: Integer number of valid coordinate rows
#'   \item total_written: Integer total rows written across all bins
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes cleaned CSV and per-bin files to disk
#'   \item Creates directories as needed
#'   \item Appends a processing summary to the run log
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom utils read.delim read.csv write.csv head flush.console
#'
#' @examples
#' \dontrun{
#'   occ <- get_ebird_occurrences("CASP")
#'   utils::head(occ$csv_copy)
#' }
#'
#' @export
get_ebird_occurrences <- function(alpha_code, project_dir = NULL) {

  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nchar(alpha_code) > 0)

  # ---- resolve rENM project directory ---------------------------------------
  project_root <- rENM_project_dir(project_dir)

  # ---- derive ebd_dataset via get_species_info -------------------------------
  .expand <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
  .mkdir  <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  .now    <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  .sep_line <- function(width = 72L) paste(rep.int("-", width), collapse = "")
  .catln  <- function(...) { cat(paste0(..., "\n")); utils::flush.console() }

  species_info <- tryCatch(
    get_species_info(alpha_code, project_dir = project_dir),
    error = function(e) stop("get_species_info('", alpha_code, "') failed: ", conditionMessage(e))
  )

  if (!("EBD.RECORDS" %in% names(species_info))) {
    stop("get_species_info() did not return an 'EBD.RECORDS' column.")
  }

  ebd_dataset <- as.character(species_info$EBD.RECORDS[[1]])

  if (!isTRUE(nchar(ebd_dataset) > 0)) {
    stop("Derived EBD.RECORDS is empty for alpha_code '", alpha_code, "'.")
  }

  # ---- helpers for parsing/columns -------------------------------------------
  .normalize <- function(x) gsub("[^a-z0-9]", "", tolower(x))

  .find_col <- function(all_names, candidates) {
    if (!length(all_names)) return(NA_character_)
    norm_all <- .normalize(all_names)
    for (cand in candidates) {
      idx <- match(.normalize(cand), norm_all, nomatch = 0)
      if (idx > 0) return(all_names[idx])
    }
    NA_character_
  }

  .write_log_section <- function(logfile, title, lines) {
    dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
    section <- c(
      "",
      .sep_line(),
      sprintf("%s", title),
      sprintf("Timestamp: %s", .now()),
      lines
    )
    cat(paste0(section, collapse = "\n"), file = logfile, append = TRUE)
    cat("\n", file = logfile, append = TRUE)
  }

  # ---- paths (derived from project root) ------------------------------------
  in_txt <- .expand(file.path(
    project_root, "data", "ebird", ebd_dataset, paste0(ebd_dataset, ".txt")
  ))

  run_dir  <- .expand(file.path(project_root, "runs", alpha_code))
  occs_dir <- file.path(run_dir, "_occs")
  tmp_dir  <- file.path(occs_dir, "tmp")

  csv_out  <- file.path(occs_dir, paste0(ebd_dataset, ".csv"))
  log_out  <- file.path(run_dir, "_log.txt")

  # ---- prep ------------------------------------------------------------------
  .mkdir(run_dir)
  .mkdir(occs_dir)
  .mkdir(tmp_dir)

  if (!file.exists(in_txt)) stop("Input EBD file not found: ", in_txt)

  # ---- read EBD .txt -> CSV --------------------------------------------------
  .catln("Reading EBD text: ", in_txt)

  ebd <- tryCatch(
    utils::read.delim(
      in_txt,
      header = TRUE,
      sep = "\t",
      quote = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) stop("Failed to read EBD text: ", conditionMessage(e))
  )

  .catln("Rows read: ", nrow(ebd), " | Columns: ", ncol(ebd))

  .catln("Writing CSV copy to: ", csv_out)
  utils::write.csv(ebd, file = csv_out, row.names = FALSE, quote = TRUE)

  # ---- re-read CSV per spec --------------------------------------------------
  .catln("Re-reading CSV: ", csv_out)

  dat <- utils::read.csv(csv_out, stringsAsFactors = FALSE, check.names = FALSE)

  # ---- locate key columns ----------------------------------------------------
  lon_col  <- .find_col(names(dat), c("LONGITUDE","longitude","decimalLongitude"))
  lat_col  <- .find_col(names(dat), c("LATITUDE","latitude","decimalLatitude"))
  date_col <- .find_col(names(dat), c("OBSERVATION.DATE","OBSERVATION DATE","observationDate","DATE","obsDt"))

  if (any(is.na(c(lon_col, lat_col, date_col)))) {
    stop(
      "Missing essential column(s). Found names include: ",
      paste(utils::head(names(dat), 20), collapse = " | ")
    )
  }

  # ---- coordinate filtering --------------------------------------------------
  .catln("Filtering to valid coordinates (duplicates retained) ...")

  dat$.lon <- suppressWarnings(as.numeric(dat[[lon_col]]))
  dat$.lat <- suppressWarnings(as.numeric(dat[[lat_col]]))

  valid <- !is.na(dat$.lon) & !is.na(dat$.lat) &
    dat$.lon >= -180 & dat$.lon <= 180 &
    dat$.lat >= -90  & dat$.lat <= 90 &
    !(dat$.lon == 0 & dat$.lat == 0)

  dat <- dat[valid, , drop = FALSE]

  # ---- year parse ------------------------------------------------------------
  .year_from_date <- function(x, pivot = as.integer(format(Sys.Date(), "%y"))) {

    x <- trimws(as.character(x))
    y <- rep(NA_integer_, length(x))

    has4 <- grepl("\\d{4}", x)

    if (any(has4)) {
      y[has4] <- as.integer(sub(".*?(\\d{4}).*", "\\1", x[has4]))
    }

    need2 <- is.na(y) & grepl("\\d{2}", x)

    if (any(need2)) {
      yy <- suppressWarnings(as.integer(sub(".*?(\\d{2})(?:\\D*)$", "\\1", x[need2])))
      y[need2] <- ifelse(!is.na(yy) & yy <= pivot, 2000L + yy, 1900L + yy)
    }

    y[!is.na(y) & (y < 1800 | y > 2100)] <- NA_integer_

    y
  }

  .bin5 <- function(y) {
    y <- as.integer(y)
    ok <- !is.na(y) & y >= 1980 & y <= 2024
    out <- rep(NA_integer_, length(y))
    out[ok] <- 1980 + 5 * floor((y[ok] - 1980) / 5)
    out
  }

  years <- .year_from_date(dat[[date_col]])

  dat$.bin <- .bin5(years)

  dat <- dat[!is.na(dat$.bin), , drop = FALSE]

  # ---- write per-bin CSVs ----------------------------------------------------
  bins <- sort(unique(dat$.bin))

  written_files <- character(0)
  bin_counts <- integer(0)

  .catln("Writing per-bin CSVs to: ", tmp_dir)

  for (b in bins) {

    sub <- dat[dat$.bin == b, c(".lon", ".lat"), drop = FALSE]

    if (nrow(sub) == 0) next

    out <- data.frame(
      species   = rep(alpha_code, nrow(sub)),
      longitude = sub$.lon,
      latitude  = sub$.lat,
      stringsAsFactors = FALSE
    )

    out_path <- file.path(tmp_dir, sprintf("of-%d.csv", b))

    utils::write.csv(out, file = out_path, row.names = FALSE, quote = TRUE)

    written_files <- c(written_files, out_path)
    bin_counts    <- c(bin_counts, nrow(out))

    .catln("  - of-", b, ".csv  (", nrow(out), " rows)")
  }

  total_written <- sum(bin_counts)

  .catln("Done writing bins. Files: ", length(written_files), " | Rows total: ", total_written)

  # ---- logging ---------------------------------------------------------------
  bins_block <- if (length(written_files)) {
    c(
      "  Bin files:",
      sprintf("    - of-%d.csv (%d rows)", bins, bin_counts)
    )
  } else {
    "  Bin files: (none)"
  }

  details <- c(
    sprintf("Alpha code:          %s", alpha_code),
    sprintf("EBD dataset:         %s", ebd_dataset),
    sprintf("Input TXT:           %s", in_txt),
    sprintf("CSV copy:            %s", csv_out),
    sprintf("Output tmp dir:      %s", tmp_dir),
    sprintf("Rows read (TXT):     %d", nrow(ebd)),
    sprintf("Rows valid coords:   %d", nrow(dat)),
    "Rows after de-dup:   (duplicates retained by request)",
    "Binning range:       1980-2024 (5-year spans)",
    sprintf("Bins written (N):    %d", length(written_files)),
    sprintf("Rows written (sum):  %d", total_written),
    bins_block
  )

  .write_log_section(
    logfile = log_out,
    title   = "Processing summary (get_ebird_occurrences)",
    lines   = details
  )

  .catln("Log updated for step: Processing summary (get_ebird_occurrences)")
  .catln("Log file: ", log_out)
  .catln(.sep_line())

  invisible(list(
    csv_copy      = csv_out,
    bin_dir       = tmp_dir,
    bin_files     = written_files,
    log_file      = log_out,
    n_read        = nrow(ebd),
    n_valid       = nrow(dat),
    total_written = total_written
  ))
}
