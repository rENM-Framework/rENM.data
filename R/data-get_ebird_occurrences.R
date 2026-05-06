#' Read and bin an eBird EBD occurrence file by species
#'
#' Reads a species' eBird Basic Dataset (EBD) text file, filters records to
#' valid coordinates, writes a cleaned CSV copy, and organizes records into
#' 5-year temporal bins for use in rENM workflows.
#'
#' @details
#' \strong{Input location}
#'
#' The EBD text file is located using the \code{EBD.RECORDS} field returned
#' by \code{\link[rENM.core]{get_species_info}}:
#'
#' \code{PROJECT_DIR/data/ebird/EBD_DATASET/EBD_DATASET.txt}
#'
#' \strong{Coordinate filtering}
#'
#' Longitude and latitude columns are detected automatically using common
#' aliases. Records are retained only if coordinates are numeric, within valid
#' ranges (longitude -180 to 180, latitude -90 to 90),
#' and not exactly (0, 0). Duplicate records are retained at this stage.
#'
#' \strong{Temporal binning}
#'
#' Observation dates are parsed to years and assigned to 5-year bins spanning
#' 1980 to 2024 (1980, 1985, ..., 2020). Records outside this range are
#' dropped. One CSV file per bin is written to \code{_occs/tmp/} with the
#' naming pattern \code{of-YEAR.csv}, containing \code{species},
#' \code{longitude}, and \code{latitude} columns.
#'
#' \strong{Outputs}
#'
#' A cleaned copy of the full EBD file is written to
#' \code{RUN_DIR/_occs/EBD_DATASET.csv}. Per-bin files go to
#' \code{RUN_DIR/_occs/tmp/}. A processing summary is appended to
#' \code{RUN_DIR/_log.txt}.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a list with elements \code{csv_copy} (path to
#'   cleaned CSV), \code{bin_dir} (tmp directory), \code{bin_files} (character
#'   vector of per-bin file paths), \code{log_file}, \code{n_read} (rows read
#'   from EBD), \code{n_valid} (rows with valid coordinates), and
#'   \code{total_written} (rows written across all bins).
#'
#' @seealso \code{\link[rENM.core]{rENM_project_dir}},
#'   \code{\link[rENM.core]{get_species_info}},
#'   \code{\link{remove_duplicate_occurrences}},
#'   \code{\link{thin_occurrences}}, \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.delim read.csv write.csv head flush.console
#'
#' @examples
#' \dontrun{
#' occ <- get_ebird_occurrences("CASP")
#' occ <- get_ebird_occurrences("CASP", project_dir = "/projects/rENM")
#' }
#'
#' @export
get_ebird_occurrences <- function(alpha_code, project_dir = NULL) {

  ## ---- validate --------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  ## ---- resolve project directory --------------------------------------------
  project_root <- rENM_project_dir(project_dir)

  ## ---- look up EBD dataset name via species metadata ------------------------
  species_info <- tryCatch(
    get_species_info(alpha_code, project_dir = project_dir),
    error = function(e) {
      stop("get_species_info('", alpha_code, "') failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (!("EBD.RECORDS" %in% names(species_info))) {
    stop("get_species_info() did not return an 'EBD.RECORDS' column.", call. = FALSE)
  }

  ebd_dataset <- trimws(as.character(species_info$EBD.RECORDS[[1L]]))
  if (!nzchar(ebd_dataset)) {
    stop("EBD.RECORDS is empty for alpha_code '", alpha_code, "'.", call. = FALSE)
  }

  ## ---- paths ----------------------------------------------------------------
  in_txt   <- .expand(file.path(project_root, "data", "ebird", ebd_dataset,
                                paste0(ebd_dataset, ".txt")))
  run_dir  <- .expand(file.path(project_root, "runs", alpha_code))
  occs_dir <- file.path(run_dir, "_occs")
  tmp_dir  <- file.path(occs_dir, "tmp")
  csv_out  <- file.path(occs_dir, paste0(ebd_dataset, ".csv"))
  log_fp   <- file.path(run_dir, "_log.txt")

  .mkdir(run_dir)
  .mkdir(occs_dir)
  .mkdir(tmp_dir)

  if (!file.exists(in_txt)) {
    stop("Input EBD file not found: ", in_txt, call. = FALSE)
  }

  ## ---- read EBD .txt --------------------------------------------------------
  .catln(.sep_line())
  .catln("get_ebird_occurrences: alpha_code=", alpha_code)
  .catln("Reading EBD text: ", in_txt)

  ebd <- tryCatch(
    utils::read.delim(in_txt, header = TRUE, sep = "\t", quote = "",
                      check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      stop("Failed to read EBD text: ", conditionMessage(e), call. = FALSE)
    }
  )

  .catln("Rows read: ", nrow(ebd), " | Columns: ", ncol(ebd))
  .catln("Writing CSV copy to: ", csv_out)
  utils::write.csv(ebd, file = csv_out, row.names = FALSE, quote = TRUE)

  ## ---- re-read CSV and locate key columns -----------------------------------
  dat <- utils::read.csv(csv_out, stringsAsFactors = FALSE, check.names = FALSE)

  lon_col  <- .find_col(names(dat), c("LONGITUDE", "longitude", "decimalLongitude"))
  lat_col  <- .find_col(names(dat), c("LATITUDE",  "latitude",  "decimalLatitude"))
  date_col <- .find_col(names(dat), c("OBSERVATION.DATE", "OBSERVATION DATE",
                                      "observationDate", "DATE", "obsDt"))

  if (any(is.na(c(lon_col, lat_col, date_col)))) {
    stop(
      "Missing essential column(s). Found: ",
      paste(utils::head(names(dat), 20L), collapse = " | "),
      call. = FALSE
    )
  }

  ## ---- coordinate filtering -------------------------------------------------
  .catln("Filtering to valid coordinates ...")

  dat$.lon <- suppressWarnings(as.numeric(dat[[lon_col]]))
  dat$.lat <- suppressWarnings(as.numeric(dat[[lat_col]]))

  valid <- !is.na(dat$.lon) & !is.na(dat$.lat) &
    dat$.lon >= -180 & dat$.lon <= 180 &
    dat$.lat >=  -90 & dat$.lat <=  90 &
    !(dat$.lon == 0 & dat$.lat == 0)

  dat <- dat[valid, , drop = FALSE]
  .catln("Valid coordinate rows: ", nrow(dat))

  ## ---- year parsing and 5-year binning --------------------------------------
  .year_from_date <- function(x) {
    x  <- trimws(as.character(x))
    y  <- rep(NA_integer_, length(x))
    h4 <- grepl("\\d{4}", x)
    if (any(h4)) y[h4] <- as.integer(sub(".*?(\\d{4}).*", "\\1", x[h4]))
    pivot <- as.integer(format(Sys.Date(), "%y"))
    h2 <- is.na(y) & grepl("\\d{2}", x)
    if (any(h2)) {
      yy    <- suppressWarnings(as.integer(sub(".*?(\\d{2})(?:\\D*)$", "\\1", x[h2])))
      y[h2] <- ifelse(!is.na(yy) & yy <= pivot, 2000L + yy, 1900L + yy)
    }
    y[!is.na(y) & (y < 1800L | y > 2100L)] <- NA_integer_
    y
  }

  .bin5 <- function(y) {
    y   <- as.integer(y)
    ok  <- !is.na(y) & y >= 1980L & y <= 2024L
    out <- rep(NA_integer_, length(y))
    out[ok] <- 1980L + 5L * floor((y[ok] - 1980L) / 5L)
    out
  }

  dat$.bin <- .bin5(.year_from_date(dat[[date_col]]))
  dat      <- dat[!is.na(dat$.bin), , drop = FALSE]

  ## ---- write per-bin CSVs ---------------------------------------------------
  bins          <- sort(unique(dat$.bin))
  written_files <- character(0)
  bin_counts    <- integer(0)

  .catln("Writing per-bin CSVs to: ", tmp_dir)

  for (b in bins) {
    sub <- dat[dat$.bin == b, c(".lon", ".lat"), drop = FALSE]
    if (!nrow(sub)) next

    out <- data.frame(
      species   = alpha_code,
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
  .catln("Done. Files: ", length(written_files), " | Total rows: ", total_written)
  .catln(.sep_line())

  ## ---- log ------------------------------------------------------------------
  bins_block <- if (length(written_files)) {
    c("Bin files:", sprintf("  - of-%d.csv (%d rows)", bins, bin_counts))
  } else {
    "Bin files: (none)"
  }

  .append_log(log_fp, "Processing summary (get_ebird_occurrences)", c(
    sprintf("Alpha code:          %s", alpha_code),
    sprintf("EBD dataset:         %s", ebd_dataset),
    sprintf("Input TXT:           %s", in_txt),
    sprintf("CSV copy:            %s", csv_out),
    sprintf("Output tmp dir:      %s", tmp_dir),
    sprintf("Rows read (TXT):     %d", nrow(ebd)),
    sprintf("Rows valid coords:   %d", nrow(dat)),
    "Binning range:       1980-2024 (5-year spans)",
    sprintf("Bins written (N):    %d", length(written_files)),
    sprintf("Rows written (sum):  %d", total_written),
    bins_block
  ))

  invisible(list(
    csv_copy      = csv_out,
    bin_dir       = tmp_dir,
    bin_files     = written_files,
    log_file      = log_fp,
    n_read        = nrow(ebd),
    n_valid       = nrow(dat),
    total_written = total_written
  ))
}
