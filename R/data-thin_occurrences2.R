#' Apply spatial thinning to occurrence records (parallel)
#'
#' Applies spatial thinning and optional record capping to per-bin occurrence
#' CSV files using parallel execution via the \pkg{future} and
#' \pkg{future.apply} packages. Files are overwritten in place. This is the
#' parallel counterpart to \code{\link{thin_occurrences}} and is suited to
#' larger datasets where per-bin files are numerous or large.
#'
#' @details
#' \strong{Spatial thinning}
#'
#' If \code{radius > 0}, a greedy Haversine nearest-neighbor filter is applied
#' per file: the first point is always retained; each subsequent point is kept
#' only if it is at least \code{radius} km from all previously retained points.
#' Earth radius is assumed to be 6371 km. Set \code{radius = 0} to skip
#' spatial thinning.
#'
#' \strong{Record cap}
#'
#' If \code{records > 0}, the thinned dataset is further randomly downsampled
#' to at most \code{records} rows per file. Set \code{records = 0} to keep all
#' rows surviving thinning.
#'
#' \strong{Parallel execution}
#'
#' Files are processed independently on separate workers using a
#' \code{future::multisession} plan. If \code{workers} is \code{NULL}, all
#' available logical cores are used. Both \pkg{future} and \pkg{future.apply}
#' must be installed; the function stops with an informative error if they are
#' not.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code
#'   (e.g., \code{"CASP"}).
#' @param radius Numeric scalar. Minimum separation distance in kilometers.
#'   Use \code{0} to disable spatial thinning. Default is 1.
#' @param records Integer scalar. Maximum records to retain per file after
#'   thinning. Use \code{0} to keep all surviving records. Default is 250.
#' @param workers Integer scalar or \code{NULL}. Number of parallel workers.
#'   If \code{NULL} (default), all available logical cores are used.
#' @param project_dir Character. Path to the rENM project root. If \code{NULL}
#'   (default), resolved via \code{\link[rENM.core]{rENM_project_dir}}
#'   (argument, \code{rENM.project_dir} option, \code{RENM_PROJECT_DIR}
#'   environment variable).
#'
#' @return Invisibly returns a data frame with one row per processed file and
#'   columns \code{file} (file name), \code{before} (rows before thinning),
#'   \code{kept_after_radius} (rows after spatial filter, or \code{NA} if
#'   disabled), \code{kept_after_records} (rows after record cap, or
#'   \code{NA} if not applied), and \code{after} (rows written to disk).
#'
#' @seealso \code{\link{thin_occurrences}},
#'   \code{\link{remove_duplicate_occurrences}},
#'   \code{\link{tidy_occurrences}}
#'
#' @importFrom utils read.csv write.csv
#' @importFrom parallel detectCores
#'
#' @examples
#' \dontrun{
#' out <- thin_occurrences2("CASP", radius = 10, records = 500)
#' out
#' thin_occurrences2("CASP", radius = 10, records = 500,
#'                   project_dir = "/projects/rENM")
#' }
#'
#' @export
thin_occurrences2 <- function(alpha_code,
                               radius     = 1,
                               records    = 250,
                               workers    = NULL,
                               project_dir = NULL) {

  ## ---- validate -------------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || radius < 0) {
    stop("`radius` must be a single non-negative number (km).", call. = FALSE)
  }
  if (!is.numeric(records) || length(records) != 1L || is.na(records) || records < 0) {
    stop("`records` must be a single non-negative number. Use 0 to keep all rows.", call. = FALSE)
  }

  records <- as.integer(records)

  ## ---- check parallel dependencies ------------------------------------------
  if (!requireNamespace("future",       quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop(
      "Packages 'future' and 'future.apply' are required for parallel thinning.\n",
      "Install them with: install.packages(c('future', 'future.apply'))",
      call. = FALSE
    )
  }

  ## ---- paths ----------------------------------------------------------------
  project_root <- rENM_project_dir(project_dir)
  run_dir <- file.path(project_root, "runs", alpha_code)
  tmp_dir <- file.path(run_dir, "_occs", "tmp")
  log_fp  <- file.path(run_dir, "_log.txt")

  if (!dir.exists(tmp_dir)) {
    stop("Occurrence tmp directory not found: ", tmp_dir, call. = FALSE)
  }

  files <- list.files(tmp_dir, pattern = "^of-\\d{4}\\.csv$", full.names = TRUE)

  ## ---- early exit if no files -----------------------------------------------
  if (!length(files)) {
    .catln("No of-<year>.csv files found in: ", tmp_dir)
    .append_log(log_fp, "Processing summary (thin_occurrences2 parallel)", c(
      sprintf("Alpha code:       %s", alpha_code),
      sprintf("Target directory: %s", tmp_dir),
      "Files discovered:  0",
      sprintf("Radius (km):      %.3f (0 = no spatial thinning)", radius),
      sprintf("Records cap:      %d (0 = keep all)", records),
      "Action:            No files to process"
    ))
    return(invisible(data.frame(
      file               = character(0),
      before             = integer(0),
      kept_after_radius  = integer(0),
      kept_after_records = integer(0),
      after              = integer(0),
      stringsAsFactors   = FALSE
    )))
  }

  ## ---- parallel setup -------------------------------------------------------
  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores(logical = TRUE))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  .catln(.sep_line())
  .catln(sprintf("thin_occurrences2: alpha_code=%s", alpha_code))
  .catln(sprintf("Files: %d | Workers: %d", length(files), workers))
  .catln(sprintf("Radius: %.3f km | Records cap: %d", radius, records))
  .catln(.sep_line())

  ## ---- Haversine filter (defined here for worker export) --------------------
  .haversine_km <- function(lon1, lat1, lon2, lat2) {
    rad  <- pi / 180
    dlat <- (lat2 - lat1) * rad
    dlon <- (lon2 - lon1) * rad
    a    <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
    2 * 6371 * asin(pmin(1, sqrt(a)))
  }

  .thin_by_radius <- function(lon, lat, radius_km) {
    n <- length(lon)
    if (n == 0L || radius_km <= 0) return(seq_len(n))
    keep_idx <- integer(0)
    for (i in seq_len(n)) {
      if (!length(keep_idx)) {
        keep_idx <- i
      } else {
        d <- .haversine_km(lon[i], lat[i], lon[keep_idx], lat[keep_idx])
        if (all(d >= radius_km)) keep_idx <- c(keep_idx, i)
      }
    }
    keep_idx
  }

  ## ---- parallel processing --------------------------------------------------
  res_list <- future.apply::future_lapply(
    files,
    future.seed = TRUE,
    FUN = function(f) {

      df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)

      required <- c("species", "longitude", "latitude")
      missing  <- setdiff(required, names(df))
      if (length(missing)) {
        stop("Missing columns in ", basename(f), ": ", paste(missing, collapse = ", "))
      }

      n_before          <- nrow(df)
      kept_after_radius <- NA_integer_

      if (radius > 0 && n_before > 1L) {
        idx_keep          <- .thin_by_radius(df$longitude, df$latitude, radius)
        df                <- df[idx_keep, , drop = FALSE]
        kept_after_radius <- nrow(df)
      }

      kept_after_records <- NA_integer_

      if (records > 0L && nrow(df) > records) {
        df                 <- df[sample(seq_len(nrow(df)), records), , drop = FALSE]
        kept_after_records <- nrow(df)
      }

      utils::write.csv(df, f, row.names = FALSE)

      list(
        file               = basename(f),
        before             = n_before,
        kept_after_radius  = kept_after_radius,
        kept_after_records = kept_after_records,
        after              = nrow(df)
      )
    }
  )

  ## ---- assemble result data frame -------------------------------------------
  out_df <- do.call(rbind, lapply(res_list, as.data.frame, stringsAsFactors = FALSE))
  rownames(out_df) <- NULL

  ## ---- log ------------------------------------------------------------------
  .append_log(log_fp, "Processing summary (thin_occurrences2 parallel)", c(
    sprintf("Alpha code:       %s", alpha_code),
    sprintf("Target directory: %s", tmp_dir),
    sprintf("Files processed:  %d", nrow(out_df)),
    sprintf("Radius (km):      %.3f (0 = no spatial thinning)", radius),
    sprintf("Records cap:      %d (0 = keep all)", records),
    sprintf("Workers:          %d", workers),
    sprintf("Total before:     %d", sum(out_df$before)),
    sprintf("Total after:      %d", sum(out_df$after))
  ))

  .catln("Done. Log updated: ", log_fp)
  .catln(.sep_line())

  invisible(out_df)
}
