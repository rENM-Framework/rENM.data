# rENM.data 0.1.0

* Initial release.
* Added `set_up_run()` to initialize the directory structure for a species run.
* Added `get_ebird_occurrences()` to read and bin eBird EBD files into
  5-year temporal bins.
* Added `remove_duplicate_occurrences()` to remove exact-coordinate duplicate
  records.
* Added `thin_occurrences()` for sequential spatial thinning using Haversine
  nearest-neighbor distance.
* Added `thin_occurrences2()` for parallel spatial thinning with an optional
  record cap.
* Added `limit_record_count()` to randomly downsample temporal bins to a
  maximum record count.
* Added `tidy_occurrences()` to finalize cleaned occurrence files by moving
  them from staging to the main run directory.
* Added `find_occurrence_extent()` to derive spatial extent from occurrence
  data using a centered percentile bounding box.
* Added `find_range_extent()` to derive spatial extent from a USGS GAP range
  polygon with optional symmetric padding.
* Added `set_extent()` to set spatial extent from explicit bounding box
  coordinates.
* Added `get_merra_variables()` to crop MERRA-2 predictor rasters to the
  species extent.
