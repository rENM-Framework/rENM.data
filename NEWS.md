---

editor_options: 
  markdown: 
    wrap: 72
---

# rENM.data 0.2.0.9000

- Ships a land/water footprint for the MERRA-2 grid as
  `inst/resources/evland-mask.tif`, and `get_merra_variables()` stages it as
  `_vars/<year>/evland.asc` when, and only when, the dataset carries no
  `evland` of its own. `create_ensemble_model()` began masking every bin to
  land in rENM.model 0.2.0, reading that mask from `_vars/<year>/evland.asc`.
  The mask had always been staged incidentally, as one more raster found
  under `data/merra/m2/<year>/`, so the change quietly turned a modeling
  decision into a data requirement. The v0.1.0 extended dataset carries
  `evland` among its thirty m2 variables and is unaffected; the v0.1.0
  example dataset carries seven and does not, so every bin failed with
  "Land mask source not found" against data that was complete for the
  version it was published with.
- The packaged file is a footprint, not data: land cells hold 1 and water is
  `NA`, which is all `raster::mask()` reads. It derives from the `evland`
  published in the v0.1.0 extended dataset, verified identical to it and to
  all nine bins of the author's working dataset, zero differing cells. The
  footprint does not vary by year, so one file serves every bin, and
  reducing the raster to its footprint takes it from 92.8 MB as `.asc` to
  129 KB. Provenance is recorded in `inst/resources/evland-mask-README.md`.
  The dataset stays authoritative where it speaks: a later MERRA release
  carrying a corrected footprint, or a study supplying its own, is used in
  preference to this copy rather than being silently overridden, and the log
  records which source was used.
- The mask is resampled onto a raster already staged for the bin rather than
  cropped to the same bounding box. It is on the same 0.0833 degree spacing
  but its global grid starts elsewhere, so cropping snapped to different
  cell edges and produced a grid five rows taller, which
  `create_ensemble_model()` rejected at `compareRaster()` with "different
  extent".

- Added `clip_occurrences()`, which removes records outside the extent box
  in `_occs/extent.txt` from each bin before thinning and capping. Records
  outside the box have no predictor values, so `sdm` drops them when a model
  is fitted, but until now they still counted toward the per-bin record cap.
  For Grace's Warbler, 126 to 181 of each bin's 250 records reached the
  models. The extent must now be set before thinning; `rENM()` runs
  `find_range_extent()` and then `clip_occurrences()` directly after
  `remove_duplicate_occurrences()`.
- `tidy_occurrences()` now checks that the temporary occurrence directory
  was actually removed. It called `unlink()` once, ignored the result and
  logged the directory as removed; inside Dropbox the sync client can hold
  files written moments earlier, so the directory was occasionally left
  behind (one run in six on 1 October). It now retries up to five times a
  second apart, warns if the directory survives, and logs what happened.
- `find_occurrence_extent()` and `set_extent()` are documented as non-standard.
  `find_range_extent()` is the method the pipeline uses and the only one the
  User Manual now documents. Both alternatives remain exported, since no USGS
  GAP range polygon exists outside CONUS or for taxa GAP does not cover, but
  their help now states two consequences. An extent derived from occurrence
  records moves as those records accumulate, and eBird effort is far from
  uniform across 1980-2024, so the region the models draw their background
  from would partly reflect sampling effort rather than an a priori hypothesis
  about accessible area. And neither writes a buffered range polygon, so
  boundary statistics have no interior or ring to compare and are skipped.

- `find_range_extent()`: the modeled extent is now derived from a true 250 km real-world buffer around the GAP range polygon rather than a symmetric percentage pad of its bounding box. The buffer is applied in EPSG:5070 (USA Contiguous Albers Equal Area Conic), so the margin is a real ground distance in every direction; the previous `pad_pct` approach padded in lon/lat degrees, which varies with latitude and corresponds to no fixed distance. On the pilot species a 2% pad worked out to roughly 14-17 km, against centroid displacements of 59-103 km already observed over the study window.
- `find_range_extent()`: the `pad_pct` argument is replaced by `buffer_km` (default `250`). The default is documented in the function help, including its derivation from Huang, Sauer & Dubayah (2017) and the caveat that it is drawn from permanent resident species.
- `find_range_extent()`: the buffered polygon is now saved to `runs/<ALPHA_CODE>/_occs/range_buffered.gpkg` (in EPSG:5070) instead of being discarded once its bounding box is known, so that boundary/buffer-ring statistics can use the polygon itself. The returned list gains a `buffered_polygon` element giving its path.
- `find_range_extent()`: `extent.txt` header comments now describe the buffer method (`buffer_km`, buffer CRS, buffered polygon filename) in place of the former `pad_pct` line, which would otherwise have described a method no longer in use. The `Upper-left:`/`Lower-right:` lines are unchanged in format and remain WGS84 (EPSG:4326).
- `find_occurrence_extent()`: default `bbox_pct` changed from 99 to 90.
- `find_occurrence_extent()`: occurrence records are now restricted to the Continental United States (CONUS) bounding box (Upper Left: -125.0 Longitude, 49.0 Latitude; Lower Right: -66.5 Longitude, 24.5 Latitude) before the percentile bounding box is computed.

# rENM.data 0.1.0

- Initial release.
- Added `set_up_run()` to initialize the directory structure for a species run.
- Added `get_ebird_occurrences()` to read and bin eBird EBD files into 5-year temporal bins.
- Added `remove_duplicate_occurrences()` to remove exact-coordinate duplicate records.
- Added `thin_occurrences()` for sequential spatial thinning using Haversine nearest-neighbor distance.
- Added `thin_occurrences2()` for parallel spatial thinning with an optional record cap.
- Added `limit_record_count()` to randomly downsample temporal bins to a maximum record count.
- Added `tidy_occurrences()` to finalize cleaned occurrence files by moving them from staging to the main run directory.
- Added `find_occurrence_extent()` to derive spatial extent from occurrence data using a centered percentile bounding box.
- Added `find_range_extent()` to derive spatial extent from a USGS GAP range polygon with optional symmetric padding.
- Added `set_extent()` to set spatial extent from explicit bounding box coordinates.
- Added `get_merra_variables()` to crop MERRA-2 predictor rasters to the species extent.
