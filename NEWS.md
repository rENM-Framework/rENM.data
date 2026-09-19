---

editor_options: 
  markdown: 
    wrap: 72
---

# rENM.data 0.2.1.9000

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
