# evland-mask.tif

A land/water footprint for the MERRA-2 grid, used to mask every temporal bin
to one modeling domain before fitting.

## What it is

Land cells carry the value 1; water is `NA`. Only the `NA` pattern is used:
`create_ensemble_model()` applies it with `raster::mask()`, which keys on `NA`
and ignores values. Grid 4320 x 1723 at 0.0833 degrees, WGS 84, extent
-180/180 and -60/83.583, matching the MERRA-2 rasters the framework reads.
2,647,830 land cells.

## Where it comes from

Derived from the MERRA-2 variable `evland` as published in the rENM Framework
v0.1.0 extended dataset (<https://doi.org/10.5281/zenodo.20765324>), by
reducing it to its `NA` footprint. It is a footprint, not data: no `evland`
value survives, and nothing here should be read as an evapotranspiration
measurement.

The footprint does not vary by year. The 1980 and 2020 rasters are identical
cell for cell, so one file serves all nine bins. Verified identical to the
`evland` in the author's working dataset across all nine bins, zero differing
cells, matching extent and resolution.

Reducing the published raster to a footprint and compressing it (DEFLATE,
horizontal predictor, 1 byte per cell) takes it from 92.8 MB as `.asc`, or
11.7 MB as `.tif`, to 129 KB.

## Why it ships with the package

Masking to land is a decision the software makes, not a property of anyone's
dataset. `create_ensemble_model()` began requiring it in rENM.model 0.2.0,
which turned a modeling choice into a data requirement: the v0.1.0 example
dataset carries seven MERRA-2 variables and no `evland`, so every bin failed
with "Land mask source not found" against data that was complete for the
version it was published with.

`get_merra_variables()` uses the dataset's own `evland` whenever the dataset
provides one, and falls back to this file when it does not, logging which was
used. The dataset stays authoritative where it speaks: a later MERRA release
carrying a corrected footprint, or a study supplying its own, is not silently
overridden by a mask frozen into a package version.
