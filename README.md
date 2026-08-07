# rENM.data

![rENM](https://img.shields.io/badge/rENM-framework-blue) ![module](https://img.shields.io/badge/module-data-informational) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20789477.svg)](https://doi.org/10.5281/zenodo.20789477)

**Data assembly and preprocessing for the rENM Framework**

## Overview

`rENM.data` provides standardized workflows for assembling, cleaning, and preparing the input datasets used in rENM analyses. It handles species occurrence data from eBird and environmental predictors from MERRA-2, producing consistent, analysis-ready inputs for downstream modeling.

This package depends on `rENM.core` for project-directory resolution and species metadata access. All functions accept an optional `project_dir` argument; see `?rENM_project_dir` for configuration options.

## Key functions

| Function | Description |
|----|----|
| `set_up_run()` | Initialize the directory structure for a species run |
| `get_ebird_occurrences()` | Read and bin an eBird EBD file into 5-year temporal bins |
| `remove_duplicate_occurrences()` | Remove exact-coordinate duplicate records |
| `thin_occurrences()` | Sequential spatial thinning (Haversine nearest-neighbor) |
| `thin_occurrences2()` | Parallel spatial thinning with optional record cap |
| `limit_record_count()` | Randomly downsample bins to a maximum record count |
| `tidy_occurrences()` | Finalize: move cleaned files from staging to main directory |
| `find_occurrence_extent()` | Derive spatial extent from occurrence data |
| `find_range_extent()` | Derive spatial extent from a USGS GAP range polygon |
| `set_extent()` | Set spatial extent from explicit coordinates |
| `get_merra_variables()` | Crop MERRA-2 predictor rasters to the species extent |

## Installation

``` r
# From GitHub
remotes::install_github("rENM-Framework/rENM.data")

# From a local source directory
remotes::install_local("rENM.data")
```

## Getting started

Set up a project directory first (see `rENM.core`), then run the occurrence preprocessing pipeline in order:

``` r
library(rENM.data)

proj <- "/path/to/your/rENM/project"

# 1. Initialize run directory structure
set_up_run("CASP", project_dir = proj)

# 2. Read eBird EBD file and bin by 5-year period
get_ebird_occurrences("CASP", project_dir = proj)

# 3. Clean occurrence records
remove_duplicate_occurrences("CASP", project_dir = proj)
thin_occurrences("CASP", thin_distance = 10, project_dir = proj)
limit_record_count("CASP", record_count = 250, project_dir = proj)

# 4. Move finalized records to main directory
tidy_occurrences("CASP", project_dir = proj)

# 5. Derive spatial extent from occurrence data
find_occurrence_extent("CASP", bbox_pct = 90, project_dir = proj)

# 6. Crop MERRA-2 predictors to the species extent
get_merra_variables("CASP", project_dir = proj)
```

For interactive work, configure the project directory once per session to avoid passing it to every function:

``` r
options(rENM.project_dir = "/path/to/your/rENM/project")

set_up_run("CASP")
get_ebird_occurrences("CASP")
# ...
```

## Occurrence preprocessing pipeline

```         
get_ebird_occurrences()
        ↓
remove_duplicate_occurrences()
        ↓
thin_occurrences()  or  thin_occurrences2()  (parallel)
        ↓
limit_record_count()          (optional)
        ↓
tidy_occurrences()
```

Each step reads from and writes to `<run_dir>/_occs/tmp/`. `tidy_occurrences()` moves the final files to `<run_dir>/_occs/` and removes the staging directory.

## Spatial extent

Three functions produce the `extent.txt` file required by `get_merra_variables()`. Use whichever fits your workflow:

- `find_occurrence_extent()` — derives extent from the species' occurrence records restricted to the Continental United States (CONUS), using a centered percentile bounding box
- `find_range_extent()` — derives extent from a USGS GAP range polygon with optional symmetric padding
- `set_extent()` — sets extent from explicit bounding box coordinates

## Role in the rENM Framework

`rENM.data` is the second stage in the pipeline:

```         
rENM.core → rENM.data → rENM.model → rENM.analysis → rENM.ai → rENM.reports
```

It produces the run directory structure, occurrence files, and cropped predictor rasters consumed by `rENM.model`.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework** — A modular system for reconstructing and analyzing long-term ecological niche dynamics.
