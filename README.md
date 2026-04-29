# rENM.data

![rENM](https://img.shields.io/badge/rENM-framework-blue) ![module](https://img.shields.io/badge/module-data-informational)

**Data assembly and preprocessing for the rENM Framework**

## Overview

`rENM.data` provides standardized workflows for assembling, cleaning, and preparing the core datasets used in rENM analyses. It handles species occurrence data and environmental predictors, ensuring consistent, analysis-ready inputs across all runs.

This package focuses exclusively on **data acquisition, filtering, and preprocessing**.

## Role in the rENM Framework

Within the modular rENM ecosystem, `rENM.data`: - Retrieves and processes **species occurrence data** (e.g., eBird) - Assembles **climate variables** (e.g., MERRA-based products) - Defines and manages **spatial extents** - Prepares **run-ready datasets** for downstream modeling

It serves as the bridge between raw data sources and the analytical pipeline.

## Key Functions

-   `get_ebird_occurrences()` — Retrieve occurrence records
-   `tidy_occurrences()` — Standardize occurrence structure
-   `remove_duplicate_occurrences()` — Clean duplicate records
-   `thin_occurrences()` — Reduce spatial sampling bias
-   `limit_record_count()` — Control dataset size
-   `get_merra_variables()` — Assemble climate predictors
-   `set_extent()` — Define spatial bounds
-   `find_occurrence_extent()` / `find_range_extent()` — Derive extents
-   `set_up_run()` — Initialize a complete, analysis-ready run

## Installation

``` r
devtools::install_local("rENM.data")
```

## Example

``` r
library(rENM.data)

# retrieve and clean occurrences
occ <- get_ebird_occurrences("Peca_cassinii")
occ <- tidy_occurrences(occ)
occ <- remove_duplicate_occurrences(occ)

# prepare predictors
vars <- get_merra_variables(year = 2000)

# initialize run
set_up_run("CASP")
```

## Relationship to Other Packages

`rENM.data` provides the standardized inputs required by all downstream components.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework**\
A modular system for reconstructing and analyzing long-term ecological niche dynamics.
