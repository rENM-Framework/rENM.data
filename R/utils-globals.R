# utils-globals.R
# Package-level declarations for rENM.data

# Suppress R CMD CHECK notes for variables used in non-standard evaluation.
# Review this list when adding functions that use column names as symbols.
utils::globalVariables(c(
  "x",
  "y"
))

#' @importFrom rENM.core rENM_project_dir get_species_info show_species show_variables
NULL
