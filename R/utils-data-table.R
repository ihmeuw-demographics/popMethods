# data.table is generally careful to minimize the scope for namespace
# conflicts (i.e., functions with the same name as in other packages);
# a more conservative approach using @importFrom should be careful to
# import any needed data.table special symbols as well, e.g., if you
# run DT[ , .N, by='grp'] in your package, you'll need to add
# @importFrom data.table .N to prevent the NOTE from R CMD check.
# See ?data.table::`special-symbols` for the list of such symbols
# data.table defines; see the 'Importing data.table' vignette for more
# advice (vignette('datatable-importing', 'data.table')).
#
#' @import data.table
NULL

# this is needed for non-standard evaluation in data.table (and some other
# packages). Multiple links suggest using `utils::globalVariables` to remove
# notes when checking the package.
# https://www.r-bloggers.com/no-visible-binding-for-global-variable/
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887
utils::globalVariables(
  c("initial", "value", "parameter", "parameters", "method", "weight",
    "original_draw", "draw", "chain_draw", "chain",
    "year", "year_start", "year_index", "year_index_start", "year_index_end",
    "sex", "sex_index", "sex_index_start", "sex_index_end",
    "age_start", "age_end", "age_index", "age_index_start", "age_index_end",
    "transformed_value", "initial_value", "component")
)
