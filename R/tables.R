#' Normals and percentile tables
#'
#' `normals_table()` builds every flavor for every date of one normals year,
#' for subjects held as rows of reference matrices (stations, zones, bands).
#' The smoother is fitted at every date except Feb 29 and evaluated at every
#' date.
#'
#' It returns `list(normals, smooth)`. `smooth` is the `normals_smooth` fit
#' (`NULL` without the smoothed flavor). `normals` is a tibble, one row per
#' flavor x date x subject, ordered by flavor (in `flavors` order), date and
#' subject (input order), with columns `flavor`, `date`, `day`, the subject
#' columns, `distribution`, `n_valid`, `ref_first_wy`, `ref_last_wy`, `mean`,
#' the family's parameters, `fit_status` and the quantile columns.
#'
#' `percentiles_table()` places observed values in those normals. It returns a
#' tibble ordered by flavor, date and subject: `flavor`, `date`, the subject
#' columns, `value`, `ptile` and `standardized` (Acklam, to 0.01). Observations
#' of unknown subjects, or on dates without normals, are dropped.
#'
#' @param ref `function(date)` returning that date's reference matrix (rows in
#'   `subjects` order, year column names), or a list of matrices named by
#'   `format(date)`.
#' @param dates The normals year's dates ([year_dates()]).
#' @param subjects Data frame of subject key columns, one row per `ref` row.
#' @param family The family for the fitted and smoothed flavors.
#' @param year_type `"wy"` or `"cy"`.
#' @param min_valid Integer, or `function(ref)` returning one. No default.
#' @param fallback Named character `c(fitted = ., smoothed = .)`, each
#'   `"empirical"` or `"none"`. No default.
#' @param smoothed_mask `"none"` or `"min_valid"`. No default.
#' @param flavors Flavors to build.
#' @param probs Quantile probabilities.
#' @param cores Forked workers for the smoother.
#' @param normals A `normals_table()$normals`.
#' @param obs Data frame of subject key columns, `date` and `value`.
#' @param ptile_digits Decimal places of `ptile`. No default.
#' @export
normals_table <- function(ref, dates, subjects, family, year_type, min_valid,
                          fallback, smoothed_mask,
                          flavors = names(NORMALS_FLAVORS),
                          probs = NORMALS_QUANTILES, cores = 1L) not_implemented()

#' @rdname normals_table
#' @export
percentiles_table <- function(normals, ref, subjects, obs, family, min_valid,
                              fallback, ptile_digits,
                              flavors = names(NORMALS_FLAVORS)) not_implemented()
