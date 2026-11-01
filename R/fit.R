#' One day's statistics and fit
#'
#' * `day_stats()`: the statistics a family's fit and smoother read, from one
#'   day's reference matrix. For the zero-inflated gamma: `n` (valid years),
#'   `z` (zero years), `l1` and `lcv` (L-moments of the positive values; `NA`
#'   with too few positive values or a degenerate sample).
#' * `fit_day()`: one day's parameters from `day_stats()` output. Rows with
#'   fewer than `min_valid` valid years are all `NA`.
#'
#' @param ref Reference matrix, n_subjects x n_years (`NA` = missing year).
#' @param stats Output of `day_stats()`.
#' @param family A `normals_family` with `stats` and `fit` slots.
#' @param min_valid Integer, length 1 or `nrow(stats)`. No default.
#' @return `day_stats()`: n x s matrix, columns `family$stats_names`.
#'   `fit_day()`: n x p matrix, columns `family$params`.
#' @export
day_stats <- function(ref, family) not_implemented()

#' @rdname day_stats
#' @export
fit_day <- function(stats, family, min_valid) not_implemented()
