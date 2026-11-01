#' Percentiles and standardized anomalies
#'
#' * `as_ptile()`: `round(p * 100, digits)`. Tables use `digits = 2`; grids
#'   use `as.integer(as_ptile(p, 0))`.
#' * `standardize()`: the standard normal quantile of `p`, clamped to
#'   +/- `z_max`. `method = "acklam"` is the published value (Acklam's
#'   approximation of `p * 100 / 100`, x 100, rounded): `as = "hundredths"`
#'   returns that integer, `as = "z"` returns it / 100. `method = "qnorm"`
#'   uses `stats::qnorm`, unrounded, with `as = "z"` only.
#'
#' @param p Probabilities in \[0, 1\].
#' @param digits Decimal places of the percentile. No default.
#' @param method `"acklam"` or `"qnorm"`. No default.
#' @param as `"z"` or `"hundredths"`. No default.
#' @param z_max Clamp.
#' @return `as_ptile()`: double. `standardize()`: double, or integer for
#'   `as = "hundredths"`.
#' @export
as_ptile <- function(p, digits) not_implemented()

#' @rdname as_ptile
#' @export
standardize <- function(p, method, as, z_max = 4) not_implemented()
