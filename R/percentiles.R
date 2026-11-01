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
as_ptile <- function(p, digits) {
  if (!is.numeric(p)) stop_arg("p", "must be numeric")
  check_whole(digits, "digits", single = TRUE)
  round(p * 100, digits)
}

#' @rdname as_ptile
#' @export
standardize <- function(p, method, as, z_max = 4) {
  if (!is.numeric(p)) stop_arg("p", "must be numeric")
  check_choice(method, c("acklam", "qnorm"), "method")
  check_choice(as, c("z", "hundredths"), "as")
  if (!is.numeric(z_max) || length(z_max) != 1L || is.na(z_max) || z_max <= 0)
    stop_arg("z_max", "must be a positive number")
  if (method == "qnorm") {
    if (as != "z") stop_arg("as", "must be \"z\" with `method = \"qnorm\"`")
    return(pmin(pmax(stats::qnorm(p), -z_max), z_max))
  }
  # mco-snodas qnorm_grid(p * 100): integer hundredths; / 100 for z.
  h <- cpp_qnorm_grid(as.double(p * 100), as.double(z_max))
  if (as == "hundredths") h else h / 100
}
