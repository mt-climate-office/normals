#' Distribution families
#'
#' A family is a `normals_family`: a named list of slots that the core
#' dispatches through (see `docs/API.md`, "The family plug-in interface").
#'
#' * `family_zig()`: the zero-inflated gamma (`"zero_inflated_gamma"`),
#'   parameters `shape`, `scale`, `p_zero`; fitted per day by L-moments
#'   (Hosking's PELGAM) and smoothed with cyclic P-splines on logit(p_zero),
#'   log(l1) and log(L-CV).
#' * `family_empirical()`: the zero-inflated type-7 ECDF of the reference
#'   values (`"empirical"`); no parameters and no smoother.
#' * `zig_control()`: the zero-inflated gamma's constants, today's values in
#'   both pipelines.
#'
#' @param control A [zig_control()] list.
#' @param min_pos_fit Fewer positive values than this and a day has no
#'   L-moments (shape and scale `NA`; the fallback signal).
#' @param k Number of cyclic basis functions.
#' @param lambda Smoothing parameters `c(nu, l1, lcv)`.
#' @param lcv_range Clamp on the smoothed L-CV.
#' @param min_pos Fewer positive values on record and only p_zero is smoothed.
#' @param min_days Fewer days with L-moments and l1 and L-CV are constant.
#' @param maxit,tol Iteration limit and convergence tolerance of the p_zero
#'   fit.
#' @return A `normals_family`, or for `zig_control()` a list.
#' @name families
NULL

#' @rdname families
#' @export
family_zig <- function(control = zig_control()) not_implemented()

#' @rdname families
#' @export
zig_control <- function(min_pos_fit = 4L, k = 24L,
                        lambda = c(nu = 100, l1 = 1000, lcv = 1000),
                        lcv_range = c(0.001, 0.999), min_pos = 10L,
                        min_days = 30L, maxit = 100L, tol = 1e-6) not_implemented()

#' @rdname families
#' @export
family_empirical <- function() not_implemented()

#' Build, check and register a family
#'
#' `new_family()` builds a `normals_family` from its slots and validates it.
#' `validate_family()` checks the slot contract. `register_family()` adds a
#' family to the session registry; `get_family()` looks one up by name
#' (`"zero_inflated_gamma"` and `"empirical"` are registered on load).
#'
#' @param name The `distribution` value written to normals tables.
#' @param params Parameter names, in column order.
#' @param zero_mass `TRUE` if P(X = 0) is part of the family.
#' @param stats_names Column names of `stats()`'s output.
#' @param stats `function(ref, control)`: one day's statistics.
#' @param fit `function(stats, min_valid, control)`: one day's parameters.
#' @param smooth `NULL`, or `list(terms, fit, eval, status)`.
#' @param cdf `function(params, x, ref, min_valid, fallback, control)`.
#' @param quantile `function(params, probs, ref, min_valid, fallback, control)`.
#' @param mean `function(params, ref, min_valid, fallback, control)`.
#' @param p_zero `function(params, ref, min_valid, control)`, or `NULL` when
#'   `zero_mass` is `FALSE`.
#' @param control The family's constants.
#' @param family,x A `normals_family` (or any object for `is_family()`).
#' @return `new_family()`, `register_family()`, `get_family()`: a
#'   `normals_family`. `validate_family()`: `family`, invisibly.
#'   `is_family()`: logical.
#' @export
new_family <- function(name, params, zero_mass, stats_names, stats, fit,
                       smooth, cdf, quantile, mean, p_zero,
                       control = list()) not_implemented()

#' @rdname new_family
#' @export
validate_family <- function(family) not_implemented()

#' @rdname new_family
#' @export
is_family <- function(x) not_implemented()

#' @rdname new_family
#' @export
register_family <- function(family) not_implemented()

#' @rdname new_family
#' @export
get_family <- function(name) not_implemented()
