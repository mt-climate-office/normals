#' One day's distribution for a set of subjects
#'
#' A `normals_dist` bundles a family, its parameters, the reference matrix
#' (for the empirical family and for fallback), `min_valid` and `fallback`.
#' Every argument is required.
#'
#' * `dist_cdf()`: P(X <= x); a zero maps to p_zero / 2 in zero-mass families.
#' * `dist_quantile()`: quantiles at `probs`, a vector shared by every row
#'   (columns named by [quantile_columns()]) or an n x m matrix of row-specific
#'   probabilities.
#' * `dist_mean()`: the mean.
#' * `dist_p_zero()`: the chance of exactly zero (zero-mass families only).
#' * `dist_n_valid()`: non-missing reference years per row.
#' * `mask_sparse()`: sets every parameter `NA` on rows with fewer than
#'   `min_valid` valid years (mco-snotel's smoothed rule).
#'
#' @param family A `normals_family`.
#' @param params n x p parameter matrix, or `NULL` for the empirical family.
#' @param ref n x y reference matrix, or `NULL` when `fallback = "none"` and
#'   the family is not empirical.
#' @param min_valid Integer (length 1 or n), or `NULL` when `ref` is `NULL`.
#' @param fallback `"empirical"` (rows without the family's parameters use the
#'   ECDF of `ref`) or `"none"` (`NA`). May be `NULL` for the empirical family.
#' @param dist A `normals_dist`.
#' @param x Values, length n.
#' @param probs Probabilities: a vector, or an n x m matrix.
#' @return `normals_dist()`: a `normals_dist`. `dist_quantile()`: n x m
#'   matrix. `mask_sparse()`: `params`. The others: a vector of length n.
#' @export
normals_dist <- function(family, params, ref, min_valid, fallback) not_implemented()

#' @rdname normals_dist
#' @export
dist_cdf <- function(dist, x) not_implemented()

#' @rdname normals_dist
#' @export
dist_quantile <- function(dist, probs) not_implemented()

#' @rdname normals_dist
#' @export
dist_mean <- function(dist) not_implemented()

#' @rdname normals_dist
#' @export
dist_p_zero <- function(dist) not_implemented()

#' @rdname normals_dist
#' @export
dist_n_valid <- function(dist) not_implemented()

#' @rdname normals_dist
#' @export
mask_sparse <- function(params, ref, min_valid) not_implemented()
