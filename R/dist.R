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
normals_dist <- function(family, params, ref, min_valid, fallback) {
  check_family(family)
  empirical <- !length(family$params)
  if (!is.null(ref)) ref <- as_double_matrix(check_matrix(ref, "ref"))
  if (empirical) {
    if (!is.null(params)) stop_arg("params", "must be NULL for a family without parameters")
    if (is.null(ref)) stop_arg("ref", "is required for a family without parameters")
    n <- nrow(ref)
    fallback <- NULL
  } else {
    check_matrix(params, "params")
    if (ncol(params) != length(family$params) ||
        (!is.null(colnames(params)) && !identical(colnames(params), family$params)))
      stop_arg("params", "must have columns ", paste(family$params, collapse = ", "))
    params <- as_double_matrix(params)
    colnames(params) <- family$params
    n <- nrow(params)
    check_choice(fallback, c("empirical", "none"), "fallback")
    if (fallback == "empirical" && is.null(ref))
      stop_arg("ref", "is required with `fallback = \"empirical\"`")
    if (!is.null(ref) && nrow(ref) != n)
      stop_arg("ref", "must have one row per row of `params` (", n, ")")
  }
  if (is.null(ref)) {
    if (!is.null(min_valid)) stop_arg("min_valid", "must be NULL when `ref` is NULL")
  } else {
    if (is.null(min_valid)) stop_arg("min_valid", "is required with `ref`")
    min_valid <- as_min_valid(min_valid, n)
  }
  structure(list(family = family, params = params, ref = ref, min_valid = min_valid,
                 fallback = fallback),
            class = "normals_dist")
}

#' @rdname normals_dist
#' @export
dist_cdf <- function(dist, x) {
  check_dist(dist)
  if (!is.numeric(x) || length(x) != dist_rows(dist)) stop_arg("x", "must be numeric, one value per row (", dist_rows(dist), ")")
  f <- dist$family
  f$cdf(dist$params, x, dist$ref, dist$min_valid, dist$fallback, f$control)
}

#' @rdname normals_dist
#' @export
dist_quantile <- function(dist, probs) {
  check_dist(dist)
  if (!is.numeric(probs) || !length(probs))
    stop_arg("probs", "must be probabilities: a vector, or a matrix with one row per row")
  if (is.matrix(probs)) {
    if (nrow(probs) != dist_rows(dist)) stop_arg("probs", "must have one row per row (", dist_rows(dist), ")")
    names_out <- NULL
  } else {
    names_out <- quantile_columns(probs)
    probs <- unname(probs)
  }
  f <- dist$family
  out <- f$quantile(dist$params, probs, dist$ref, dist$min_valid, dist$fallback, f$control)
  colnames(out) <- names_out
  out
}

#' @rdname normals_dist
#' @export
dist_mean <- function(dist) {
  check_dist(dist)
  f <- dist$family
  f$mean(dist$params, dist$ref, dist$min_valid, dist$fallback, f$control)
}

#' @rdname normals_dist
#' @export
dist_p_zero <- function(dist) {
  check_dist(dist)
  f <- dist$family
  if (!f$zero_mass) stop("normals: family \"", f$name, "\" has no zero mass", call. = FALSE)
  if (!length(f$params) && is.null(dist$ref)) stop_arg("dist", "has no reference matrix")
  f$p_zero(dist$params, dist$ref, dist$min_valid, f$control)
}

#' @rdname normals_dist
#' @export
dist_n_valid <- function(dist) {
  check_dist(dist)
  if (is.null(dist$ref)) stop_arg("dist", "has no reference matrix")
  as.integer(rowSums(!is.na(dist$ref)))
}

#' @rdname normals_dist
#' @export
mask_sparse <- function(params, ref, min_valid) {
  check_matrix(params, "params")
  check_matrix(ref, "ref")
  if (nrow(ref) != nrow(params)) stop_arg("ref", "must have one row per row of `params`")
  mv <- as_min_valid(min_valid, nrow(params))
  # mco-snotel calc_station_normals() (the smoothed flavor).
  params[rowSums(!is.na(ref)) < mv, ] <- NA_real_
  params
}

check_dist <- function(dist) {
  if (!inherits(dist, "normals_dist")) stop_arg("dist", "must be a normals_dist")
  invisible(dist)
}

dist_rows <- function(dist) if (is.null(dist$params)) nrow(dist$ref) else nrow(dist$params)
