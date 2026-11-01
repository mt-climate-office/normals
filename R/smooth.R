#' Smoothing through the year
#'
#' * `smooth_basis()`: the cyclic cubic P-spline basis (mgcv `bs = "cp"`,
#'   knots 0.5 - 365.5, no constraint absorbed) at cycle positions `at`, and
#'   its penalty: `list(X, S, k)`.
#' * `stack_stats()`: a list of `day_stats()` matrices, one per day, as the
#'   stat-major list `smooth_fit()` takes (one n x n_days matrix per
#'   statistic).
#' * `smooth_fit()`: fits the family's smoother to every row. Rows are split
#'   into `cores` contiguous chunks and fitted in forked workers; the result
#'   does not depend on `cores`. Returns a `normals_smooth`: `coef`, `edf`,
#'   `iter`, `status`, `family`, `at`, `k`, `lambda`.
#' * `smooth_eval()`: the parameters at cycle positions `at`, a list (by
#'   parameter) of n x `length(at)` matrices.
#' * `smooth_params()`: the n x p parameter matrix for position `j` of a
#'   `smooth_eval()` result.
#' * `fit_status_label()`: status codes as spec section 8 `fit_status` labels.
#' * `smooth_tables()`: today's coefficient and diagnostics tables (layout
#'   v0), `list(coef, diagnostics)`.
#'
#' @param at Cycle positions ([cycle_day()]).
#' @param k Number of basis functions.
#' @param stats_by_day List of `day_stats()` matrices, one per day.
#' @param stats Stat-major daily statistics (see `stack_stats()`).
#' @param family A `normals_family` with a `smooth` slot.
#' @param cores Forked workers.
#' @param fit A `normals_smooth`.
#' @param evaluated A `smooth_eval()` result.
#' @param j Position index.
#' @param status Integer status codes.
#' @param subjects Data frame of subject key columns, one row per fit row.
#' @name smoothing
NULL

#' @rdname smoothing
#' @export
smooth_basis <- function(at, k) {
  if (!is.numeric(at) || !length(at) || anyNA(at) || any(at < 0.5 | at > 365.5))
    stop_arg("at", "must be cycle positions in [0.5, 365.5]")
  check_whole(k, "k", positive = TRUE, single = TRUE)
  # mco-snodas zig_spline_basis(), unchanged. Knots depend only on the cycle
  # range, so the basis does not depend on `at`; mgcv's warning when `at`
  # covers little of the cycle is harmless.
  sm <- suppressWarnings(mgcv::smoothCon(
    mgcv::s(day, bs = "cp", k = k),
    data  = data.frame(day = as.numeric(at)),
    knots = list(day = c(0.5, 365.5)),   # day 365 wraps to day 1
    absorb.cons = FALSE
  ))[[1]]
  list(X = sm$X, S = sm$S[[1]], k = k)
}

#' @rdname smoothing
#' @export
stack_stats <- function(stats_by_day) {
  if (!is.list(stats_by_day) || !length(stats_by_day) ||
      !all(vapply(stats_by_day, is.matrix, logical(1))))
    stop_arg("stats_by_day", "must be a list of day_stats() matrices")
  vars <- colnames(stats_by_day[[1]])
  n    <- nrow(stats_by_day[[1]])
  if (is.null(vars) || !all(vapply(stats_by_day, \(m) identical(colnames(m), vars) && nrow(m) == n,
                                   logical(1))))
    stop_arg("stats_by_day", "must hold matrices with the same rows and columns")
  # Filled day by day, as the pipelines do today.
  out <- lapply(stats::setNames(vars, vars), \(v) matrix(NA_real_, n, length(stats_by_day)))
  for (j in seq_along(stats_by_day))
    for (v in vars) out[[v]][, j] <- stats_by_day[[j]][, v]
  out
}

#' @rdname smoothing
#' @export
smooth_fit <- function(stats, at, family, cores = 1L) {
  check_family(family)
  if (is.null(family$smooth))
    stop("normals: family \"", family$name, "\" has no smoother", call. = FALSE)
  if (!is.list(stats) || !all(family$stats_names %in% names(stats)))
    stop_arg("stats", "must be a list with ", paste(family$stats_names, collapse = ", "),
             " (stack_stats())")
  n <- nrow(stats[[family$stats_names[1]]])
  for (v in family$stats_names) {
    m <- stats[[v]]
    if (!is.matrix(m) || nrow(m) != n || ncol(m) != length(at))
      stop_arg("stats", "must hold n x length(at) matrices (", v, " does not)")
  }
  check_whole(cores, "cores", positive = TRUE, single = TRUE)
  basis <- smooth_basis(at, family$control$k)
  st  <- lapply(stats[family$stats_names], as_double_matrix)
  res <- family$smooth$fit(st, basis, as.integer(cores), family$control)
  structure(
    list(coef = res$coef, edf = res$edf, iter = res$iter, status = as.integer(res$status),
         family = family, at = at, k = family$control$k, lambda = family$control$lambda),
    class = "normals_smooth")
}

#' @rdname smoothing
#' @export
smooth_eval <- function(fit, at) {
  check_smooth(fit)
  X <- smooth_basis(at, fit$k)$X
  # One position at a time, as both pipelines evaluate today (one basis row
  # per date): coef %*% t(X) for one row and for many need not agree to the
  # last bit under every BLAS.
  fam <- fit$family
  per <- lapply(seq_len(nrow(X)), \(j)
    fam$smooth$eval(fit$coef, X[j, , drop = FALSE], fit$status, fam$control))
  stats::setNames(lapply(fam$params, \(p) {
    out <- matrix(NA_real_, nrow(fit$coef), nrow(X))
    for (j in seq_along(per)) out[, j] <- per[[j]][[p]][, 1]
    out
  }), fam$params)
}

#' @rdname smoothing
#' @export
smooth_params <- function(evaluated, j) {
  if (!is.list(evaluated) || is.null(names(evaluated)) || !all(vapply(evaluated, is.matrix, logical(1))))
    stop_arg("evaluated", "must be a smooth_eval() result")
  check_whole(j, "j", positive = TRUE, single = TRUE)
  if (j > ncol(evaluated[[1]])) stop_arg("j", "is beyond the evaluated positions")
  out <- do.call(cbind, lapply(evaluated, \(m) m[, j]))
  if (!is.matrix(out)) out <- matrix(out, 1L)
  colnames(out) <- names(evaluated)
  out
}

#' @rdname smoothing
#' @export
fit_status_label <- function(status, family) {
  check_family(family)
  if (is.null(family$smooth)) stop("normals: family \"", family$name, "\" has no smoother", call. = FALSE)
  codes <- family$smooth$status
  out <- names(codes)[match(status, codes)]
  if (any(is.na(out) & !is.na(status)))
    stop_arg("status", "has codes unknown to family \"", family$name, "\"")
  out
}

#' @rdname smoothing
#' @export
smooth_tables <- function(fit, subjects) {
  check_smooth(fit)
  check_subjects(subjects, nrow(fit$coef))
  # mco-snotel calc_station_normals() (layout v0).
  k <- fit$k
  terms <- names(fit$lambda)
  coef <- fit$coef
  colnames(coef) <- sprintf("%s_%02d", rep(terms, each = k), rep(seq_len(k), length(terms)))
  keys <- tibble::as_tibble(as.list(subjects))
  diag_cols <- c(
    list(status = as.integer(fit$status), iter = as.integer(fit$iter)),
    stats::setNames(as.list(unname(fit$lambda)), paste0("lambda_", terms)),
    stats::setNames(lapply(seq_along(terms), \(i) fit$edf[, i]), paste0("edf_", terms)))
  list(
    coef        = tibble::as_tibble(c(as.list(keys), as.list(as.data.frame(coef)))),
    diagnostics = tibble::as_tibble(c(as.list(keys), diag_cols))
  )
}

check_smooth <- function(fit) {
  if (!inherits(fit, "normals_smooth")) stop_arg("fit", "must be a normals_smooth (smooth_fit())")
  invisible(fit)
}

check_subjects <- function(subjects, n) {
  if (!is.data.frame(subjects) || !ncol(subjects))
    stop_arg("subjects", "must be a data frame of key columns")
  if (nrow(subjects) != n) stop_arg("subjects", "must have one row per subject (", n, ")")
  invisible(subjects)
}

utils::globalVariables("day")
