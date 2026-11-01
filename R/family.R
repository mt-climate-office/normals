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
family_zig <- function(control = zig_control()) {
  if (!is.list(control) || !identical(names(control), names(zig_control())))
    stop_arg("control", "must be a zig_control() list")
  new_family(
    name        = "zero_inflated_gamma",
    params      = c("shape", "scale", "p_zero"),
    zero_mass   = TRUE,
    stats_names = c("n", "z", "l1", "lcv"),
    stats       = zig_stats,
    fit         = zig_fit_day,
    smooth      = list(terms  = ZIG_SMOOTH_TERMS,
                       fit    = zig_smooth_fit,
                       eval   = zig_smooth_eval,
                       status = ZIG_STATUS),
    cdf         = zig_dist_cdf,
    quantile    = zig_dist_quantile,
    mean        = zig_dist_mean,
    p_zero      = zig_dist_p_zero,
    control     = control
  )
}

#' @rdname families
#' @export
zig_control <- function(min_pos_fit = 4L, k = 24L,
                        lambda = c(nu = 100, l1 = 1000, lcv = 1000),
                        lcv_range = c(0.001, 0.999), min_pos = 10L,
                        min_days = 30L, maxit = 100L, tol = 1e-6) {
  check_whole(min_pos_fit, "min_pos_fit", positive = TRUE, single = TRUE)
  check_whole(k, "k", positive = TRUE, single = TRUE)
  if (!is.numeric(lambda) || length(lambda) != 3L || anyNA(lambda) || any(lambda < 0) ||
      !setequal(names(lambda), c("nu", "l1", "lcv")))
    stop_arg("lambda", "must be c(nu = , l1 = , lcv = ), non-negative")
  if (!is.numeric(lcv_range) || length(lcv_range) != 2L || anyNA(lcv_range) ||
      lcv_range[1] <= 0 || lcv_range[1] >= lcv_range[2])
    stop_arg("lcv_range", "must be c(lower, upper) with 0 < lower < upper")
  check_whole(min_pos, "min_pos", single = TRUE)
  check_whole(min_days, "min_days", single = TRUE)
  check_whole(maxit, "maxit", positive = TRUE, single = TRUE)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0)
    stop_arg("tol", "must be a positive number")
  list(min_pos_fit = as.integer(min_pos_fit), k = as.integer(k),
       lambda = lambda[c("nu", "l1", "lcv")], lcv_range = lcv_range,
       min_pos = as.integer(min_pos), min_days = as.integer(min_days),
       maxit = as.integer(maxit), tol = tol)
}

#' @rdname families
#' @export
family_empirical <- function() {
  new_family(
    name        = "empirical",
    params      = character(),
    zero_mass   = TRUE,
    stats_names = character(),
    stats       = NULL,
    fit         = NULL,
    smooth      = NULL,
    cdf         = emp_cdf,
    quantile    = emp_quantile,
    mean        = emp_mean,
    p_zero      = emp_p_zero,
    control     = list()
  )
}

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
                       control = list()) {
  family <- structure(
    list(name = name, params = params, zero_mass = zero_mass,
         stats_names = stats_names, stats = stats, fit = fit, smooth = smooth,
         cdf = cdf, quantile = quantile, mean = mean, p_zero = p_zero,
         control = control),
    class = "normals_family")
  validate_family(family)
}

#' @rdname new_family
#' @export
validate_family <- function(family) {
  bad <- function(...) stop("normals: invalid family: ", ..., call. = FALSE)
  if (!inherits(family, "normals_family") || !is.list(family)) bad("not a normals_family")
  slots <- c("name", "params", "zero_mass", "stats_names", "stats", "fit", "smooth",
             "cdf", "quantile", "mean", "p_zero", "control")
  missing <- setdiff(slots, names(family))
  if (length(missing)) bad("missing slot(s) ", paste(missing, collapse = ", "))
  f <- family
  if (!is.character(f$name) || length(f$name) != 1L || is.na(f$name) || !nzchar(f$name))
    bad("`name` must be a string")
  if (!is.character(f$params) || anyNA(f$params) || anyDuplicated(f$params))
    bad("`params` must be distinct names")
  reserved <- c("flavor", "date", "day", "distribution", "n_valid", "ref_first_wy", "ref_last_wy",
                "mean", "fit_status", "value", "ptile", "standardized")
  if (any(f$params %in% reserved | grepl("^q[0-9]{4}$", f$params)))
    bad("`params` must not reuse a table column name (", paste(intersect(f$params, reserved), collapse = ", "), ")")
  if (!is.logical(f$zero_mass) || length(f$zero_mass) != 1L || is.na(f$zero_mass))
    bad("`zero_mass` must be TRUE or FALSE")
  if (!is.character(f$stats_names) || anyNA(f$stats_names) || anyDuplicated(f$stats_names))
    bad("`stats_names` must be distinct names")
  if (!is.list(f$control)) bad("`control` must be a list")
  has_formals <- function(fn, args) is.function(fn) && identical(names(formals(fn)), args)
  check_fn <- function(slot, args, optional = FALSE) {
    fn <- f[[slot]]
    if (optional && is.null(fn)) return(invisible())
    if (!has_formals(fn, args))
      bad("`", slot, "` must be function(", paste(args, collapse = ", "), ")")
  }
  check_fn("stats", c("ref", "control"), optional = TRUE)
  check_fn("fit", c("stats", "min_valid", "control"), optional = TRUE)
  if (is.null(f$fit) != (length(f$params) == 0L))
    bad("`params` must be non-empty exactly when `fit` is given")
  if (!is.null(f$fit) && (is.null(f$stats) || !length(f$stats_names)))
    bad("a family with `fit` needs `stats` and `stats_names`")
  if (!is.null(f$smooth)) {
    s <- f$smooth
    if (!is.list(s) || !all(c("terms", "fit", "eval", "status") %in% names(s)))
      bad("`smooth` must be NULL or list(terms, fit, eval, status)")
    if (is.null(f$fit)) bad("a family with `smooth` needs `fit`")
    if (!is.data.frame(s$terms) || !all(c("term", "scale", "method", "weight") %in% names(s$terms)))
      bad("`smooth$terms` must be a data frame with term, scale, method, weight")
    if (!has_formals(s$fit, c("stats", "basis", "cores", "control")))
      bad("`smooth$fit` must be function(stats, basis, cores, control)")
    if (!has_formals(s$eval, c("coef", "X", "status", "control")))
      bad("`smooth$eval` must be function(coef, X, status, control)")
    if (!is.integer(s$status) || is.null(names(s$status)) || anyNA(s$status))
      bad("`smooth$status` must be a named integer vector")
  }
  check_fn("cdf", c("params", "x", "ref", "min_valid", "fallback", "control"))
  check_fn("quantile", c("params", "probs", "ref", "min_valid", "fallback", "control"))
  check_fn("mean", c("params", "ref", "min_valid", "fallback", "control"))
  if (f$zero_mass) check_fn("p_zero", c("params", "ref", "min_valid", "control"))
  else if (!is.null(f$p_zero)) bad("`p_zero` must be NULL when `zero_mass` is FALSE")
  invisible(family)
}

#' @rdname new_family
#' @export
is_family <- function(x) inherits(x, "normals_family")

#' @rdname new_family
#' @export
register_family <- function(family) {
  validate_family(family)
  old <- .families[[family$name]]
  if (!is.null(old) && !identical(old, family, ignore.environment = TRUE))
    stop("normals: a different family named \"", family$name, "\" is already registered",
         call. = FALSE)
  assign(family$name, family, envir = .families)
  invisible(family)
}

#' @rdname new_family
#' @export
get_family <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) stop_arg("name", "must be a string")
  out <- .families[[name]]
  if (is.null(out)) stop("normals: no family named \"", name, "\" is registered", call. = FALSE)
  out
}

# The session registry; zero_inflated_gamma and empirical are added on load.
.families <- new.env(parent = emptyenv())


# ── Zero-inflated gamma: the slots ──────────────────────────────────────────
# Defined once at top level, so every family_zig() with the same control is
# identical() (the registry relies on it).

ZIG_STATUS <- c(ok = 0L, no_snow = 1L, p_zero_only = 2L, no_data = 3L,
                p_zero_not_converged = 4L, constant = 5L)

ZIG_SMOOTH_TERMS <- data.frame(
  term   = c("p_zero", "l1", "lcv"),
  scale  = c("logit", "log", "log"),
  method = c("binomial", "wls", "wls"),
  weight = c("years (z of n)", "positive years", "positive years"),
  stringsAsFactors = FALSE
)

# One day's L-moment statistics (mco-snodas zig_lmom() / zig_spline_suffstats()).
zig_stats <- function(ref, control) {
  cpp_zig_lmom(as_double_matrix(as.matrix(ref)), control$min_pos_fit)
}

# Parameters from the statistics (zig_params_from_lmom()).
zig_fit_day <- function(stats, min_valid, control) {
  out <- cpp_zig_params_from_lmom(as_double_matrix(stats), min_valid)
  colnames(out) <- c("shape", "scale", "p_zero")
  out
}

# The smoother: mco-snodas zig_spline_chunk_stats() + zig_spline_fit_parallel()
# + zig_spline_fit_stats(), unchanged.
zig_smooth_fit <- function(stats, basis, cores, control) {
  fit_one <- function(chunk) {
    cpp_zig_spline_fit(
      n_obs = chunk$n, n_zero = chunk$z, l1 = chunk$l1, lcv = chunk$lcv,
      basis = basis$X, penalty = basis$S,
      lambda = unname(control$lambda[c("nu", "l1", "lcv")]),
      min_pos = control$min_pos, min_days = control$min_days,
      maxit = control$maxit, tol = control$tol)
  }
  fits <- fit_chunks(chunk_rows(stats[c("n", "z", "l1", "lcv")], cores), fit_one, cores)
  fits <- unname(fits)
  list(
    coef   = do.call(rbind, lapply(fits, `[[`, "coef")),
    edf    = do.call(rbind, lapply(fits, `[[`, "edf")),
    iter   = unlist(lapply(fits, `[[`, "iter")),
    status = unlist(lapply(fits, `[[`, "status"))
  )
}

# Parameters on basis rows X: mco-snodas zig_spline_params(), unchanged.
zig_smooth_eval <- function(coef, X, status, control) {
  k <- ncol(X)
  eta <- \(cols) coef[, cols, drop = FALSE] %*% t(X)
  p_zero <- stats::plogis(eta(seq_len(k)))
  l1     <- exp(eta(k + seq_len(k)))
  lcv    <- pmin(pmax(exp(eta(2L * k + seq_len(k))), control$lcv_range[1]), control$lcv_range[2])
  shape  <- matrix(cpp_zig_shape_from_lcv(as.vector(lcv)), nrow(lcv))
  if (!is.null(status))
    p_zero[status == 1L, ] <- 1
  list(shape = shape, scale = l1 / shape, p_zero = p_zero)
}

no_ref <- function() matrix(NA_real_, 0L, 0L)

zig_dist_cdf <- function(params, x, ref, min_valid, fallback, control) {
  fb <- identical(fallback, "empirical")
  cpp_zig_p(params, as.double(x), if (fb) ref else no_ref(), if (fb) min_valid else 1L, fb)
}

zig_dist_quantile <- function(params, probs, ref, min_valid, fallback, control) {
  if (identical(fallback, "empirical")) {
    pm <- if (is.matrix(probs)) probs else matrix(probs, 1L)
    return(cpp_zig_q(params, as_double_matrix(pm), ref, min_valid, TRUE))
  }
  zig_smooth_q(params, probs)
}

zig_dist_mean <- function(params, ref, min_valid, fallback, control) {
  if (identical(fallback, "empirical")) return(cpp_zig_mean(params, ref, min_valid, TRUE))
  zig_smooth_mean(params)
}

zig_dist_p_zero <- function(params, ref, min_valid, control) unname(params[, 3])

# Quantiles and mean with no ECDF fallback: mco-snodas zig_smooth_q() and
# zig_smooth_mean() (R/calc_snodas_zig_spline_normals.R), unchanged but for
# row-specific probabilities (an n x m matrix, read column by column).
# Where shape and scale are NA, quantiles above p_zero and the mean are NA;
# p_zero >= 1 (no snow on record) gives 0.
zig_smooth_q <- function(params, probs) {
  shape <- params[, 1]; scale <- params[, 2]; p0 <- params[, 3]
  cols <- if (is.matrix(probs)) seq_len(ncol(probs)) else seq_along(probs)
  out <- vapply(cols, \(j) {
    p <- if (is.matrix(probs)) probs[, j] else probs[j]
    p_adj <- (p - p0) / (1 - p0)
    q <- suppressWarnings(stats::qgamma(pmin(pmax(p_adj, 0), 1), shape = shape, scale = scale))
    q[p_adj <= 0] <- 0
    q[p_adj >= 1] <- NA_real_
    q[p0 >= 1]    <- 0          # no snow on record
    q[is.na(p0)]  <- NA_real_
    q
  }, numeric(nrow(params)))
  matrix(out, nrow(params))
}

zig_smooth_mean <- function(params) {
  ifelse(params[, 3] >= 1, 0, (1 - params[, 3]) * params[, 1] * params[, 2])
}


# ── Empirical: the slots ────────────────────────────────────────────────────

emp_cdf <- function(params, x, ref, min_valid, fallback, control) {
  cpp_ecdf_p(ref, as.double(x), min_valid)
}

emp_quantile <- function(params, probs, ref, min_valid, fallback, control) {
  pm <- if (is.matrix(probs)) probs else matrix(probs, 1L)
  cpp_ecdf_q(ref, as_double_matrix(pm), min_valid)
}

emp_mean <- function(params, ref, min_valid, fallback, control) {
  cpp_ecdf_mean(ref, min_valid)
}

emp_p_zero <- function(params, ref, min_valid, control) {
  n <- rowSums(!is.na(ref))
  ifelse(n >= min_valid, rowSums(ref == 0, na.rm = TRUE) / pmax(1L, n), NA_real_)
}


# ── Row chunks for forked workers ───────────────────────────────────────────

# Split stat-major statistics into contiguous row chunks, one per worker
# (mco-snodas zig_spline_chunk_stats()).
chunk_rows <- function(stats, n_chunks) {
  n <- nrow(stats[[1]])
  k <- max(1L, min(n_chunks, n))
  # cut() can't make one interval from one value, so one chunk is built directly.
  rows <- if (k == 1L) list(seq_len(n)) else split(seq_len(n), cut(seq_len(n), k, labels = FALSE))
  lapply(rows, \(r) lapply(stats, \(m) m[r, , drop = FALSE]))
}

# Fit chunks, forked when cores > 1 (mco-snodas zig_spline_fit_parallel()).
fit_chunks <- function(chunks, fit_one, cores) {
  fits <- if (cores <= 1L || length(chunks) == 1L) lapply(chunks, fit_one)
          else parallel::mclapply(chunks, fit_one, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(fits, \(f) is.null(f) || inherits(f, "try-error"), logical(1))
  if (any(bad))
    stop("normals: the smoothed fit failed in ", sum(bad), " worker(s)",
         if (inherits(fits[[which(bad)[1]]], "try-error")) paste0(": ", fits[[which(bad)[1]]])
         else " (killed, likely out of memory)", call. = FALSE)
  fits
}

.onLoad <- function(libname, pkgname) {
  register_family(family_zig())
  register_family(family_empirical())
}
