# Test-only shim for the outlook tests, until the core (R/dist.R, R/family.R,
# R/calendar.R, R/constants.R) lands. Every core function the outlook calls
# that is still a "not implemented" stub is replaced, inside one test, by a
# plain-R version of today's behavior:
#
#   normals_dist()   a list(family, params, ref, min_valid, fallback)
#   dist_cdf()       mco-snodas zig_p(): x = 0 -> p_zero / 2; pgamma above
#                    zero; the zero-inflated type-7 ECDF of `ref` where shape
#                    or scale is NA ("empirical" fallback)
#   dist_quantile()  zig_q(): p_adj = (p - p_zero) / (1 - p_zero); <= 0 -> 0,
#                    >= 1 -> NA; qgamma, or the ECDF fallback
#   dist_p_zero()    the p_zero parameter
#   dist_n_valid()   non-missing reference years per row
#   cycle_day()      normals_cycle_day() (Feb 29 = 151.5 / 59.5)
#   quantile_columns()  sprintf("q%04d", round(p * 1000))
#
# Once the core is implemented nothing is replaced, and the same tests run
# against it. The shim is used by tests/testthat/test-outlook.R only.

shim_family_zig <- function() {
  structure(list(name = "zero_inflated_gamma", params = c("shape", "scale", "p_zero"),
                 zero_mass = TRUE),
            class = "normals_family")
}

shim_ecdf_buf <- function(h, min_valid) {
  v <- h[!is.na(h)]
  if (length(v) < min_valid) return(NULL)
  list(pos = sort(v[v != 0]), n_zero = sum(v == 0), n_valid = length(v))
}

shim_ecdf_p <- function(b, v) {
  n_pos <- length(b$pos); pz <- b$n_zero / b$n_valid
  if (v == 0) return(pz / 2)
  if (n_pos == 0L) return(1)
  lo <- sum(b$pos <= v) - 1L          # 0-based, as std::upper_bound - 1
  p_adj <- if (lo < 0L) 0 else if (lo >= n_pos - 1L) 1 else {
    span <- b$pos[lo + 2L] - b$pos[lo + 1L]
    frac <- if (span > 0) (v - b$pos[lo + 1L]) / span else 0
    (lo + frac) / (n_pos - 1L)
  }
  pz + (1 - pz) * p_adj
}

shim_ecdf_q <- function(b, p) {
  n_pos <- length(b$pos); pz <- b$n_zero / b$n_valid
  if (n_pos == 0L) return(if (p <= 1) 0 else NA_real_)
  p_adj <- (p - pz) / (1 - pz)
  if (p_adj < 0) return(0)
  if (p_adj > 1) return(NA_real_)
  rank <- p_adj * (n_pos - 1L)
  lo <- max(0L, min(as.integer(rank), n_pos - 1L))
  hi <- min(lo + 1L, n_pos - 1L)
  frac <- max(0, min(rank - lo, 1))
  b$pos[lo + 1L] + frac * (b$pos[hi + 1L] - b$pos[lo + 1L])
}

shim_normals_dist <- function(family, params, ref, min_valid, fallback) {
  stopifnot(identical(colnames(params), family$params), nrow(params) == nrow(ref),
            identical(fallback, "empirical"))
  structure(list(family = family, params = params, ref = ref, min_valid = min_valid,
                 fallback = fallback),
            class = "normals_dist")
}

shim_dist_n_valid <- function(dist) as.integer(rowSums(!is.na(dist$ref)))

shim_dist_p_zero <- function(dist) unname(dist$params[, "p_zero"])

shim_dist_cdf <- function(dist, x) {
  p <- dist$params; n <- nrow(p); mv <- rep_len(dist$min_valid, n)
  vapply(seq_len(n), \(i) {
    sh <- p[i, 1]; sc <- p[i, 2]; p0 <- p[i, 3]; v <- x[i]
    if (is.na(p0) || is.na(v)) return(NA_real_)
    if (v == 0) return(p0 / 2)
    if (!is.na(sh) && !is.na(sc)) return(p0 + (1 - p0) * stats::pgamma(v, sh, scale = sc))
    b <- shim_ecdf_buf(dist$ref[i, ], mv[i])
    if (is.null(b)) NA_real_ else shim_ecdf_p(b, v)
  }, 0)
}

shim_dist_quantile <- function(dist, probs) {
  p <- dist$params; n <- nrow(p); mv <- rep_len(dist$min_valid, n)
  u <- if (is.matrix(probs)) probs else matrix(probs, n, length(probs), byrow = TRUE)
  out <- matrix(NA_real_, n, ncol(u))
  for (i in seq_len(n)) {
    sh <- p[i, 1]; sc <- p[i, 2]; p0 <- p[i, 3]
    if (is.na(p0)) next
    emp <- is.na(sh) || is.na(sc)
    b <- if (emp) shim_ecdf_buf(dist$ref[i, ], mv[i])
    for (k in seq_len(ncol(u))) {
      p_adj <- (u[i, k] - p0) / (1 - p0)
      out[i, k] <- if (p_adj <= 0) 0 else if (p_adj >= 1) NA_real_
                   else if (!emp) stats::qgamma(p_adj, sh, scale = sc)
                   else if (is.null(b)) NA_real_ else shim_ecdf_q(b, u[i, k])
    }
  }
  if (!is.matrix(probs)) colnames(out) <- shim_quantile_columns(probs)
  out
}

shim_cycle_day <- function(date, year_type) {
  y  <- as.integer(format(date, "%Y")); m <- as.integer(format(date, "%m"))
  if (year_type == "wy") {
    yr <- y + (m >= 10L)
    start <- as.Date(sprintf("%d-10-01", yr - 1L))
    feb29 <- as.Date(sprintf("%d-02-29", yr), optional = TRUE)
  } else {
    yr <- y
    start <- as.Date(sprintf("%d-01-01", yr))
    feb29 <- as.Date(sprintf("%d-02-29", yr), optional = TRUE)
  }
  pos  <- as.numeric(date - start) + 1
  leap <- !is.na(feb29)
  pos[leap & date >  feb29[leap]] <- pos[leap & date > feb29[leap]] - 1
  pos[leap & date == feb29[leap]] <- pos[leap & date == feb29[leap]] - 0.5
  pos
}

shim_quantile_columns <- function(probs = NORMALS_QUANTILES) {
  sprintf("q%04d", round(unname(probs) * 1000))
}

is_stub <- function(name) {
  identical(body(get(name, envir = asNamespace("normals"))), quote(not_implemented()))
}

# Replaces the core stubs for the calling test.
local_outlook_shim <- function(env = parent.frame()) {
  impl <- list(normals_dist = shim_normals_dist, dist_cdf = shim_dist_cdf,
               dist_quantile = shim_dist_quantile, dist_p_zero = shim_dist_p_zero,
               dist_n_valid = shim_dist_n_valid, cycle_day = shim_cycle_day,
               quantile_columns = shim_quantile_columns)
  impl <- impl[vapply(names(impl), is_stub, TRUE)]
  if (length(impl)) {
    do.call(testthat::local_mocked_bindings, c(impl, list(.package = "normals", .env = env)))
  }
  invisible(length(impl))
}

# The ZIG family: the core's once it exists.
outlook_test_family <- function() {
  if (is_stub("family_zig")) shim_family_zig() else family_zig()
}

# A dist with ZIG parameters (shape, scale, p_zero) and its reference matrix.
zig_dist <- function(params, ref, min_valid) {
  colnames(params) <- c("shape", "scale", "p_zero")
  normals_dist(outlook_test_family(), params, ref, min_valid, "empirical")
}
