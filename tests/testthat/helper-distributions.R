# scripts/test_distributions.R (mco-snodas 9b45d74, vendored in mco-snotel),
# Parts 1-4, as testthat. The checks and tolerances are unchanged; today's
# function names come from helper-legacy.R. The R reference implementations
# below are copied verbatim.

# Evaluate expr after set.seed(seed), restoring the session's RNG state.
withr_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()))
  set.seed(seed)
  expr
}

TOL <- 1e-10   # machine-precision threshold for double comparisons

r_ecdf_p <- function(hist_mat, obs) {
  vapply(seq_len(nrow(hist_mat)), function(i) {
    x      <- hist_mat[i, ]
    v      <- obs[i]
    if (is.na(v) || sum(!is.na(x)) < min(10L, length(x))) return(NA_real_)
    x      <- x[!is.na(x)]
    n      <- length(x)
    pos    <- sort(x[x > 0])
    n_pos  <- length(pos)
    p_zero <- sum(x == 0) / n
    if (v == 0)     return(p_zero / 2)
    if (n_pos == 0) return(1.0)
    lo    <- findInterval(v, pos)    # last 1-based index where pos[lo] <= v
    p_adj <- if (lo == 0L) {
      0.0
    } else if (lo >= n_pos) {
      1.0
    } else {
      span <- pos[lo + 1L] - pos[lo]
      frac <- if (span > 0) (v - pos[lo]) / span else 0
      ((lo - 1L) + frac) / (n_pos - 1)
    }
    p_zero + (1 - p_zero) * p_adj
  }, numeric(1))
}

# Inverse CDF: zero-inflated type-7, matching C++ ecdf_q().
r_ecdf_q <- function(hist_mat, probs) {
  t(apply(hist_mat, 1, function(x) {
    n      <- length(x)
    pos    <- sort(x[x > 0])
    n_pos  <- length(pos)
    p_zero <- sum(x == 0) / n
    if (n_pos == 0L) return(ifelse(probs <= 1, 0, NA_real_))
    vapply(probs, function(p) {
      p_adj <- (p - p_zero) / (1 - p_zero)
      if (p_adj < 0) return(0)
      if (p_adj > 1) return(NA_real_)
      quantile(pos, p_adj, type = 7, names = FALSE)
    }, numeric(1))
  }))
}

ECDF_PROBS  <- c(0, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3,
                 0.5, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.98, 0.99, 1.0)
SAMPLE_SIZES <- c(5L, 10L, 22L, 50L, 100L)
ZERO_FRACS   <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9)
ECDF_POOL_DBL <- withr_seed(42, sort(round(rgamma(max(SAMPLE_SIZES), shape = 4, rate = 0.2), 4)))
ECDF_POOL_INT <- as.integer(round(ECDF_POOL_DBL * 25.4))   # inches -> mm
pool_ecdf <- function(dt) if (dt == "double") ECDF_POOL_DBL else ECDF_POOL_INT
make_raw_ecdf <- function(n_zero, pool) c(rep(pool[1L] * 0L, n_zero), pool)

r_zig_p_param <- function(shape, scale, p_zero, v)
  p_zero + (1 - p_zero) * pgamma(v, shape = shape, scale = scale)

r_zig_q <- function(shape, scale, p_zero, probs)
  vapply(probs, function(p) {
    p_adj <- (p - p_zero) / (1 - p_zero)
    if (p_adj <= 0) return(0.0)
    if (p_adj >= 1) return(NA_real_)
    qgamma(p_adj, shape = shape, scale = scale)
  }, numeric(1))

r_zig_p_fallback <- function(raw, obs)
  vapply(obs, function(v) {
    mat <- matrix(raw, nrow = 1L)
    r_ecdf_p(mat, v)
  }, numeric(1))

r_zig_q_fallback <- function(raw, probs) {
  # Matches zig_q convention: p_adj <= 0 (inclusive) maps to 0.
  # r_ecdf_q uses strict p_adj < 0, so we handle the boundary here.
  n <- length(raw); pos <- sort(raw[raw > 0]); n_pos <- length(pos)
  p_zero <- sum(raw == 0) / n
  vapply(probs, function(p) {
    if (n_pos == 0) return(if (p <= 1) 0.0 else NA_real_)
    p_adj <- (p - p_zero) / (1 - p_zero)
    if (p_adj <= 0) return(0.0)
    if (p_adj >= 1) return(NA_real_)
    quantile(pos, p_adj, type = 7, names = FALSE)
  }, numeric(1))
}

ZIG_PROBS  <- ECDF_PROBS
ZIG_SIZES  <- c(10L, 22L, 50L, 100L)
ZIG_POOL <- withr_seed(42, sort(round(rgamma(max(ZIG_SIZES), shape = 3, scale = 10), 2)))  # mm

r_lmom <- function(x) {
  v <- sort(x[!is.na(x) & x > 0]); n <- length(v)
  out <- c(n = sum(!is.na(x)), z = sum(x == 0, na.rm = TRUE), l1 = NA, lcv = NA)
  if (n < 4) return(out)
  b0 <- mean(v); b1 <- sum((seq_len(n) - 1) / (n - 1) * v) / n
  l2 <- 2 * b1 - b0
  if (l2 > 0 && b0 > 0) out[c("l1", "lcv")] <- c(b0, l2 / b0)
  out
}
