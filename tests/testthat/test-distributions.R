# scripts/test_distributions.R, Parts 1-4 (see helper-distributions.R).

# ── Part 1: ECDF ──────────────────────────────────────────────────────────────

run_ecdf_scenario <- function(n, fz, dt) {
  pool   <- pool_ecdf(dt)[seq_len(n)]
  n_zero <- min(round(fz * n), n)
  n_pos  <- n - n_zero
  raw    <- if (n_pos == 0L) rep(pool[1L] * 0L, n) else make_raw_ecdf(n_zero, pool[seq_len(n_pos)])
  mat    <- matrix(raw, nrow = 1L)
  cpp_q <- ecdf_stats(mat)$quantiles[1L, ]
  r_q   <- r_ecdf_q(mat, ECDF_PROBS)[1L, ]
  max_q <- max(abs(cpp_q - r_q), na.rm = TRUE)
  positives <- sort(unique(raw[raw > 0]))
  above_max <- if (length(positives)) max(positives) * 2 else if (dt == "integer") 1L else 1.0
  obs_vals  <- if (dt == "integer") as.integer(c(0, positives, above_max)) else
                                    as.double( c(0, positives, above_max))
  mat_rep   <- mat[rep(1L, length(obs_vals)), , drop = FALSE]
  max_p     <- max(abs(ecdf_p(mat_rep, obs_vals) - r_ecdf_p(mat_rep, obs_vals)), na.rm = TRUE)
  max_rt <- NA_real_
  if (length(positives)) {
    obs_rt <- if (dt == "integer") as.integer(positives) else as.double(positives)
    p_vals <- ecdf_p(mat[rep(1L, length(positives)), , drop = FALSE], obs_rt)
    q_vals <- vapply(p_vals, function(p) ecdf_q(mat, p)[1L, 1L], numeric(1))
    max_rt <- max(abs(q_vals - as.double(positives)), na.rm = TRUE)
  }
  list(max_q = max_q, max_p = max_p, max_rt = max_rt)
}

test_that("1-1 to 1-3: ecdf_q and ecdf_p match type 7; ecdf_p -> ecdf_q round-trips", {
  for (n in SAMPLE_SIZES) for (fz in ZERO_FRACS) for (dt in c("double", "integer")) {
    res <- run_ecdf_scenario(n, fz, dt)
    lab <- sprintf("n=%d fz=%.2f %s", n, fz, dt)
    expect_lt(res$max_q, TOL, label = paste(lab, "q"))
    expect_lt(res$max_p, TOL, label = paste(lab, "p"))
    if (!is.na(res$max_rt)) expect_lt(res$max_rt, TOL, label = paste(lab, "round trip"))
  }
})

test_that("1-4: all-zero records", {
  for (n in SAMPLE_SIZES) for (dt in c("double", "integer")) {
    mat <- matrix(if (dt == "integer") integer(n) else double(n), nrow = 1L)
    expect_true(all(ecdf_stats(mat)$quantiles[1L, ] == 0, na.rm = TRUE))
    expect_true(is.na(ecdf_q(mat, 1.5)[1L, 1L]))
    expect_equal(ecdf_p(mat, if (dt == "integer") 0L else 0.0), 0.5)
    expect_equal(ecdf_p(mat, if (dt == "integer") 100L else 10.0), 1.0)
    expect_equal(ecdf_mean(mat), 0.0)
  }
})

test_that("1-5: integer and double input give identical output", {
  for (n in SAMPLE_SIZES) for (fz in ZERO_FRACS) {
    n_zero <- min(round(fz * n), n); n_pos <- n - n_zero
    pool   <- ECDF_POOL_INT[seq_len(n)]
    raw_i  <- if (n_pos == 0L) integer(n) else as.integer(make_raw_ecdf(n_zero, pool[seq_len(n_pos)]))
    mat_i  <- matrix(raw_i, nrow = 1L); mat_d <- matrix(as.double(raw_i), nrow = 1L)
    expect_identical(max(abs(ecdf_stats(mat_d)$quantiles - ecdf_stats(mat_i)$quantiles), na.rm = TRUE), 0)
    obs_i <- c(0L, as.integer(sort(unique(raw_i[raw_i > 0]))))
    nr <- length(obs_i)
    expect_identical(max(abs(ecdf_p(mat_d[rep(1L, nr), , drop = FALSE], as.double(obs_i)) -
                             ecdf_p(mat_i[rep(1L, nr), , drop = FALSE], obs_i)), na.rm = TRUE), 0)
  }
})

test_that("1-6: ecdf_p(0) = p_zero / 2 as p_zero -> 1", {
  for (nz in c(seq(0L, 90L, 10L), 95L, 99L, 100L)) {
    mat <- matrix(c(rep(0.0, nz), rep(10.0, 100L - nz)), nrow = 1L)
    expect_equal(ecdf_p(mat, 0.0), nz / 100 / 2)
  }
})

# ── Part 2: ZIG ───────────────────────────────────────────────────────────────

run_zig_scenario <- function(n, fz) {
  n_zero <- min(round(fz * n), n); n_pos <- n - n_zero
  raw    <- c(rep(0.0, n_zero), ZIG_POOL[seq_len(max(n_pos, 0L))])
  mat    <- matrix(raw, nrow = 1L)
  params <- zig_fit(mat)
  shape  <- params[1L, 1L]; scale <- params[1L, 2L]; pz <- params[1L, 3L]
  param  <- !is.na(shape) && !is.na(scale)
  pos    <- sort(unique(raw[raw > 0]))
  obs    <- c(0.0, pos, if (length(pos)) max(pos) * 2 else 100.0)
  nr     <- length(obs)
  cpp_p  <- zig_p(params[rep(1L, nr), , drop = FALSE], obs, mat[rep(1L, nr), , drop = FALSE])
  r_p    <- if (param) vapply(obs, function(v) if (v == 0) pz / 2 else r_zig_p_param(shape, scale, pz, v),
                              numeric(1)) else r_zig_p_fallback(raw, obs)
  cpp_q  <- zig_q(params, ZIG_PROBS, mat)[1L, ]
  r_q    <- if (param) r_zig_q(shape, scale, pz, ZIG_PROBS) else r_zig_q_fallback(raw, ZIG_PROBS)
  max_rt <- NA_real_
  if (param && length(pos)) {
    p_rt <- zig_p(params[rep(1L, length(pos)), , drop = FALSE], pos, mat[rep(1L, length(pos)), , drop = FALSE])
    q_rt <- vapply(p_rt, function(p) zig_q(params, p, mat)[1L, 1L], numeric(1))
    max_rt <- max(abs(q_rt - pos), na.rm = TRUE)
  }
  list(max_p = max(abs(cpp_p - r_p), na.rm = TRUE), max_q = max(abs(cpp_q - r_q), na.rm = TRUE),
       na_match = identical(is.na(cpp_q), is.na(r_q)), max_rt = max_rt)
}

test_that("2-1 to 2-3: zig_p and zig_q match the R formulas and the ECDF fallback; round trip", {
  for (n in ZIG_SIZES) for (fz in ZERO_FRACS) {
    res <- run_zig_scenario(n, fz)
    lab <- sprintf("n=%d fz=%.2f", n, fz)
    expect_lt(res$max_p, TOL, label = paste(lab, "p"))
    expect_lt(res$max_q, TOL, label = paste(lab, "q"))
    expect_true(res$na_match, label = paste(lab, "NA pattern"))
    if (!is.na(res$max_rt)) expect_lt(res$max_rt, TOL, label = paste(lab, "round trip"))
  }
})

test_that("2-4: ZAGA equivalence (gamlss.dist)", {
  skip_if_not_installed("gamlss.dist")
  raw_z <- c(rep(0.0, 6), ZIG_POOL[seq_len(16)])
  par_z <- zig_fit(matrix(raw_z, nrow = 1L))
  sh <- par_z[1L, 1L]; sc <- par_z[1L, 2L]; pz_z <- par_z[1L, 3L]
  mu_z <- sh * sc; sig_z <- 1 / sqrt(sh)
  pos_z <- sort(unique(raw_z[raw_z > 0]))
  expect_lt(max(abs(vapply(pos_z, function(v) r_zig_p_param(sh, sc, pz_z, v), numeric(1)) -
                    gamlss.dist::pZAGA(pos_z, mu = mu_z, sigma = sig_z, nu = pz_z))), TOL)
  probs_pos <- ZIG_PROBS[ZIG_PROBS > pz_z & ZIG_PROBS < 1.0]
  expect_lt(max(abs(r_zig_q(sh, sc, pz_z, probs_pos) -
                    vapply(probs_pos, function(p) gamlss.dist::qZAGA(p, mu = mu_z, sigma = sig_z, nu = pz_z),
                           numeric(1)))), TOL)
})

test_that("2-5: all-zero records", {
  for (n in ZIG_SIZES) {
    mat <- matrix(rep(0.0, n), nrow = 1L); par <- zig_fit(mat)
    expect_equal(zig_p(par, 0.0, mat), 0.5)
    expect_equal(zig_p(par, 10.0, mat), 1.0)
    expect_identical(zig_q(par, 0.0, mat)[1L, 1L], 0)
    expect_identical(zig_q(par, 0.5, mat)[1L, 1L], 0)
    expect_identical(zig_q(par, 1.0, mat)[1L, 1L], 0)
  }
})

test_that("2-6: integer and double input give identical output", {
  for (n in ZIG_SIZES) for (fz in ZERO_FRACS) {
    n_zero <- min(round(fz * n), n); n_pos <- n - n_zero
    pool_i <- as.integer(round(ZIG_POOL[seq_len(n)] * 1))
    raw_i  <- if (n_pos == 0L) integer(n) else as.integer(c(rep(0L, n_zero), pool_i[seq_len(n_pos)]))
    mat_i  <- matrix(raw_i, nrow = 1L); mat_d <- matrix(as.double(raw_i), nrow = 1L)
    par_i  <- zig_fit(mat_i); par_d <- zig_fit(mat_d)
    expect_identical(max(abs(par_i - par_d), na.rm = TRUE), 0)
    obs_i <- c(0L, as.integer(sort(unique(raw_i[raw_i > 0])))); nr <- length(obs_i)
    expect_identical(max(abs(zig_p(par_d[rep(1L, nr), , drop = FALSE], as.double(obs_i), mat_d[rep(1L, nr), , drop = FALSE]) -
                             zig_p(par_i[rep(1L, nr), , drop = FALSE], obs_i, mat_i[rep(1L, nr), , drop = FALSE])),
                         na.rm = TRUE), 0)
    expect_identical(max(abs(zig_q(par_d, ZIG_PROBS, mat_d) - zig_q(par_i, ZIG_PROBS, mat_i)), na.rm = TRUE), 0)
  }
})

test_that("2-7: zig_p(0) = p_zero / 2 as p_zero -> 1", {
  for (nz in c(seq(0L, 90L, 10L), 95L, 99L, 100L)) {
    mat <- matrix(c(rep(0.0, nz), rep(30.0, 100L - nz)), nrow = 1L)
    par <- zig_fit(mat)
    expect_equal(zig_p(par, 0.0, mat), par[1L, 3L] / 2)
  }
})

test_that("2-8: q0000 = 0; q1000 = NA with positives, 0 when all zero", {
  for (n in ZIG_SIZES) for (fz in ZERO_FRACS) {
    n_zero <- min(round(fz * n), n); n_pos <- n - n_zero
    mat <- matrix(c(rep(0.0, n_zero), ZIG_POOL[seq_len(max(n_pos, 0L))]), nrow = 1L)
    par <- zig_fit(mat)
    expect_identical(zig_q(par, 0.0, mat)[1L, 1L], 0)
    q1 <- zig_q(par, 1.0, mat)[1L, 1L]
    if (n_pos == 0L) expect_identical(q1, 0) else expect_true(is.na(q1))
  }
})

# ── Part 3: NA years ──────────────────────────────────────────────────────────

all_outputs <- function(x, obs) {
  mat <- matrix(x, nrow = 1L)
  par <- zig_fit(mat)
  list(ecdf_q = ecdf_q(mat, ZIG_PROBS), ecdf_p = ecdf_p(mat, obs), ecdf_mean = ecdf_mean(mat),
       zig_fit = par, zig_q = zig_q(par, ZIG_PROBS, mat), zig_p = zig_p(par, obs, mat),
       zig_mean = zig_mean(par, mat))
}
is_all_na <- function(out) all(vapply(out, \(v) all(is.na(v)), logical(1)))

test_that("3-1: a row with NA years equals the row with those years dropped", {
  withr_seed(42, {
    rgamma(100, shape = 3, scale = 10)   # the script draws ZIG_POOL first
    for (n in c(22L, 30L)) for (fz in c(0, 0.3, 0.9)) for (n_na in c(1L, 5L)) {
      n_zero <- round(fz * n)
      x      <- sample(c(rep(0.0, n_zero), ZIG_POOL[seq_len(n - n_zero)]))
      na_at  <- sample(n, n_na)
      x_na   <- replace(x, na_at, NA_real_)
      obs    <- c(0, median(ZIG_POOL))
      expect_equal(all_outputs(x_na, obs[1]), all_outputs(x[-na_at], obs[1]))
      expect_equal(all_outputs(x_na, obs[2]), all_outputs(x[-na_at], obs[2]))
    }
  })
})

test_that("3-2: too few valid years -> NA; a complete short record is used", {
  x22 <- ZIG_POOL[seq_len(22)]
  cases <- list(
    list(x = replace(x22, 1:12, NA), na = FALSE),
    list(x = replace(x22, 1:13, NA), na = TRUE),
    list(x = x22[1:5],               na = FALSE),
    list(x = replace(x22[1:5], 1, NA), na = TRUE),
    list(x = rep(NA_real_, 22),      na = TRUE))
  for (cs in cases) {
    out <- all_outputs(cs$x, 10)
    if (cs$na) expect_true(is_all_na(out))
    else expect_false(anyNA(c(out$ecdf_mean, out$zig_fit[1L, 3L], out$ecdf_q[1L, 10L])))
  }
})

# ── Part 4: ZIG chain and ZIG-smooth ──────────────────────────────────────────

test_that("4-1: zig_lmom matches the reference L-moments", {
  rows <- list(ZIG_POOL[1:22], c(0, 0, ZIG_POOL[1:20]), replace(ZIG_POOL[1:22], c(3, 9), NA),
               c(rep(0, 19), ZIG_POOL[1:3]), rep(0, 22), c(rep(0, 18), 5, 5, 5, 5))
  mat4 <- do.call(rbind, rows)
  lm4  <- zig_lmom(mat4)
  expect_equal(unname(lm4), unname(t(vapply(rows, r_lmom, numeric(4)))), tolerance = TOL)
  expect_identical(colnames(lm4), c("n", "z", "l1", "lcv"))
  skip_if_not_installed("lmom")
  L <- lmom::samlmu(ZIG_POOL[1:22])
  expect_equal(unname(lm4[1, 3:4]), unname(c(L[[1]], L[[2]] / L[[1]])), tolerance = TOL)
})

test_that("4-2: zig_fit == zig_params_from_lmom(zig_lmom); shape == zig_shape_from_lcv", {
  mat42 <- withr_seed(7, t(replicate(200, {
    x <- sample(c(rep(0, sample(0:20, 1)), ZIG_POOL))[1:22]
    replace(x, sample(22, sample(0:14, 1)), NA)
  })))
  chain <- zig_params_from_lmom(zig_lmom(mat42), ncol(mat42))
  expect_identical(chain, zig_fit(mat42))
  # And today's one-pass C++ zig_fit(), kept internal for this check.
  expect_identical(chain, normals:::cpp_zig_fit(mat42, 4L, legacy_mv(mat42)))
  lm42 <- zig_lmom(mat42); has <- !is.na(chain[, 1])
  expect_identical(zig_shape_from_lcv(unname(lm42[has, "lcv"])), chain[has, 1])
  expect_equal(chain[has, 1] * chain[has, 2], unname(lm42[has, "l1"]), tolerance = TOL)
})

SPLINE_B <- zig_spline_basis()
arr_stats <- function(arr) {          # arr: n_pix x 365 x n_years
  ss <- lapply(seq_len(dim(arr)[2]), \(d) zig_spline_suffstats(arr[, d, , drop = TRUE] |> matrix(dim(arr)[1])))
  lapply(c(n = "n", z = "z", l1 = "l1", lcv = "lcv"), \(v) matrix(sapply(ss, \(m) m[, v]), dim(arr)[1]))
}
curve <- \(fit, part, i = 1) drop(SPLINE_B$X %*% fit$coef[i, (part - 1) * SPLINE_B$k + seq_len(SPLINE_B$k)])

test_that("4-3: a stationary record gives constant curves at the pooled values", {
  arr <- array(rep(ZIG_POOL[1:22], each = 365), c(1, 365, 22))
  fit <- zig_spline_fit_stats(arr_stats(arr), SPLINE_B)
  l   <- r_lmom(ZIG_POOL[1:22])
  expect_equal(exp(curve(fit, 2)), rep(l[["l1"]], 365), tolerance = 1e-8)
  expect_equal(exp(curve(fit, 3)), rep(l[["lcv"]], 365), tolerance = 1e-8)
  expect_identical(fit$status, 0L)
  expect_lt(plogis(max(curve(fit, 1))), 1e-3)
})

test_that("4-4: the fit recovers known seasonal l1 and L-CV curves", {
  day <- 1:365
  true_l1  <- 200 + 150 * sin(2 * pi * day / 365)
  true_lcv <- 0.25 + 0.1 * cos(2 * pi * day / 365)
  true_sh  <- zig_shape_from_lcv(true_lcv)
  arr <- array(NA_real_, c(1, 365, 60))
  withr_seed(4, for (d in day) arr[1, d, ] <- rgamma(60, shape = true_sh[d], scale = true_l1[d] / true_sh[d]))
  fit <- zig_spline_fit_stats(arr_stats(arr), SPLINE_B)
  expect_identical(fit$status, 0L)
  expect_lt(max(abs(exp(curve(fit, 2)) / true_l1 - 1)), 0.05)
  expect_lt(max(abs(exp(curve(fit, 3)) / true_lcv - 1)), 0.08)
})

test_that("4-5: fit status codes", {
  mk <- function(f) { a <- array(0, c(1, 365, 22)); for (d in 1:365) a[1, d, ] <- f(d); a }
  cases <- list(
    list(arr = mk(\(d) if (d > 60 && d < 270) ZIG_POOL[1:22] else rep(0, 22)), status = 0L),
    list(arr = mk(\(d) rep(0, 22)), status = 1L),
    list(arr = mk(\(d) if (d == 150) c(ZIG_POOL[1:3], rep(0, 19)) else rep(0, 22)), status = 2L),
    list(arr = mk(\(d) if (d %% 30 == 0) c(ZIG_POOL[1:3], rep(0, 19)) else rep(0, 22)), status = 2L),
    list(arr = mk(\(d) rep(NA, 22)), status = 3L),
    list(arr = mk(\(d) if (d >= 150 && d < 160) ZIG_POOL[1:22] else rep(0, 22)), status = 5L))
  for (cs in cases) {
    st <- arr_stats(cs$arr); st$n[is.na(st$n)] <- 0
    fit <- zig_spline_fit_stats(st, SPLINE_B)
    expect_identical(fit$status, cs$status)
    if (cs$status == 5L) {
      expect_lt(diff(range(curve(fit, 2))), 1e-10)
      expect_lt(diff(range(curve(fit, 3))), 1e-10)
    }
    if (cs$status %in% c(1L, 2L, 3L)) expect_true(all(is.na(fit$coef[1, SPLINE_B$k + 1:(2 * SPLINE_B$k)])))
  }
})

test_that("4-6: zig_spline_params: Hosking shape, scale = l1 / shape, clamp, status 1 -> p_zero 1", {
  k <- SPLINE_B$k
  coefs <- rbind(c(rep(-3, k), rep(log(300), k), rep(log(0.2), k)),
                 c(rep(-3, k), rep(log(50), k),  rep(log(5), k)),
                 c(rep(0, k),  rep(NA, 2 * k)))
  pr <- zig_spline_params(coefs, SPLINE_B, status = c(0L, 0L, 1L))
  expect_equal(pr$shape[1, ], rep(zig_shape_from_lcv(0.2), 365))
  expect_equal(pr$shape[1, ] * pr$scale[1, ], rep(300, 365))
  expect_equal(pr$p_zero[1, ], rep(plogis(-3), 365))
  expect_equal(pr$shape[2, ], rep(zig_shape_from_lcv(ZIG_SMOOTH_LCV_RANGE[2]), 365))
  expect_true(all(pr$p_zero[3, ] == 1))
  # smooth_eval() is the same evaluation, one position at a time.
  fit <- structure(list(coef = coefs, status = c(0L, 0L, 1L), k = k, family = legacy_zig),
                   class = "normals_smooth")
  ev <- smooth_eval(fit, 1:365)
  expect_identical(ev, pr)
})
