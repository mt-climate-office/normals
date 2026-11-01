#' The outlook: the likely value for the rest of the year
#'
#' Today's value is placed in its day's distribution (z0, clamped), and that
#' percentile is carried forward on the standard-normal scale, fading toward
#' the normals:
#'
#' * rho(h) = exp(-(h / length)^alpha), h days ahead (`outlook_persistence()`);
#' * k = min(exp(-beta rho min(z0, 0)), k_max) widens the spread below normal
#'   (`outlook_spread()`);
#' * q_p = Q_target(pnorm(rho z0 + qnorm(p) k sqrt(1 - rho^2))) x anchor^rho
#'   (`outlook_quantiles()`), with the anchor scaling a value beyond z = +/-4
#'   (`outlook_start()`).
#'
#' z is the exact `stats::qnorm()` of the percentile, unrounded, clamped to
#' +/- `z_max` (not the Acklam approximation of the published standardized
#' values).
#'
#' For zero-mass families whose parameters include `p_zero`, a p_zero below
#' `p_zero_floor` / n_valid is taken as zero for subjects with a positive value
#' today (a smoothed fit never reaches exactly zero, and a record-low subject
#' would otherwise be placed on that sliver and read as zero). A subject at
#' zero today keeps the fitted chance, so it starts in the zero mass.
#' `outlook_p_zero()` gives the chance of exactly zero on each target date.
#'
#' The alpha table (persistence shape by cycle day) and beta (spread) are per
#' dataset and are passed in: `outlook_alpha_at()` interpolates an alpha table
#' the caller has already filtered to one period and flavor. `outlook_lengths()`
#' fits the persistence length per subject and date from the reference years.
#' `outlook()` composes the pieces for one issue date. `outlook_crps()` and
#' `outlook_pit()` score hindcasts.
#'
#' The issue and target distributions need their reference matrix (`ref`):
#' n_valid for the floor, the empirical fallback, and (in `outlook_lengths()`)
#' the reference years' own values come from it.
#'
#' @param k_max,fit_lags,pool_days,min_pairs,min_lags,z_max,l_grid,p_zero_floor,anchor_range,rho_floor
#'   The method's constants (today's values).
#' @param quantile_digits,p_zero_digits,length_digits,alpha_digits Roundings
#'   of the published columns.
#' @param control An `outlook_control()` list.
#' @param table Data frame with `day` (cycle day) and `alpha`.
#' @param dates,target_dates `Date` vectors.
#' @param year_type `"wy"` or `"cy"`.
#' @param lead Days ahead.
#' @param length Persistence length (days), per subject or scalar.
#' @param alpha Persistence shape, per subject or scalar (`outlook()`), or one
#'   per date (`outlook_lengths()`).
#' @param z0 Starting z.
#' @param rho Persistence.
#' @param beta Spread factor.
#' @param issue The issue date's `normals_dist`.
#' @param target A target date's `normals_dist`.
#' @param targets List of `normals_dist`, one per target date.
#' @param x Values on the issue date, length n.
#' @param start An `outlook_start()` result.
#' @param probs Quantile probabilities.
#' @param issue_date The issue date.
#' @param subjects Data frame of subject key columns, one row per subject.
#' @param dists List of `normals_dist`, one per date of the normals year.
#' @param q n x m quantiles at `probs`.
#' @param y Outcomes, length n.
#' @param u Uniform draws for zero outcomes, length n.
#' @return `outlook_control()`: a list. `outlook_start()`: `list(z0, anchor,
#'   floor)`. `outlook_quantiles()`: an n x m matrix, unrounded. `outlook()`:
#'   a tibble (`date`, subject columns, `q0000` ... `q1000`, `p_zero`), or
#'   `NULL` when there is no target date or no value today.
#'   `outlook_lengths()`: a tibble (`date`, subject columns, `outlook_length`,
#'   `outlook_alpha`). The others: double vectors.
#' @name outlook
NULL

#' @rdname outlook
#' @export
outlook_control <- function(k_max = 3,
                            fit_lags = c(3, 7, 14, 21, 30, 45, 60, 90, 120),
                            pool_days = 15L, min_pairs = 8L, min_lags = 3L,
                            z_max = 4,
                            l_grid = exp(seq(log(2), log(3000), length.out = 400)),
                            p_zero_floor = 0.5, anchor_range = c(0.1, 10),
                            rho_floor = 1e-12, quantile_digits = 1,
                            p_zero_digits = 3, length_digits = 2,
                            alpha_digits = 3) {
  num1 <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      stop("normals: `", nm, "` must be a single number", call. = FALSE)
    }
  }
  for (nm in c("k_max", "pool_days", "min_pairs", "min_lags", "z_max",
               "p_zero_floor", "rho_floor", "quantile_digits", "p_zero_digits",
               "length_digits", "alpha_digits")) {
    num1(get(nm), nm)
  }
  if (!is.numeric(fit_lags) || !length(fit_lags) || anyNA(fit_lags)) {
    stop("normals: `fit_lags` must be a numeric vector", call. = FALSE)
  }
  if (!is.numeric(l_grid) || !length(l_grid) || anyNA(l_grid)) {
    stop("normals: `l_grid` must be a numeric vector", call. = FALSE)
  }
  if (!is.numeric(anchor_range) || length(anchor_range) != 2L || anyNA(anchor_range)) {
    stop("normals: `anchor_range` must be two numbers", call. = FALSE)
  }
  list(k_max = k_max, fit_lags = fit_lags, pool_days = pool_days,
       min_pairs = min_pairs, min_lags = min_lags, z_max = z_max,
       l_grid = l_grid, p_zero_floor = p_zero_floor,
       anchor_range = anchor_range, rho_floor = rho_floor,
       quantile_digits = quantile_digits, p_zero_digits = p_zero_digits,
       length_digits = length_digits, alpha_digits = alpha_digits)
}

#' @rdname outlook
#' @export
outlook_alpha_at <- function(table, dates, year_type) {
  if (!is.data.frame(table) || !all(c("day", "alpha") %in% names(table))) {
    stop("normals: `table` must be a data frame with columns `day` and `alpha`", call. = FALSE)
  }
  if (!nrow(table)) stop("normals: `table` has no rows", call. = FALSE)
  if (!inherits(dates, "Date")) stop("normals: `dates` must be a Date vector", call. = FALSE)
  table <- table[order(table$day), ]
  stats::approx(table$day, table$alpha, xout = cycle_day(dates, year_type), rule = 2)$y
}

#' @rdname outlook
#' @export
outlook_persistence <- function(lead, length, alpha) {
  exp(-(lead / length)^alpha)
}

#' @rdname outlook
#' @export
outlook_spread <- function(z0, rho, beta, control = outlook_control()) {
  pmin(exp(-beta * rho * pmin(z0, 0)), control$k_max)
}

#' @rdname outlook
#' @export
outlook_start <- function(issue, x, control = outlook_control()) {
  check_outlook_dist(issue, "issue")
  n <- dist_rows(issue)
  if (!is.numeric(x) || length(x) != n) {
    stop("normals: `x` must be numeric with one value per row of `issue` (", n, ")", call. = FALSE)
  }
  floor <- !is.na(x) & x > 0
  d0 <- outlook_floored(issue, floor, control)
  z0 <- outlook_z(dist_cdf(d0, x), control)
  # The ratio of today's value to the distribution's quantile at today's
  # (clamped) percentile: about 1 within the distribution, beyond it when the
  # value is past z = +/- z_max. Kept on the side of today's anomaly (a quantile
  # that underflows at -z_max must not scale a low subject up) and within
  # `anchor_range`; 1 where either is zero or missing.
  u0 <- stats::pnorm(z0)
  q0 <- outlook_q(d0, matrix(u0, ncol = 1))[, 1]
  r  <- x / q0
  r  <- ifelse(is.finite(r) & r > 0 & !is.na(z0), r, 1)
  r  <- ifelse(z0 > 0, pmax(r, 1), pmin(r, 1))
  anchor <- pmin(pmax(r, control$anchor_range[1]), control$anchor_range[2])
  list(z0 = z0, anchor = anchor, floor = floor)
}

#' @rdname outlook
#' @export
outlook_quantiles <- function(target, start, rho, beta,
                              probs = NORMALS_QUANTILES,
                              control = outlook_control()) {
  check_outlook_dist(target, "target")
  check_outlook_start(start, dist_rows(target))
  zq <- outlook_zq(probs, control)
  sg <- outlook_sigma(start$z0, rho, beta, control)
  u  <- stats::pnorm(rho * start$z0 + outer(sg, zq))
  q  <- outlook_q(outlook_floored(target, start$floor, control), u) * start$anchor^rho
  colnames(q) <- quantile_columns(probs)
  q
}

#' @rdname outlook
#' @export
outlook_p_zero <- function(target, start, rho, beta,
                           control = outlook_control()) {
  check_outlook_dist(target, "target")
  if (!isTRUE(target$family$zero_mass)) {
    stop("normals: `target`'s family has no zero mass", call. = FALSE)
  }
  check_outlook_start(start, dist_rows(target))
  sg <- outlook_sigma(start$z0, rho, beta, control)
  p0 <- dist_p_zero(outlook_floored(target, start$floor, control))
  ifelse(p0 > 0, stats::pnorm((stats::qnorm(p0) - rho * start$z0) / sg), 0)
}

#' @rdname outlook
#' @export
outlook <- function(issue, targets, x, issue_date, target_dates, subjects,
                    length, alpha, beta, probs = NORMALS_QUANTILES,
                    control = outlook_control()) {
  check_outlook_dist(issue, "issue")
  n <- dist_rows(issue)
  if (!inherits(issue_date, "Date") || base::length(issue_date) != 1L || is.na(issue_date)) {
    stop("normals: `issue_date` must be a single Date", call. = FALSE)
  }
  if (!inherits(target_dates, "Date")) {
    stop("normals: `target_dates` must be a Date vector", call. = FALSE)
  }
  if (!is.list(targets) || inherits(targets, "normals_dist") ||
      base::length(targets) != base::length(target_dates)) {
    stop("normals: `targets` must be a list of normals_dist, one per `target_dates`", call. = FALSE)
  }
  if (!is.data.frame(subjects) || nrow(subjects) != n) {
    stop("normals: `subjects` must be a data frame with one row per row of `issue` (", n, ")",
         call. = FALSE)
  }
  if (!is.numeric(length) || !base::length(length) %in% c(1L, n)) {
    stop("normals: `length` must be numeric, length 1 or ", n, call. = FALSE)
  }
  if (!is.numeric(alpha) || !base::length(alpha) %in% c(1L, n)) {
    stop("normals: `alpha` must be numeric, length 1 or ", n, call. = FALSE)
  }
  if (!is.numeric(beta) || base::length(beta) != 1L || is.na(beta)) {
    stop("normals: `beta` must be a single number", call. = FALSE)
  }
  if (!base::length(target_dates) || all(is.na(x))) return(NULL)
  for (t in seq_along(targets)) {
    check_outlook_dist(targets[[t]], "targets")
    if (dist_rows(targets[[t]]) != n) {
      stop("normals: every element of `targets` must have `issue`'s ", n, " rows", call. = FALSE)
    }
  }
  zero_mass <- isTRUE(issue$family$zero_mass)
  length <- rep_len(length, n); alpha <- rep_len(alpha, n)

  start <- outlook_start(issue, x, control)
  keep  <- !is.na(start$z0) & !is.na(length)
  nt <- base::length(target_dates)
  q  <- vector("list", nt); p0 <- vector("list", nt)
  for (t in seq_len(nt)) {
    rho <- outlook_persistence(as.numeric(target_dates[t] - issue_date), length, alpha)
    q[[t]] <- round(outlook_quantiles(targets[[t]], start, rho, beta, probs, control),
                    control$quantile_digits)[keep, , drop = FALSE]
    if (zero_mass) {
      p0[[t]] <- round(outlook_p_zero(targets[[t]], start, rho, beta, control),
                       control$p_zero_digits)[keep]
    }
  }
  k  <- sum(keep)
  ix <- rep(which(keep), nt)
  out <- c(list(date = rep(target_dates, each = k)),
           as.list(subjects[ix, , drop = FALSE]),
           as.list(as.data.frame(do.call(rbind, q))))
  if (zero_mass) out$p_zero <- unlist(p0, use.names = FALSE)
  out <- tibble::as_tibble(out)
  ord <- do.call(order, c(unname(as.list(subjects[ix, , drop = FALSE])),
                          list(out$date, method = "radix")))
  out[ord, ]
}

#' @rdname outlook
#' @export
outlook_lengths <- function(dists, dates, subjects, alpha,
                            control = outlook_control()) {
  if (!is.list(dists) || inherits(dists, "normals_dist") || !base::length(dists)) {
    stop("normals: `dists` must be a list of normals_dist", call. = FALSE)
  }
  if (!inherits(dates, "Date") || base::length(dates) != base::length(dists)) {
    stop("normals: `dates` must be a Date vector, one per element of `dists`", call. = FALSE)
  }
  if (base::length(dates) > 1L && !all(diff(as.numeric(dates)) == 1)) {
    stop("normals: `dates` must be consecutive days", call. = FALSE)
  }
  for (d in dists) check_outlook_dist(d, "dists")
  n <- dist_rows(dists[[1]])
  if (any(vapply(dists, dist_rows, 1L) != n)) {
    stop("normals: every element of `dists` must have the same rows", call. = FALSE)
  }
  if (!is.data.frame(subjects) || nrow(subjects) != n) {
    stop("normals: `subjects` must be a data frame with one row per row of `dists` (", n, ")",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || !base::length(alpha) %in% c(1L, base::length(dates))) {
    stop("normals: `alpha` must be numeric, length 1 or one per date", call. = FALSE)
  }
  L <- base::length(dates)
  alpha <- rep_len(alpha, L)
  cr <- outlook_correlations(dists, control)
  len <- outlook_fit_lengths(cr, alpha, control)
  out <- c(list(date = rep(dates, each = n)),
           as.list(subjects[rep(seq_len(n), L), , drop = FALSE]),
           list(outlook_length = round(as.vector(len), control$length_digits),
                outlook_alpha = rep(round(alpha, control$alpha_digits), each = n)))
  tibble::as_tibble(out)
}

#' @rdname outlook
#' @export
outlook_crps <- function(q, y, probs) {
  if (!is.matrix(q) || ncol(q) != base::length(probs) || nrow(q) != base::length(y)) {
    stop("normals: `q` must be a length(y) x length(probs) matrix", call. = FALSE)
  }
  2 * rowMeans(((y < q) - matrix(probs, base::length(y), base::length(probs), byrow = TRUE)) * (q - y))
}

#' @rdname outlook
#' @export
outlook_pit <- function(target, start, rho, beta, y, u,
                        control = outlook_control()) {
  check_outlook_dist(target, "target")
  check_outlook_start(start, dist_rows(target))
  sg <- outlook_sigma(start$z0, rho, beta, control)
  dt <- outlook_floored(target, start$floor, control)
  # On the anchored scale: the outcome divided by the anchor factor.
  ya <- y / start$anchor^rho
  p0 <- if (isTRUE(target$family$zero_mass)) dist_p_zero(dt) else rep(0, dist_rows(target))
  ifelse(ya > 0, stats::pnorm((stats::qnorm(dist_cdf(dt, ya)) - rho * start$z0) / sg),
         u * ifelse(p0 > 0, stats::pnorm((stats::qnorm(p0) - rho * start$z0) / sg), 0))
}


# Internal ----------------------------------------------------------------

check_outlook_dist <- function(d, arg) {
  if (!inherits(d, "normals_dist")) {
    stop("normals: `", arg, "` must be a normals_dist", call. = FALSE)
  }
}

check_outlook_start <- function(start, n) {
  if (!is.list(start) || !all(c("z0", "anchor", "floor") %in% names(start)) ||
      base::length(start$z0) != n) {
    stop("normals: `start` must be an outlook_start() result for the same rows", call. = FALSE)
  }
}

dist_rows <- function(d) base::length(dist_n_valid(d))

# z of a percentile: exact qnorm, clamped (API D6).
outlook_z <- function(p, control) {
  pmin(pmax(stats::qnorm(p), -control$z_max), control$z_max)
}

# z of the quantile probabilities; 0 and 1 are read at -z_max and z_max.
outlook_zq <- function(probs, control) {
  stats::qnorm(pmin(pmax(unname(probs), stats::pnorm(-control$z_max)),
                    stats::pnorm(control$z_max)))
}

outlook_sigma <- function(z0, rho, beta, control) {
  sqrt(pmax(control$rho_floor, 1 - rho^2)) * outlook_spread(z0, rho, beta, control)
}

# The distribution with p_zero below `p_zero_floor` / n_valid taken as zero on
# rows where `floor` is TRUE. Only for zero-mass families with a `p_zero`
# parameter; others are returned unchanged.
outlook_floored <- function(d, floor, control) {
  fam <- d$family
  if (!isTRUE(fam$zero_mass) || !"p_zero" %in% fam$params) return(d)
  p  <- d$params
  p0 <- p[, "p_zero"]
  floor <- rep_len(floor, base::length(p0))
  floor[is.na(floor)] <- FALSE
  p0[floor & !is.na(p0) & p0 < control$p_zero_floor / pmax(1L, dist_n_valid(d))] <- 0
  p[, "p_zero"] <- p0
  normals_dist(fam, p, d$ref, d$min_valid, d$fallback)
}

# Quantiles at row-specific probabilities u (n x m), as mco-snotel's outlook
# reads them (outlook:R/outlook.R:113-130), which differs from
# dist_quantile() in two corners:
# * a row that falls back to the empirical distribution with no reference year
#   is NA (the parametric rows never read the record);
# * on a parametric row, a probability past the top of the distribution
#   (dist_quantile() NA; for ZIG, (u - p_zero) / (1 - p_zero) >= 1 in floating
#   point) is +Inf, today's qgamma(1). For the published q1000 (read at
#   z = z_max) this happens only with a large beta and a start near -z_max:
#   rho z0 + z_max sigma k must reach about 8.3, which beta = 0.25 cannot.
outlook_q <- function(d, u) {
  bad <- is.na(u)
  u[bad] <- 0.5
  q <- unname(dist_quantile(d, u))
  n <- nrow(q)
  if (base::length(d$family$params)) {
    complete <- rowSums(is.na(d$params)) == 0
    par <- complete
    if (isTRUE(d$family$zero_mass)) par <- par & dist_p_zero(d) < 1
    q[is.na(q) & par[row(q)] & !bad] <- Inf
    empty <- !complete & dist_n_valid(d) == 0L
  } else {
    empty <- dist_n_valid(d) == 0L
  }
  q[bad] <- NA_real_
  q[rep_len(empty, n), ] <- NA_real_
  q
}

# Lag correlations (n x dates x lags) of the reference years' z: each
# reference year's value placed in its date's distribution (z clamped), a
# reference year at zero in the zero mass (as a subject at zero today is), and
# z's correlation across years between each date and `fit_lags` days later.
outlook_correlations <- function(dists, control) {
  Z <- lapply(dists, \(d) {
    h <- d$ref
    if (is.null(h)) stop("normals: every element of `dists` needs its `ref`", call. = FALSE)
    df <- outlook_floored(d, TRUE, control)
    z <- matrix(NA_real_, nrow(h), ncol(h))
    for (y in seq_len(ncol(h))) {
      fl <- !(h[, y] %in% 0)
      p  <- ifelse(fl, dist_cdf(df, h[, y]), dist_cdf(d, h[, y]))
      z[, y] <- stats::qnorm(p)
    }
    z[is.na(h)] <- NA
    pmin(pmax(z, -control$z_max), control$z_max)
  })
  L <- base::length(dists); lags <- control$fit_lags
  cr <- array(NA_real_, c(nrow(Z[[1]]), L, base::length(lags)))
  for (k in seq_along(lags)) for (t in seq_len(L - lags[k])) {
    cr[, t, k] <- outlook_rowcor(Z[[t]], Z[[t + lags[k]]], control)
  }
  cr
}

# Row-wise correlation of two n x y matrices over the years both have; NA with
# fewer than `min_pairs` years or no spread.
outlook_rowcor <- function(a, b, control) {
  ok <- !is.na(a) & !is.na(b); a[!ok] <- NA; b[!ok] <- NA
  n  <- rowSums(ok)
  da <- a - rowMeans(a, na.rm = TRUE); db <- b - rowMeans(b, na.rm = TRUE)
  sa <- sqrt(rowSums(da^2, na.rm = TRUE)); sb <- sqrt(rowSums(db^2, na.rm = TRUE))
  r  <- rowSums(da * db, na.rm = TRUE) / (sa * sb)
  r[n < control$min_pairs | sa < 1e-6 | sb < 1e-6] <- NA
  r
}

# Persistence lengths from lag correlations cr (n x dates x lags): pooled over
# dates within `pool_days`, fitted by exp(-(h / l)^alpha) by least squares on
# `l_grid`. Dates with fewer than `min_lags` lags take the subject's median.
outlook_fit_lengths <- function(cr, alpha, control) {
  L <- dim(cr)[2]; alpha <- rep_len(alpha, L)
  grid <- control$l_grid; lags <- control$fit_lags
  len <- matrix(NA_real_, dim(cr)[1], L)
  for (t in seq_len(L)) {
    w <- max(1L, t - control$pool_days):min(L, t + control$pool_days)
    r <- rowMeans(aperm(cr[, w, , drop = FALSE], c(1, 3, 2)), na.rm = TRUE, dims = 2)
    r[!is.finite(r)] <- NA
    kern <- exp(-outer(grid, lags, \(l, h) (h / l)^alpha[t]))
    sse <- matrix(0, nrow(r), base::length(grid))
    for (k in seq_along(lags)) {
      d <- outer(r[, k], kern[, k], `-`)^2
      d[is.na(d)] <- 0
      sse <- sse + d
    }
    fit <- grid[max.col(-sse, ties.method = "first")]
    fit[rowSums(!is.na(r)) < control$min_lags] <- NA
    len[, t] <- fit
  }
  med <- apply(len, 1, stats::median, na.rm = TRUE)
  len[is.na(len)] <- med[row(len)[is.na(len)]]
  len
}
