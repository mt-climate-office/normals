# mco-snodas zonal normals (calc_zonal_normals(), annual_fit_zonal()):
# API.md 6.2, the matrix-level calls, per G1 date and level.

# A raw table's reference columns (before the water year) and its latest
# current-year column (NULL if none), as calc_zonal_normals() splits them.
zonal_hist <- function(x, wy_start) {
  dts  <- as.Date(names(x), format = "%Y-%m-%d", optional = TRUE)
  cols <- names(x)[!is.na(dts)]
  dts  <- dts[!is.na(dts)]
  cur  <- cols[dts >= wy_start]
  list(hist = as.matrix(x[, cols[dts < wy_start], drop = FALSE]),
       current = if (length(cur)) as.numeric(x[[cur[length(cur)]]]))
}

zonal_setup <- function() {
  fx <- golden("zonal")
  dates <- year_dates(fx$wy, "wy")
  list(fx = fx, zig = family_zig(), emp = family_empirical(), P = NORMALS_QUANTILES,
       mv_of = function(r) min_valid_years(ncol(r), cap = 10),     # D3, mco-snodas
       dates = dates, wy_start = year_start(fx$wy, "wy"),
       fit_dates = dates[day_key(dates) != "02-29"],
       levels = c(zone = "swe", hypsomeswe = "hypsomeswe"))
}

# One flavor's table for one date and level: swe, ptile, standardized, mean,
# quantiles and (for ZIG flavors) the parameters.
zonal_check <- function(ex, dist, params, x, keys, lab, P) {
  p <- dist_cdf(dist, x)
  expect_golden(x, ex$swe, label = paste(lab, "swe"))
  expect_golden(as_ptile(p, 2), ex$swe_ptile, label = paste(lab, "ptile"))
  expect_golden(standardize(p, "acklam", "z"), ex$swe_standardized, label = paste(lab, "standardized"))
  expect_golden(dist_mean(dist), ex$mean, label = paste(lab, "mean"))
  expect_golden(dist_quantile(dist, P), as.matrix(ex[names(P)]), label = paste(lab, "quantiles"))
  if (!is.null(params)) {
    expect_golden(unname(params[, "shape"]), ex$zig_shape, label = paste(lab, "shape"))
    expect_golden(unname(params[, "scale"]), ex$zig_scale, label = paste(lab, "scale"))
    expect_golden(unname(params[, "p_zero"]), ex$zig_p_zero, label = paste(lab, "p_zero"))
  }
  expect_identical(as.character(ex$zone), keys)
}

test_that("zonal ECDF and ZIG normals, percentiles and daily statistics match mco-snodas", {
  z <- zonal_setup(); fx <- z$fx
  stats_swe <- stack_stats(lapply(z$fit_dates, \(d) day_stats(zonal_hist(fx$raw_swe[[format(d)]], z$wy_start)$hist, z$zig)))
  expect_golden(stats_swe, fx$expected$stats_swe, label = "zone daily statistics")
  for (ds in names(fx$expected$tables)) {
    ex_d <- fx$expected$tables[[ds]]
    raw <- list(zone = fx$raw_swe[[ds]], hypsomeswe = fx$raw_hyps[[ds]])
    for (lv in names(z$levels)) {
      hc <- zonal_hist(raw[[lv]], z$wy_start); h <- hc$hist; mv <- z$mv_of(h)
      keys <- as.character(raw[[lv]]$zone)
      zonal_check(ex_d$ecdf[[lv]], normals_dist(z$emp, NULL, h, mv, NULL), NULL, hc$current, keys,
                  paste(ds, lv, "ecdf"), z$P)
      pz <- fit_day(day_stats(h, z$zig), z$zig, mv)
      zonal_check(ex_d$zig[[lv]], normals_dist(z$zig, pz, h, mv, "empirical"), pz, hc$current, keys,
                  paste(ds, lv, "zig"), z$P)
      # The smoothed flavor's percentiles from its stored parameters (today's
      # nightly path; no linear algebra).
      es <- ex_d$zig_smooth[[lv]]
      ps <- cbind(shape = es$zig_shape, scale = es$zig_scale, p_zero = es$zig_p_zero)
      p <- dist_cdf(normals_dist(z$zig, ps, h, mv, "empirical"), hc$current)
      expect_golden(as_ptile(p, 2), es$swe_ptile, label = paste(ds, lv, "zig_smooth ptile"))
      expect_golden(standardize(p, "acklam", "z"), es$swe_standardized,
                    label = paste(ds, lv, "zig_smooth standardized"))
    }
  }
})

test_that("zonal ZIG-smooth fits and normals match mco-snodas", {
  skip_smoothed_off_platform()
  z <- zonal_setup(); fx <- z$fx
  stats_swe <- stack_stats(lapply(z$fit_dates, \(d) day_stats(zonal_hist(fx$raw_swe[[format(d)]], z$wy_start)$hist, z$zig)))
  fits <- list(zone = smooth_fit(stats_swe, cycle_day(z$fit_dates, "wy"), z$zig),
               hypsomeswe = smooth_fit(fx$stats_hyps, cycle_day(z$fit_dates, "wy"), z$zig))
  for (lv in names(z$levels)) {
    ex <- fx$expected$fit[[z$levels[[lv]]]]
    got <- fits[[lv]]
    expect_golden(got$coef, ex$coef, label = paste(lv, "coef"))
    expect_golden(got$edf, ex$edf, label = paste(lv, "edf"))
    expect_golden(got$iter, ex$iter, label = paste(lv, "iter"))
    expect_golden(got$status, ex$status, label = paste(lv, "status"))
  }
  evals <- lapply(fits, \(f) smooth_eval(f, cycle_day(z$dates, "wy")))
  for (ds in names(fx$expected$tables)) {
    j <- match(as.Date(ds), z$dates)
    raw <- list(zone = fx$raw_swe[[ds]], hypsomeswe = fx$raw_hyps[[ds]])
    for (lv in names(z$levels)) {
      hc <- zonal_hist(raw[[lv]], z$wy_start); h <- hc$hist; mv <- z$mv_of(h)
      # D1a zonal = "empirical", D2 = no mask.
      ps <- smooth_params(evals[[lv]], j)
      zonal_check(fx$expected$tables[[ds]]$zig_smooth[[lv]], normals_dist(z$zig, ps, h, mv, "empirical"),
                  ps, hc$current, as.character(raw[[lv]]$zone), paste(ds, lv, "zig_smooth"), z$P)
    }
  }
})
