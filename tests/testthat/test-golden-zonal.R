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

test_that("zonal ECDF, ZIG and ZIG-smooth normals and percentiles match mco-snodas", {
  fx  <- golden("zonal")
  zig <- family_zig()
  emp <- family_empirical()
  P   <- NORMALS_QUANTILES
  mv_of <- function(r) min_valid_years(ncol(r), cap = 10)    # D3, mco-snodas
  dates <- year_dates(fx$wy, "wy")
  wy_start <- year_start(fx$wy, "wy")
  fit_dates <- dates[day_key(dates) != "02-29"]
  levels <- c(zone = "swe", hypsomeswe = "hypsomeswe")

  # The daily statistics and the smoothed fit, each level.
  stats_swe <- stack_stats(lapply(fit_dates, \(d) day_stats(zonal_hist(fx$raw_swe[[format(d)]], wy_start)$hist, zig)))
  expect_golden(stats_swe, fx$expected$stats_swe, label = "zone daily statistics")
  fits <- list(zone = smooth_fit(stats_swe, cycle_day(fit_dates, "wy"), zig),
               hypsomeswe = smooth_fit(fx$stats_hyps, cycle_day(fit_dates, "wy"), zig))
  for (lv in names(levels)) {
    ex <- fx$expected$fit[[levels[[lv]]]]
    got <- fits[[lv]]
    expect_golden(got$coef, ex$coef, label = paste(lv, "coef"))
    expect_golden(got$edf, ex$edf, label = paste(lv, "edf"))
    expect_golden(got$iter, ex$iter, label = paste(lv, "iter"))
    expect_golden(got$status, ex$status, label = paste(lv, "status"))
  }
  evals <- lapply(fits, \(f) smooth_eval(f, cycle_day(dates, "wy")))

  for (ds in names(fx$expected$tables)) {
    d <- as.Date(ds); j <- match(d, dates)
    ex_d <- fx$expected$tables[[ds]]
    raw <- list(zone = fx$raw_swe[[ds]], hypsomeswe = fx$raw_hyps[[ds]])
    for (lv in names(levels)) {
      lab <- paste(ds, lv)
      hc <- zonal_hist(raw[[lv]], wy_start); h <- hc$hist; x <- hc$current
      mv <- mv_of(h)
      check <- function(ex, dist, params, flavor) {
        p <- dist_cdf(dist, x)
        expect_golden(x, ex$swe, label = paste(lab, flavor, "swe"))
        expect_golden(as_ptile(p, 2), ex$swe_ptile, label = paste(lab, flavor, "ptile"))
        expect_golden(standardize(p, "acklam", "z"), ex$swe_standardized, label = paste(lab, flavor, "standardized"))
        expect_golden(dist_mean(dist), ex$mean, label = paste(lab, flavor, "mean"))
        q <- dist_quantile(dist, P)
        expect_golden(q, as.matrix(ex[names(P)]), label = paste(lab, flavor, "quantiles"))
        if (!is.null(params)) {
          expect_golden(unname(params[, "shape"]), ex$zig_shape, label = paste(lab, flavor, "shape"))
          expect_golden(unname(params[, "scale"]), ex$zig_scale, label = paste(lab, flavor, "scale"))
          expect_golden(unname(params[, "p_zero"]), ex$zig_p_zero, label = paste(lab, flavor, "p_zero"))
        }
        expect_identical(as.character(ex$zone), as.character(raw[[lv]]$zone))
      }
      # ECDF
      check(ex_d$ecdf[[lv]], normals_dist(emp, NULL, h, mv, NULL), NULL, "ecdf")
      # ZIG: the parameters serve the quantiles, mean and percentiles.
      pz <- fit_day(day_stats(h, zig), zig, mv)
      check(ex_d$zig[[lv]], normals_dist(zig, pz, h, mv, "empirical"), pz, "zig")
      # ZIG-smooth: D1a zonal = "empirical", D2 = no mask.
      ps <- smooth_params(evals[[lv]], j)
      check(ex_d$zig_smooth[[lv]], normals_dist(zig, ps, h, mv, "empirical"), ps, "zig_smooth")
    }
  }
})
