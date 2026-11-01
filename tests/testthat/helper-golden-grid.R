# The gridded golden comparison, shared by the 20 x 20 fixture and the
# 50 x 50 release asset: mco-snodas calc_snodas_normals(), annual_fit_block()
# and annual_materialize(), as API.md 6.2's matrix-level calls.

grid_setup <- function(fx) {
  cube <- grid_cube(fx$nr, fx$nc)
  expect_identical(grid_checksum(cube), fx$checksum)
  W <- fx$wy
  wy_start <- year_start(W, "wy")
  # P1: the same day in the 30 years before, newest first (as terra stacks
  # the COGs), from SNODAS's first day; the April 2015 COG is missing.
  hist_mat <- function(d) {
    rd <- reference_dates(d, W:(W - 30L), "wy")
    rd <- rd[rd >= as.Date("2003-09-30") & rd < wy_start]
    rd <- rd[!(format(d, "%Y-%m") == "2024-04" & format(rd, "%Y") == "2015")]
    grid_at(cube, rd)
  }
  list(cube = cube, dates = year_dates(W, "wy"), hist_mat = hist_mat,
       mv_of = function(r) min_valid_years(ncol(r), cap = 10))
}

# Everything that does not pass through the spline fit: ECDF and ZIG normals,
# the daily statistics, every flavor's percentiles (the smoothed ones from
# the float32 parameters the pipeline reads, a fixture input, P3), and the
# smoothed mean and quantiles from those parameters (D1a grid: no fallback).
grid_golden_check <- function(fx) {
  g <- grid_setup(fx)
  zig <- family_zig(); emp <- family_empirical(); P <- NORMALS_QUANTILES
  for (ds in names(fx$expected$per_date)) {
    e <- fx$expected$per_date[[ds]]
    d <- as.Date(ds); lab <- ds
    h <- g$hist_mat(d); x <- grid_at(g$cube, d)[, 1]; mv <- g$mv_of(h)
    expect_identical(ncol(h), e$n_layers)
    d_e <- normals_dist(emp, NULL, h, mv, NULL)
    expect_golden(dist_mean(d_e), e$ecdf$mean, label = paste(lab, "ecdf mean"))
    expect_golden(dist_quantile(d_e, P), e$ecdf$quantiles, label = paste(lab, "ecdf quantiles"))
    lm <- day_stats(h, zig)
    expect_golden(lm, e$zig$lmom, label = paste(lab, "lmom"))
    pz <- fit_day(lm, zig, mv)
    expect_golden(unname(pz), e$zig$params, label = paste(lab, "zig params"))
    d_z <- normals_dist(zig, pz, h, mv, "empirical")
    expect_golden(unname(dist_quantile(d_z, P)), e$zig$quantiles, label = paste(lab, "zig quantiles"))
    expect_golden(dist_mean(d_z), e$zig$mean, label = paste(lab, "zig mean"))
    d_s <- normals_dist(zig, e$smooth_f32, NULL, NULL, "none")
    expect_golden(dist_mean(d_s), e$smooth$mean, label = paste(lab, "smooth mean from float32"))
    expect_golden(unname(dist_quantile(d_s, P)), unname(e$smooth$quantiles),
                  label = paste(lab, "smooth quantiles from float32"))
    # Percentiles (nightly): D1b = ECDF fallback, against the float32 params;
    # D4 integer percent, D5 hundredths.
    p <- list(ecdf = dist_cdf(d_e, x),
              zig = dist_cdf(normals_dist(zig, e$zig_f32, h, mv, "empirical"), x),
              zig_smooth = dist_cdf(normals_dist(zig, e$smooth_f32, h, mv, "empirical"), x))
    for (fl in names(p)) {
      expect_golden(as.integer(as_ptile(p[[fl]], 0)), e$ptile[[fl]], label = paste(lab, fl, "ptile"))
      expect_golden(standardize(p[[fl]], "acklam", "hundredths"), e$standardized[[fl]],
                    label = paste(lab, fl, "standardized"))
    }
  }
}

# The spline fit and the smoothed parameters it gives (fixture platform only).
grid_golden_check_smoothed <- function(fx) {
  skip_smoothed_off_platform()
  g <- grid_setup(fx)
  zig <- family_zig()
  fit_dates <- g$dates[day_key(g$dates) != "02-29"]
  stats <- stack_stats(lapply(fit_dates, \(d) day_stats(g$hist_mat(d), zig)))
  fit <- smooth_fit(stats, cycle_day(fit_dates, "wy"), zig)
  ex <- fx$expected$fit
  expect_golden(fit$coef, ex$coef, label = "coef")
  expect_golden(fit$edf, ex$edf, label = "edf")
  expect_golden(fit$iter, unname(ex$iter), label = "iter")   # names: chunk artifacts
  expect_golden(fit$status, unname(ex$status), label = "status")
  expect_true(all(c(0L, 1L, 2L, 3L, 5L) %in% fit$status))
  ev <- smooth_eval(fit, cycle_day(g$dates, "wy"))
  for (ds in names(fx$expected$per_date)) {
    ps <- smooth_params(ev, match(as.Date(ds), g$dates))
    expect_golden(unname(ps), fx$expected$per_date[[ds]]$smooth$params, label = paste(ds, "smooth params"))
  }
}
