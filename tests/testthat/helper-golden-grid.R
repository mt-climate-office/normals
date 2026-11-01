# The gridded golden comparison, shared by the 20 x 20 fixture and the
# 50 x 50 release asset: mco-snodas calc_snodas_normals(), annual_fit_block()
# and annual_materialize(), as API.md 6.2's matrix-level calls.

grid_golden_check <- function(fx) {
  cube <- grid_cube(fx$nr, fx$nc)
  expect_identical(grid_checksum(cube), fx$checksum)
  zig <- family_zig()
  emp <- family_empirical()
  P   <- NORMALS_QUANTILES
  W   <- fx$wy
  dates <- year_dates(W, "wy")
  wy_start <- year_start(W, "wy")
  # P1: the same day in the 30 years before, newest first (as terra stacks
  # the COGs), from SNODAS's first day; the April 2015 COG is missing.
  hist_mat <- function(d) {
    rd <- reference_dates(d, W:(W - 30L), "wy")
    rd <- rd[rd >= as.Date("2003-09-30") & rd < wy_start]
    rd <- rd[!(format(d, "%Y-%m") == "2024-04" & format(rd, "%Y") == "2015")]
    grid_at(cube, rd)
  }
  mv_of <- function(r) min_valid_years(ncol(r), cap = 10)

  fit_dates <- dates[day_key(dates) != "02-29"]
  stats <- stack_stats(lapply(fit_dates, \(d) day_stats(hist_mat(d), zig)))
  fit <- smooth_fit(stats, cycle_day(fit_dates, "wy"), zig)
  ex <- fx$expected$fit
  expect_golden(fit$coef, ex$coef, label = "coef")
  expect_golden(fit$edf, ex$edf, label = "edf")
  expect_golden(fit$iter, unname(ex$iter), label = "iter")   # names: chunk artifacts
  expect_golden(fit$status, unname(ex$status), label = "status")
  expect_true(all(c(0L, 1L, 2L, 3L, 5L) %in% fit$status))
  ev <- smooth_eval(fit, cycle_day(dates, "wy"))

  for (ds in names(fx$expected$per_date)) {
    e <- fx$expected$per_date[[ds]]
    d <- as.Date(ds); lab <- ds
    h <- hist_mat(d); x <- grid_at(cube, d)[, 1]; mv <- mv_of(h)
    expect_identical(ncol(h), e$n_layers)
    # ECDF
    d_e <- normals_dist(emp, NULL, h, mv, NULL)
    expect_golden(dist_mean(d_e), e$ecdf$mean, label = paste(lab, "ecdf mean"))
    expect_golden(dist_quantile(d_e, P), e$ecdf$quantiles, label = paste(lab, "ecdf quantiles"))
    # ZIG (annual: double parameters)
    lm <- day_stats(h, zig)
    expect_golden(lm, e$zig$lmom, label = paste(lab, "lmom"))
    pz <- fit_day(lm, zig, mv)
    expect_golden(unname(pz), e$zig$params, label = paste(lab, "zig params"))
    d_z <- normals_dist(zig, pz, h, mv, "empirical")
    expect_golden(unname(dist_quantile(d_z, P)), e$zig$quantiles, label = paste(lab, "zig quantiles"))
    expect_golden(dist_mean(d_z), e$zig$mean, label = paste(lab, "zig mean"))
    # ZIG-smooth: the date's parameters, then (P3) the float32 tile the
    # pipeline reads; D1a grid = no fallback for the quantiles and mean.
    ps <- smooth_params(ev, match(d, dates))
    expect_golden(unname(ps), e$smooth$params, label = paste(lab, "smooth params"))
    d_s <- normals_dist(zig, e$smooth_f32, NULL, NULL, "none")
    expect_golden(dist_mean(d_s), e$smooth$mean, label = paste(lab, "smooth mean"))
    expect_golden(unname(dist_quantile(d_s, P)), unname(e$smooth$quantiles), label = paste(lab, "smooth quantiles"))
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
