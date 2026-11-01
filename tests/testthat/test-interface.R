# The interface (API.md 7, item 3): required arguments, shapes, plug-in
# families, and that the smoothed fit does not depend on `cores`.

toy_ref <- function(n = 4, y = 20, seed = 1) {
  withr_seed(seed, matrix(round(pmax(0, rgamma(n * y, 2, scale = 50) - 40)), n, y,
                          dimnames = list(NULL, 2001:(2000 + y))))
}

test_that("required arguments have no defaults", {
  ref <- toy_ref(); zig <- family_zig()
  expect_error(normals_dist(zig, NULL, ref, 10))
  expect_error(fit_day(day_stats(ref, zig), zig))
  expect_error(standardize(0.5, "acklam"))
  expect_error(as_ptile(0.5))
  expect_error(normals_table(ref = list(), dates = year_dates(2024, "wy"),
                             subjects = data.frame(id = 1:4), family = zig, year_type = "wy",
                             min_valid = 10, fallback = c(fitted = "empirical", smoothed = "empirical")))
})

test_that("arguments are checked, naming the argument", {
  ref <- toy_ref(); zig <- family_zig()
  p <- fit_day(day_stats(ref, zig), zig, 10)
  expect_error(normals_dist(zig, p, ref[1:3, ], 10, "empirical"), "`ref`")
  expect_error(normals_dist(zig, p, NULL, NULL, "empirical"), "`ref`")
  expect_error(normals_dist(zig, p[, 1:2], ref, 10, "empirical"), "`params`")
  expect_error(normals_dist(zig, p, ref, 0, "empirical"), "`min_valid`")
  expect_error(normals_dist(zig, p, ref, 10, "maybe"), "`fallback`")
  expect_error(normals_dist(family_empirical(), p, ref, 10, NULL), "`params`")
  d <- normals_dist(zig, p, ref, 10, "empirical")
  expect_named(d, c("family", "params", "ref", "min_valid", "fallback"))
  expect_identical(colnames(d$params), c("shape", "scale", "p_zero"))
  expect_error(dist_cdf(d, 1:3), "`x`")
  expect_error(dist_quantile(d, matrix(0.5, 3, 2)), "`probs`")
  expect_error(standardize(0.5, "qnorm", "hundredths"), "`as`")
  expect_error(day_stats(ref, family_empirical()), "no daily statistics")
})

test_that("row-specific probabilities equal one row at a time", {
  ref <- toy_ref(n = 6); zig <- family_zig()
  p <- fit_day(day_stats(ref, zig), zig, 10)
  p[2, 1:2] <- NA   # an ECDF fallback row
  u <- withr_seed(3, matrix(runif(6 * 4), 6))
  u[1, 1] <- 1      # p_adj >= 1: NA, as zig_q()
  for (fb in c("empirical", "none")) {
    d <- normals_dist(zig, p, ref, 10, fb)
    q <- dist_quantile(d, u)
    for (i in 1:6) {
      di <- normals_dist(zig, p[i, , drop = FALSE], ref[i, , drop = FALSE], 10, fb)
      expect_identical(q[i, ], unname(dist_quantile(di, u[i, ])[1, ]))
    }
  }
  expect_true(is.na(dist_quantile(normals_dist(zig, p, ref, 10, "empirical"), u)[1, 1]))
  e <- normals_dist(family_empirical(), NULL, ref, 10, NULL)
  q <- dist_quantile(e, u)
  for (i in 1:6) expect_identical(q[i, ], unname(dist_quantile(
    normals_dist(family_empirical(), NULL, ref[i, , drop = FALSE], 10, NULL), u[i, ])[1, ]))
})

test_that("p_zero, n_valid and the sparse mask", {
  ref <- toy_ref(); ref[1, 1:15] <- NA
  zig <- family_zig()
  p <- fit_day(day_stats(ref, zig), zig, 10)
  expect_identical(dist_p_zero(normals_dist(zig, p, ref, 10, "empirical")), unname(p[, 3]))
  e <- normals_dist(family_empirical(), NULL, ref, 10, NULL)
  expect_identical(dist_n_valid(e), c(5L, 20L, 20L, 20L))
  pe <- dist_p_zero(e)
  expect_true(is.na(pe[1]))
  expect_identical(pe[-1], rowSums(ref[-1, ] == 0) / 20)
  expect_true(all(is.na(mask_sparse(p, ref, 10)[1, ])))
  expect_identical(mask_sparse(p, ref, 10)[-1, ], p[-1, ])
})

test_that("standardize and as_ptile", {
  p <- c(NA, 0, 0.025, 0.5, 0.975, 1)
  h <- standardize(p, "acklam", "hundredths")
  expect_type(h, "integer")
  expect_identical(h, c(NA, -400L, -196L, 0L, 196L, 400L))
  expect_identical(standardize(p, "acklam", "z"), h / 100)
  expect_identical(standardize(p, "qnorm", "z"), pmin(pmax(qnorm(p), -4), 4))
  expect_identical(as_ptile(0.123456, 2), 12.35)
})

test_that("records: matrix, lookups, reference matrices and spans", {
  x <- data.frame(id = c("b", "a", "a", "b"), date = as.Date(c("2020-01-01", "2020-01-02", "2021-01-02", "2021-01-01")),
                  v = c(1, 2, 3, 4))
  r <- record_matrix(x, "id", "date", "v")
  expect_s3_class(r, "normals_record")
  expect_identical(r$subjects, c("a", "b"))
  expect_identical(record_at(r, as.Date(c("2019-12-31", "2020-01-02"))), matrix(c(NA, NA, 2, NA), 2))
  h <- reference_record(r, as.Date("2021-01-02"), 2019:2021, "cy")
  expect_identical(colnames(h), c("2019", "2020", "2021"))
  expect_identical(reference_span(h), list(first = c(2020L, NA), last = c(2021L, NA)))
  expect_error(record_matrix(rbind(x, x[1, ]), "id", "date", "v"), "more than one row")
  expect_error(record_matrix(x, "id", "when", "v"), "`date`")
})

# A plug-in family: normal, no zero mass, no smoother (API.md 4).
family_normal <- function() {
  new_family(
    name = "normal", params = c("mu", "sigma"), zero_mass = FALSE, stats_names = c("n", "mean", "sd"),
    stats = function(ref, control) cbind(n = rowSums(!is.na(ref)), mean = rowMeans(ref, na.rm = TRUE),
                                         sd = apply(ref, 1, stats::sd, na.rm = TRUE)),
    fit = function(stats, min_valid, control) {
      out <- stats[, c("mean", "sd"), drop = FALSE]
      colnames(out) <- c("mu", "sigma")
      out[stats[, "n"] < min_valid, ] <- NA
      out
    },
    smooth = NULL,
    cdf = function(params, x, ref, min_valid, fallback, control) stats::pnorm(x, params[, 1], params[, 2]),
    quantile = function(params, probs, ref, min_valid, fallback, control) {
      if (is.matrix(probs)) matrix(stats::qnorm(probs, params[, 1], params[, 2]), nrow(params))
      else vapply(probs, \(p) stats::qnorm(p, params[, 1], params[, 2]), numeric(nrow(params)))
    },
    mean = function(params, ref, min_valid, fallback, control) unname(params[, 1]),
    p_zero = NULL)
}

test_that("a plug-in family validates and runs through both tables", {
  fam <- family_normal()
  expect_true(is_family(fam))
  register_family(fam)
  expect_identical(get_family("normal"), fam)
  expect_error(register_family(new_family("normal", "m", FALSE, "n", function(ref, control) ref,
                                          function(stats, min_valid, control) stats, NULL,
                                          fam$cdf, fam$quantile, fam$mean, NULL)), "already registered")
  expect_error(dist_p_zero(normals_dist(fam, cbind(mu = 1, sigma = 1), NULL, NULL, "none")), "zero mass")
  clash <- fam; clash$params <- c("mean", "sd")
  expect_error(validate_family(clash), "table column")
  bad <- fam; bad$cdf <- function(params, x) x
  expect_error(validate_family(bad), "cdf")

  dates <- year_dates(2021, "cy")
  set <- lapply(dates, \(d) toy_ref(n = 3, y = 15, seed = as.numeric(d)) + 10)
  names(set) <- format(dates)
  subjects <- data.frame(station = c("x", "y", "z"))
  res <- normals_table(set, dates, subjects, fam, "cy", min_valid = 10,
                       fallback = c(fitted = "none"), smoothed_mask = "none",
                       flavors = c("empirical", "fitted"), probs = c(0.1, 0.5, 0.9))
  expect_null(res$smooth)
  n <- res$normals
  expect_identical(nrow(n), 2L * 365L * 3L)
  expect_identical(names(n), c("flavor", "date", "day", "station", "distribution", "n_valid",
                               "ref_first_wy", "ref_last_wy", "mean", "mu", "sigma", "fit_status",
                               "q0100", "q0500", "q0900"))
  expect_identical(unique(n$distribution[n$flavor == "fitted"]), "normal")
  obs <- data.frame(station = c("y", "q"), date = dates[c(10, 10)], value = c(30, 1))
  pc <- percentiles_table(n, set, subjects, obs, fam, 10, c(fitted = "none"), 2,
                          flavors = c("empirical", "fitted"))
  expect_identical(pc$station, c("y", "y"))
  row <- n$flavor == "fitted" & n$date == dates[10] & n$station == "y"
  expect_identical(pc$ptile[2], as_ptile(stats::pnorm(30, n$mu[row], n$sigma[row]), 2))
  rm("normal", envir = normals:::.families)
})

test_that("the zero-inflated gamma and empirical families are registered on load", {
  expect_identical(get_family("zero_inflated_gamma"), family_zig())
  expect_identical(get_family("empirical"), family_empirical())
  expect_error(get_family("nope"), "no family")
  expect_identical(fit_status_label(0:5, family_zig()),
                   c("ok", "no_snow", "p_zero_only", "no_data", "p_zero_not_converged", "constant"))
})

test_that("the smoothed fit does not depend on cores", {
  skip_on_os("windows")
  zig <- family_zig()
  dates <- year_dates(2023, "wy")
  ref <- lapply(seq_along(dates), \(j) {
    s <- pmax(0, sin(pi * (j - 30) / 250))
    toy_ref(n = 7, y = 20, seed = j) * s
  })
  stats <- stack_stats(lapply(ref, \(h) day_stats(h, zig)))
  f1 <- smooth_fit(stats, cycle_day(dates, "wy"), zig, cores = 1L)
  f2 <- smooth_fit(stats, cycle_day(dates, "wy"), zig, cores = 2L)
  expect_identical(f1, f2)
  tabs <- smooth_tables(f1, data.frame(id = 1:7))
  expect_identical(names(tabs$coef)[c(1, 2, 26, 50, 73)], c("id", "nu_01", "l1_01", "lcv_01", "lcv_24"))
  expect_identical(names(tabs$diagnostics),
                   c("id", "status", "iter", "lambda_nu", "lambda_l1", "lambda_lcv", "edf_nu", "edf_l1", "edf_lcv"))
})
