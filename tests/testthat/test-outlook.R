# The outlook (R/outlook.R).
#
# 1. The 12 checks of mco-snotel outlook:scripts/test_outlook.R (cases 1-7),
#    on synthetic normals with known gamma parameters, tolerances unchanged.
# 2. Golden fixtures: mco-snotel's calc_station_outlook() and
#    outlook_lengths() at bb8d06a on synthetic stations
#    (data-raw/outlook-fixtures.R), reproduced exactly.
#
# Until the core lands, helper-outlook-shim.R stands in for its stubs.

probs  <- unname(NORMALS_QUANTILES)
q_cols <- names(NORMALS_QUANTILES)

wy_dates <- function(wy) seq(as.Date(sprintf("%d-10-01", wy - 1L)), as.Date(sprintf("%d-09-30", wy)), by = "day")

# Synthetic stations: 30 reference years, constant gamma normals every day
# (today's make_case()). Returns a function giving the outlook from an issue
# date with values x.
make_case <- function(wy, length_days, alpha = 1.3, shape = 4, scale = 5, p_zero = 0.01,
                      stations = c("1:XX:SNTL", "2:XX:SNTL")) {
  n <- length(stations)
  ref <- matrix(stats::rgamma(n * 30, shape, scale = scale), n, 30)
  dates <- wy_dates(wy)
  function(x, issue, beta = 0, which = stations) {
    i <- match(which, stations)
    issue <- as.Date(issue)
    tgt <- dates[dates > issue]
    di <- zig_dist(matrix(c(shape, scale, p_zero), length(i), 3, byrow = TRUE),
                   ref[i, , drop = FALSE], 10L)
    outlook(di, rep(list(di), length(tgt)), x, issue, tgt,
            data.frame(triplet = which), length_days, alpha, beta)
  }
}

test_that("1. rho ~ 0: quantiles = the normals' (to 0.1 in)", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2027L, length_days = 0.01)
  o <- k(2, "2026-12-01", which = "1:XX:SNTL")
  row <- o[o$date == as.Date("2027-01-15"), ]
  inner <- 2:18
  # p_zero 0.01 is below half a reference year (0.5 / 30), so a station with snow uses 0.
  ref <- stats::qgamma(probs[inner], 4, scale = 5)
  expect_lte(max(abs(unlist(row[q_cols[inner]]) - ref)), 0.051)
})

test_that("2. rho ~ 1: median = today's value (to 0.1 in)", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2027L, length_days = 1e6)
  o <- k(25, "2026-12-01", which = "1:XX:SNTL")
  expect_true(all(abs(o$q0500 - 25) <= 0.051))
})

test_that("3. quantiles never decrease with probability; p_zero in [0, 1]", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2027L, length_days = 60)
  o <- k(c(0, 30), "2027-01-10")
  m <- as.matrix(o[q_cols])
  expect_true(all(apply(m, 1, \(r) all(diff(r) >= -1e-9))))
  expect_true(all(o$p_zero >= 0 & o$p_zero <= 1))
})

test_that("4. more today never lowers the outlook; a station at zero starts at zero", {
  local_outlook_shim(); set.seed(1)
  swe <- c(0, 0.1, 1, 5, 20, 40)
  st <- paste0(seq_along(swe), ":XX:SNTL")
  k <- make_case(2027L, length_days = 60, p_zero = 0.001, stations = st)   # floored for snowy stations
  o <- k(swe, "2027-01-10")
  first <- o[o$date == as.Date("2027-01-11"), ]
  first <- first[match(st, first$triplet), ]
  expect_true(all(diff(first$q0500) >= -1e-9))
  expect_true(all(diff(first$p_zero) <= 1e-9))
  expect_identical(first$q0500[1], 0)
})

test_that("5. far above the normals: the range starts at today's value", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2027L, length_days = 60)
  o <- k(200, "2027-01-10", which = "1:XX:SNTL")
  expect_lt(abs(o$q0500[1] / 200 - 1), 0.1)
})

test_that("6. leap years: WY2028 includes Feb 29; issued Sep 30 gives no outlook", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2028L, length_days = 60)
  o <- k(10, "2028-02-01", which = "1:XX:SNTL")
  expect_true(as.Date("2028-02-29") %in% o$date)
  expect_null(k(10, "2028-09-30", which = "1:XX:SNTL"))
})

test_that("7. beta widens a low start and fades with rho", {
  local_outlook_shim(); set.seed(1)
  k <- make_case(2027L, length_days = 60)
  w <- \(b, d) { o <- k(3, "2027-01-10", beta = b, which = "1:XX:SNTL"); r <- o[o$date == d, ]; r$q0900 - r$q0100 }
  expect_gt(w(0.3, as.Date("2027-01-24")), w(0, as.Date("2027-01-24")))
  expect_lte(abs(w(0.3, as.Date("2027-07-28")) - w(0, as.Date("2027-07-28"))), 0.25)
})


# Golden fixtures -----------------------------------------------------------

golden <- function() readRDS(test_path("fixtures", "outlook-golden.rds"))

golden_dists <- function(cs, min_valid) {
  ref <- cs$ref_tenths / 10
  lapply(seq_along(cs$dates), \(t) {
    h <- ref[, , t]; colnames(h) <- cs$years
    zig_dist(cs$params[, , t], h, min_valid)
  })
}

alpha_30yr <- function() {
  tab <- utils::read.csv(system.file("extdata", "outlook_alpha.csv", package = "normals"),
                         comment.char = "#", stringsAsFactors = FALSE)
  tab[tab$period == "30yr" & tab$flavor == "zig_smooth", ]
}

test_that("outlook_lengths() reproduces mco-snotel's outlook_lengths() exactly", {
  local_outlook_shim()
  fx <- golden()
  for (cs in fx$cases) {
    dists <- golden_dists(cs, fx$min_valid)
    # Today's dates are sort(unique(normals$date)), stored as double.
    dates <- as.Date(as.numeric(cs$dates))
    got <- outlook_lengths(dists, dates, data.frame(triplet = fx$stations),
                           outlook_alpha_at(alpha_30yr(), dates, "wy"))
    expect_identical(got, cs$lengths[c("date", "triplet", "outlook_length", "outlook_alpha")])
  }
})

test_that("outlook() reproduces mco-snotel's calc_station_outlook() exactly", {
  local_outlook_shim()
  fx <- golden()
  dists <- lapply(fx$cases, golden_dists, min_valid = fx$min_valid)
  for (o in fx$outlooks) {
    cs <- fx$cases[[as.character(o$wy)]]
    d <- dists[[as.character(o$wy)]]
    t0 <- match(o$issue, cs$dates)
    tt <- seq.int(t0 + 1L, length(cs$dates))
    at <- cs$lengths[cs$lengths$date == o$issue, ]
    at <- at[match(fx$stations, at$triplet), ]
    got <- outlook(d[[t0]], d[tt], o$obs, o$issue, cs$dates[tt],
                   data.frame(triplet = fx$stations), at$outlook_length, at$outlook_alpha, o$beta)
    names(got)[names(got) == "p_zero"] <- "p_snow_free"   # the Phase 1 edge rename
    expect_identical(got, o$expected, label = paste("outlook from", o$issue, "beta", o$beta))
  }
})


# Pieces --------------------------------------------------------------------

test_that("the pieces follow their formulas", {
  expect_identical(outlook_persistence(30, 60, 1.3), exp(-(30 / 60)^1.3))
  expect_identical(outlook_spread(c(-2, 0, 1), 0.5, 0.25), c(exp(0.25), 1, 1))
  expect_identical(outlook_spread(-4, 1, 2), 3)
  expect_identical(outlook_crps(matrix(c(1, 2), 1), 1.5, c(0.25, 0.75)),
                   2 * mean(c((0 - 0.25) * (1 - 1.5), (1 - 0.75) * (2 - 1.5))))
  ctl <- outlook_control()
  expect_identical(ctl$k_max, 3)
  expect_identical(ctl$fit_lags, c(3, 7, 14, 21, 30, 45, 60, 90, 120))
  expect_error(outlook_control(k_max = "a"), "k_max")
})

test_that("outlook_alpha_at() interpolates on the cycle day, constant past the ends", {
  local_outlook_shim()
  tab <- data.frame(day = c(365, 1, 151, 152), alpha = c(2, 1, 1.5, 1.7))
  d <- as.Date(c("2027-10-01", "2028-02-28", "2028-02-29", "2028-03-01", "2028-09-30"))
  expect_identical(outlook_alpha_at(tab, d, "wy"),
                   stats::approx(c(1, 151, 152, 365), c(1, 1.5, 1.7, 2),
                                 xout = c(1, 151, 151.5, 152, 365), rule = 2)$y)
  expect_error(outlook_alpha_at(data.frame(x = 1), d, "wy"), "table")
})

test_that("outlook_pit() matches today's hindcast formula", {
  local_outlook_shim(); set.seed(2)
  n <- 4
  ref <- matrix(stats::rgamma(n * 30, 3, scale = 4), n, 30)
  ref[1, 1:5] <- 0
  par <- cbind(c(3, 3, NA, 3), c(4, 4, NA, 4), c(0.2, 0.01, 0.3, 1))
  d <- zig_dist(par, ref, 10L)
  x <- c(5, 0, 2, 0); y <- c(8, 0, 1, 3); u <- c(0.3, 0.6, 0.9, 0.1)
  s <- outlook_start(d, x)
  rho <- outlook_persistence(20, 40, 1.2)
  sg <- sqrt(pmax(1e-12, 1 - rho^2)) * outlook_spread(s$z0, rho, 0.25)
  dt <- normals:::outlook_floored(d, s$floor, outlook_control())
  ya <- y / s$anchor^rho
  p0 <- dist_p_zero(dt)
  want <- ifelse(ya > 0, stats::pnorm((stats::qnorm(dist_cdf(dt, ya)) - rho * s$z0) / sg),
                 u * ifelse(p0 > 0, stats::pnorm((stats::qnorm(p0) - rho * s$z0) / sg), 0))
  expect_identical(outlook_pit(d, s, rho, 0.25, y, u), want)
  expect_identical(s$floor, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("inputs are checked", {
  local_outlook_shim()
  d <- zig_dist(matrix(c(4, 5, 0.1), 1), matrix(1, 1, 30), 10L)
  expect_error(outlook_start(list(), 1), "normals_dist")
  expect_error(outlook_start(d, c(1, 2)), "`x`")
  expect_error(outlook(d, list(d), 1, as.Date("2027-01-01"), as.Date(c("2027-01-02", "2027-01-03")),
                       data.frame(triplet = "a"), 60, 1.3, 0), "targets")
  expect_error(outlook(d, list(d), 1, as.Date("2027-01-01"), as.Date("2027-01-02"),
                       data.frame(triplet = c("a", "b")), 60, 1.3, 0), "subjects")
  expect_error(outlook_lengths(list(d, d), as.Date(c("2027-01-01", "2027-01-03")),
                               data.frame(triplet = "a"), 1.3), "consecutive")
  expect_null(outlook(d, list(), 1, as.Date("2027-09-30"), as.Date(character()),
                      data.frame(triplet = "a"), 60, 1.3, 0))
})
