#!/usr/bin/env Rscript
# Golden outlook fixtures: today's mco-snotel outlook (calc_station_outlook()
# and outlook_lengths(), mco-snotel `outlook` at bb8d06a) on synthetic
# stations, for tests/testthat/test-outlook.R to reproduce exactly.
#
# Run once, from an export of mco-snotel's outlook branch, never by the 2.0
# code:
#   git -C mco-snotel archive bb8d06a R src config | tar -x -C SNOTEL_DIR
#   Rscript data-raw/outlook-fixtures.R SNOTEL_DIR
# Writes tests/testthat/fixtures/outlook-golden.rds.
#
# Stations (zig_smooth normals given directly, seasonal by cycle day; the
# reference record is a persistent daily series drawn from those normals):
#   1  snowy; p_zero = 1 and no shape/scale when snow-free (status 1)
#   2  shallow; shape/scale NA where p_zero > 0.5 (ECDF fallback rows)
#   3  deep, p_zero 0.005 in winter (below half a year: the floor)
#   4  snowy, 8 reference years missing (n_valid 22)
#   5  no normals (every parameter NA)
#   6  shape/scale NA every day (all ECDF fallback), p_zero given
# Issue dates: WY2027 Dec 1, Feb 1, Apr 15 (beta 0.25, and Feb 1 with beta
# 0.6 too); WY2028 Feb 1 (leap year).

args <- commandArgs(TRUE)
if (length(args) < 1) stop("Usage: outlook-fixtures.R SNOTEL_DIR")
SNOTEL <- normalizePath(args[1])
OUT <- normalizePath(file.path("tests", "testthat", "fixtures"), mustWork = TRUE)

suppressPackageStartupMessages(library(magrittr))
owd <- setwd(SNOTEL)
for (f in c("utils", "zig_spline_normals", "normals", "outlook")) source(file.path("R", paste0(f, ".R")))
for (f in c("ecdf", "zig")) Rcpp::sourceCpp(file.path("src", paste0(f, ".cpp")))

stations <- sprintf("%d:MT:SNTL", 1:6)
n <- length(stations)

# Seasonal zig_smooth parameters for station s on cycle day cd.
env <- function(cd) pmax(0, sin(pi * (cd - 20) / 260))
params_for <- function(s, cd) {
  e  <- env(cd)
  m  <- c(1, 0.6, 1.5, 1, 1, 0.8)[s]
  p0 <- stats::plogis(6 * (0.3 - e))
  sh <- 2 + 6 * e
  sc <- m * (0.5 + 5 * e)
  if (s == 1) { bare <- e == 0; p0[bare] <- 1; sh[bare] <- NA; sc[bare] <- NA }
  if (s == 2) { fb <- p0 > 0.5; sh[fb] <- NA; sc[fb] <- NA }
  if (s == 3) p0 <- ifelse(e > 0.3, 0.005, p0)
  if (s == 5) { p0[] <- NA; sh[] <- NA; sc[] <- NA }
  if (s == 6) { sh[] <- NA; sc[] <- NA }
  cbind(shape = sh, scale = sc, p_zero = p0)
}
# Values on the 0.1 grid, made from integer tenths.
q_value <- function(u, p) {
  sh <- p[, "shape"]; sc <- p[, "scale"]; p0 <- p[, "p_zero"]
  sh[is.na(sh)] <- 3; sc[is.na(sc)] <- 2 * (1 - pmin(p0[is.na(sc)], 0.99)); p0[is.na(p0)] <- 0.3
  adj <- (u - p0) / (1 - p0)
  v <- ifelse(adj <= 0, 0, stats::qgamma(pmin(adj, 1 - 1e-9), sh, scale = sc))
  round(v * 10) / 10
}

set.seed(20261101)
# A persistent record: per station and water year, an AR(1) latent z (phi
# 0.97 a day) mapped through that day's normals.
first_wy <- 1997L; last_wy <- 2027L
days <- seq(as.Date(sprintf("%d-10-01", first_wy - 1L)), as.Date(sprintf("%d-09-30", last_wy)), by = "day")
rec <- vector("list", n)
for (s in seq_len(n)) {
  v <- numeric(length(days))
  for (wy in first_wy:last_wy) {
    i <- which(water_year(days) == wy)
    z <- numeric(length(i)); z[1] <- stats::rnorm(1)
    for (t in seq_along(i)[-1]) z[t] <- 0.97 * z[t - 1] + sqrt(1 - 0.97^2) * stats::rnorm(1)
    p <- params_for(if (s == 5) 1 else s, normals_cycle_day(days[i]))
    v[i] <- q_value(stats::pnorm(z), p)
  }
  if (s == 4) v[water_year(days) %in% c(1999, 2003, 2004, 2010, 2015, 2016, 2021, 2024)] <- NA
  rec[[s]] <- tibble::tibble(triplet = stations[s], date = days, wteq = as.integer(round(v * 10)) / 10)
}
rec <- dplyr::bind_rows(rec)

build <- function(wy) {
  years <- (wy - 30L):(wy - 1L)
  sm <- swe_matrix(dplyr::filter(rec, date < as.Date(sprintf("%d-10-01", wy - 1L))), stations)
  dates <- normals_dates(wy)
  ref <- lapply(dates, \(d) reference_record(sm, d, years))
  par <- lapply(dates, \(d) do.call(rbind, lapply(seq_len(n), \(s) params_for(s, normals_cycle_day(d)))))
  normals <- dplyr::bind_rows(lapply(seq_along(dates), \(t) tibble::tibble(
    flavor = OUTLOOK_FLAVOR, date = dates[t], triplet = stations,
    n_valid = as.integer(rowSums(!is.na(ref[[t]]))),
    zig_shape = par[[t]][, 1], zig_scale = par[[t]][, 2], zig_p_zero = par[[t]][, 3])))
  lengths <- outlook_lengths(sm, normals, years, "30yr")
  normals <- dplyr::left_join(normals, lengths, by = c("flavor", "date", "triplet"))
  list(wy = wy, years = years, sm = sm, dates = dates, ref = ref, params = par,
       normals = normals, lengths = lengths)
}

obs_for <- list(
  "2026-12-01" = c(1.2, 0, 0.1, 30, 2, 0.4),
  "2027-02-01" = c(4.5, 0, 0.3, 0.2, 5, 1.1),
  "2027-04-15" = c(9, 2, 50, NA, 4, 0),
  "2028-02-01" = c(6, 0.5, 15, 3, 2, 0)
)
issues <- list(
  list(wy = 2027L, issue = "2026-12-01", beta = 0.25),
  list(wy = 2027L, issue = "2027-02-01", beta = 0.25),
  list(wy = 2027L, issue = "2027-02-01", beta = 0.6),
  list(wy = 2027L, issue = "2027-04-15", beta = 0.25),
  list(wy = 2028L, issue = "2028-02-01", beta = 0.25)
)

cases <- list(); outlooks <- list()
for (wy in c(2027L, 2028L)) {
  b <- build(wy)
  cases[[as.character(wy)]] <- list(
    wy = wy, years = b$years, dates = b$dates,
    # Reference values as integer tenths (exact on reload), n x years x dates.
    ref_tenths = array(as.integer(round(unlist(b$ref) * 10)), c(n, length(b$years), length(b$dates))),
    params = array(unlist(b$params), c(n, 3L, length(b$dates))),
    lengths = b$lengths)
  stopifnot(identical(as.vector(cases[[as.character(wy)]]$ref_tenths) / 10, unlist(b$ref)))
  for (k in seq_along(issues)) {
    is <- issues[[k]]
    if (is$wy != wy) next
    obs <- tibble::tibble(triplet = stations, wteq = obs_for[[is$issue]])
    o <- calc_station_outlook(b$sm, b$normals, b$years, obs, as.Date(is$issue), beta = is$beta)
    outlooks[[length(outlooks) + 1L]] <- list(wy = wy, issue = as.Date(is$issue), beta = is$beta,
                                              obs = obs$wteq, expected = o)
  }
}
setwd(owd)
fx <- list(stations = stations, min_valid = 10L, cases = cases, outlooks = outlooks,
           source = list(repo = "mt-climate-office/mco-snotel", ref = "outlook",
                         commit = "bb8d06a5703af8db21b8dc7d2f4ab3e20ef44917",
                         files = c("R/outlook.R", "R/normals.R", "R/utils.R", "src/zig.cpp",
                                   "src/ecdf.cpp", "src/ecdf_impl.h", "config/outlook_alpha.csv"),
                         r = R.version.string, made = format(Sys.Date())))
saveRDS(fx, file.path(OUT, "outlook-golden.rds"), version = 3, compress = "xz")
message("Wrote ", file.path(OUT, "outlook-golden.rds"))
