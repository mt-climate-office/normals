#!/usr/bin/env Rscript
# Hindcast of the outlook (R/outlook.R), for calibrating its persistence shape
# alpha and its spread factor beta (data-raw/outlook-calibrate.R). Ported from
# mco-snotel outlook:scripts/outlook-hindcast.R (bb8d06a) onto the normals 2.0
# API; the method is unchanged. Not run in CI.
#
# Out of sample: for each test year Y, the smoothed normals (one period) and
# the persistence lengths come only from the record before Y. Outlooks are
# issued on the 1st and 15th of every month and scored 7-120 days on, plus the
# normals alone ("climatology") as the baseline. Two modes:
#   alpha  every alpha on a grid (aX.X), beta from SPREAD_CSV
#   beta   alpha from ALPHA_CSV, every beta on a grid (bX.XX)
#
# Each row: method, flavor, triplet, wy, issue (MM-DD), lead, obs, s0 and z0
# (the value and z on the issue date), crps (outlook_crps(), from the
# 5th-95th percentiles in steps of 5) and pit (outlook_pit(); a zero outcome
# is placed at random within the forecast's chance of zero).
#
# The record is SNOTEL's (mco-snotel archive/snotel.parquet: triplet, date,
# wteq); periods are SNOTEL's (30yr, snodas, por). Another dataset adapts the
# reading and the period rules at the top.
#
# Usage (from the package root, with normals installed):
#   Rscript data-raw/outlook-hindcast.R ARCHIVE STATIONS OUT_DIR PERIOD [N_STATIONS] [CORES] [MODE] [ALPHA_CSV] [SPREAD_CSV]
#     ARCHIVE    local copy of archive/snotel.parquet
#     STATIONS   local copy of stations/stations.parquet
#     PERIOD     30yr, snodas or por
#     MODE       alpha (default) or beta
#     ALPHA_CSV  default inst/extdata/outlook_alpha.csv
#     SPREAD_CSV default inst/extdata/outlook_spread.csv
# Writes OUT_DIR/hindcast_{PERIOD}.parquet (alpha) or hindcast_beta_{PERIOD}.parquet.

suppressPackageStartupMessages({ library(normals); library(magrittr) })

args <- commandArgs(TRUE)
if (length(args) < 4) stop("Usage: outlook-hindcast.R ARCHIVE STATIONS OUT_DIR PERIOD [N_STATIONS] [CORES] [MODE] [ALPHA_CSV] [SPREAD_CSV]")
ARCHIVE <- args[1]; STATIONS <- args[2]; OUT <- args[3]; PERIOD <- args[4]
N_ST  <- if (length(args) >= 5) as.integer(args[5]) else 120L
CORES <- if (length(args) >= 6) as.integer(args[6]) else 8L
MODE  <- if (length(args) >= 7) args[7] else "alpha"
ALPHA_CSV  <- if (length(args) >= 8) args[8] else file.path("inst", "extdata", "outlook_alpha.csv")
SPREAD_CSV <- if (length(args) >= 9) args[9] else file.path("inst", "extdata", "outlook_spread.csv")
stopifnot(PERIOD %in% c("30yr", "snodas", "por"), MODE %in% c("alpha", "beta"))

ALPHAS   <- seq(0.3, 2.4, by = 0.1)
BETAS    <- seq(0, 0.8, by = 0.05)
LEADS    <- c(7, 14, 30, 60, 90, 120)
ISSUE_MD <- sprintf("%02d-%02d", rep(c(10:12, 1:9), each = 2), rep(c(1, 15), 12))
PROBS    <- seq(0.05, 0.95, by = 0.05)
FLAVOR   <- "zig_smooth"     # the published name of the smoothed flavor
CONTROL  <- outlook_control()
FAMILY   <- family_zig()
set.seed(20261031)

read_tab <- function(path) utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE)
alpha_tab <- read_tab(ALPHA_CSV)
alpha_tab <- alpha_tab[alpha_tab$period == PERIOD & alpha_tab$flavor == FLAVOR, ]
spread_beta <- local({
  b <- if (file.exists(SPREAD_CSV)) read_tab(SPREAD_CSV) else data.frame(period = character(), beta = numeric())
  b <- b$beta[b$period == PERIOD]
  if (length(b)) b[1] else 0
})

# SNOTEL's period rules (mco-snotel R/normals.R period_years()).
period_of <- function(wy, first_wy) {
  y <- switch(PERIOD,
              `30yr` = period_years("rolling", wy, length = 30L),
              snodas = period_years("rolling", wy, length = 30L, earliest = 2004L),
              por    = period_years("por", wy, first_year = first_wy))
  y[y >= first_wy]
}

# Active stations with records from WY1987 or earlier, sampled across states.
st <- arrow::read_parquet(STATIONS) %>% dplyr::filter(active, begin_date <= as.Date("1986-10-01"))
pick <- st %>% dplyr::group_by(state) %>%
  dplyr::slice_sample(n = ceiling(N_ST / dplyr::n_distinct(st$state))) %>% dplyr::ungroup()
stations <- sort(head(pick$triplet[sample(nrow(pick))], N_ST))
rec <- arrow::open_dataset(ARCHIVE) %>% dplyr::select(triplet, date, wteq) %>%
  dplyr::filter(!is.na(wteq), triplet %in% stations) %>% dplyr::collect()
last <- max(rec$date)
last_wy <- year_of(last, "wy") - (format(last, "%m-%d") != "09-30")
TEST_YEARS <- (last_wy - 10L):last_wy
message(length(stations), " stations, ", PERIOD, ", test years WY", min(TEST_YEARS), "-", max(TEST_YEARS))

run_year <- function(Y) {
  t0 <- Sys.time()
  rm <- record_matrix(dplyr::filter(rec, date < year_start(Y, "wy")), "triplet", "date", "wteq",
                      subjects = stations)
  years <- period_of(Y, year_of(rm$origin, "wy"))
  mv <- min_valid_years(length(years), cap = 10)
  dates <- year_dates(Y, "wy"); L <- length(dates)
  refs <- lapply(dates, \(d) reference_record(rm, d, years, "wy"))
  names(refs) <- format(dates)
  subjects <- data.frame(triplet = stations)
  nt <- normals_table(refs, dates, subjects, FAMILY, "wy", mv,
                      fallback = c(fitted = "empirical", smoothed = "empirical"),
                      smoothed_mask = "min_valid", flavors = "smoothed")$normals
  dists <- lapply(seq_len(L), \(t) {
    x <- nt[nt$date == dates[t], ]
    x <- x[match(stations, x$triplet), ]
    normals_dist(FAMILY, as.matrix(x[FAMILY$params]), refs[[t]], mv, "empirical")
  })
  # The lag correlations once, the lengths for each grid point.
  cr <- normals:::outlook_correlations(dists, CONTROL)
  fit_len <- \(a) normals:::outlook_fit_lengths(cr, a, CONTROL)
  grid <- if (MODE == "alpha")
    lapply(ALPHAS, \(a) list(label = sprintf("a%.1f", a), alpha = rep(a, L), len = fit_len(a),
                             beta = spread_beta))
  else {
    a <- outlook_alpha_at(alpha_tab, dates, "wy"); len <- fit_len(a)
    lapply(BETAS, \(b) list(label = sprintf("b%.2f", b), alpha = a, len = len, beta = b))
  }
  act <- dplyr::filter(rec, date >= dates[1], date <= dates[L])
  y_of <- function(d) act$wteq[act$date == d][match(stations, act$triplet[act$date == d])]

  out <- list()
  for (md in ISSUE_MD) {
    d0 <- as.Date(sprintf("%d-%s", if (md >= "10") Y - 1L else Y, md)); i0 <- match(d0, dates)
    if (is.na(i0)) next
    s0 <- y_of(d0)
    for (h in LEADS) {
      if (i0 + h > L) next
      y <- y_of(dates[i0 + h])
      ok <- !is.na(s0) & !is.na(y)
      if (!any(ok)) next
      start <- outlook_start(dists[[i0]], s0, CONTROL)
      tgt <- dists[[i0 + h]]
      # The target as the outlook reads it (the zero-chance floor where s0 > 0).
      tf <- normals:::outlook_floored(tgt, start$floor, CONTROL)
      p0 <- dist_p_zero(tf)
      cdf_y <- dist_cdf(tf, y)
      good <- ok & !is.na(start$z0) & !is.na(cdf_y) & !is.na(p0)
      if (!any(good)) next
      add <- function(method, q, pit) out[[length(out) + 1L]] <<- tibble::tibble(
        method = method, flavor = FLAVOR, triplet = stations[good], wy = Y, issue = md, lead = h,
        obs = y[good], s0 = s0[good], z0 = start$z0[good],
        crps = outlook_crps(q[good, , drop = FALSE], y[good], PROBS), pit = pit[good])
      rnd <- stats::runif(length(y))
      add("climatology",
          normals:::outlook_q(tf, matrix(PROBS, length(y), length(PROBS), byrow = TRUE)),
          ifelse(y > 0, cdf_y, rnd * p0))
      for (g in grid) {
        rho <- outlook_persistence(h, g$len[, i0], g$alpha[i0])
        add(g$label, outlook_quantiles(tgt, start, rho, g$beta, PROBS, CONTROL),
            outlook_pit(tgt, start, rho, g$beta, y, rnd, CONTROL))
      }
    }
  }
  message(PERIOD, " WY", Y, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  dplyr::bind_rows(out)
}

scores <- parallel::mclapply(TEST_YEARS, run_year, mc.cores = CORES, mc.preschedule = FALSE)
bad <- vapply(scores, \(x) inherits(x, "try-error"), TRUE)
if (any(bad)) stop("WY", TEST_YEARS[which(bad)[1]], " failed: ", conditionMessage(attr(scores[[which(bad)[1]]], "condition")))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
path <- file.path(OUT, paste0(if (MODE == "beta") "hindcast_beta_" else "hindcast_", PERIOD, ".parquet"))
arrow::write_parquet(dplyr::bind_rows(scores), path, compression = "zstd")
message("Wrote ", path)
