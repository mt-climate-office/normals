#!/usr/bin/env Rscript
# Calibrate the outlook from data-raw/outlook-hindcast.R output. Ported from
# mco-snotel outlook:scripts/outlook-calibrate.R (bb8d06a); the method is
# unchanged, only the paths differ. Not run in CI.
#
#   alpha  (hindcast mode alpha) the persistence shape: writes
#          OUT_CSV (default inst/extdata/outlook_alpha.csv).
#   beta   (hindcast mode beta) the spread's dependence on the starting z:
#          writes OUT_CSV (default inst/extdata/outlook_spread.csv). beta is
#          the grid value whose 80% and 50% ranges come closest to their
#          targets across starting classes (bare, D4, D3-D2, D1-D0, normal,
#          W0-W1, W2-W3, W4), each class weighted equally, summed over leads.
#
# For each period (and each flavor in the hindcast; the pipeline uses only
# zig_smooth), alpha follows a cyclic Fourier series over the
# water year (two yearly harmonics), fitted to the hindcast's coverage error:
# at each issue date and alpha, the squared error of the central 50% and 80%
# ranges (from pit, every case, zero or not) summed over leads. The curve
# minimises that error summed over the issue dates, interpolating each date's
# error across the alpha grid, so dates where alpha barely matters (late
# summer, when nearly every station is bare) pull little on it.
#
# Prints, leaving each test year out, the coverage and CRPS skill of the curve
# against picking alpha month by month.
#
# The two interact a little: calibrate beta, then alpha with that beta, and
# check the classes again with a final beta-mode hindcast.
#
# The tables are a dataset's calibration: the pipeline keeps its own copy (for
# SNOTEL, mco-snotel config/) and passes alpha and beta to the outlook.
#
# Usage (from the package root):
#   Rscript data-raw/outlook-calibrate.R HINDCAST_DIR [alpha|beta] [OUT_CSV]
#     HINDCAST_DIR  holds hindcast_{period}.parquet (alpha) or
#                   hindcast_beta_{period}.parquet (beta) for 30yr, snodas, por

suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(TRUE)
D    <- args[1]
WHAT <- if (length(args) >= 2) args[2] else "alpha"
stopifnot(WHAT %in% c("alpha", "beta"))
OUT  <- if (length(args) >= 3) args[3] else
  file.path("inst", "extdata", if (WHAT == "alpha") "outlook_alpha.csv" else "outlook_spread.csv")
K   <- 2L   # yearly harmonics
PERIODS <- c("30yr", "snodas", "por")

issue_day <- function(issue) {   # cycle day (Oct 1 = 1) of "MM-DD" in a non-leap year
  d <- as.Date(paste0(ifelse(substr(issue, 1, 2) >= "10", "2021-", "2022-"), issue))
  as.integer(d - as.Date("2021-10-01")) + 1L
}
fourier <- function(day, coef) {
  w <- 2 * pi * (day - 1) / 365; a <- coef[1]
  for (k in seq_len((length(coef) - 1) / 2)) a <- a + coef[2 * k] * cos(k * w) + coef[2 * k + 1] * sin(k * w)
  pmin(2.4, pmax(0.3, a))
}
# Coverage error per issue date and alpha, from counts.
err_surface <- function(cnt) cnt %>% group_by(issue, alpha, lead) %>% summarise(across(c(n, k80, k50), sum), .groups = "drop") %>%
  mutate(e = (k80 / n - 0.8)^2 + (k50 / n - 0.5)^2) %>% group_by(issue, alpha) %>% summarise(e = sum(e), .groups = "drop")
fit_curve <- function(es) {
  E <- tidyr::pivot_wider(es, id_cols = issue, names_from = alpha, values_from = e) %>% arrange(issue_day(issue))
  days <- issue_day(E$issue); M <- as.matrix(E[, -1]); grid <- as.numeric(colnames(M))
  obj <- function(coef) { a <- fourier(days, coef)
    sum(vapply(seq_along(days), \(i) stats::approx(grid, M[i, ], xout = a[i], rule = 2)$y, 0)) }
  starts <- list(c(1.2, rep(0, 2 * K)), c(1.3, -0.5, 0, rep(0, 2 * K - 2)))
  fits <- lapply(starts, \(s) stats::optim(s, obj, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10)))
  fits[[which.min(vapply(fits, `[[`, 0, "value"))]]$par
}

# Starting class of each case, from the issue date's SWE and z.
CLASSES <- c("Bare", "D4", "D3-D2", "D1-D0", "Normal", "W0-W1", "W2-W3", "W4")
start_class <- function(s0, z0) {
  pt <- 100 * stats::pnorm(z0)
  factor(dplyr::case_when(s0 == 0 ~ "Bare", pt < 2 ~ "D4", pt < 10 ~ "D3-D2", pt < 30 ~ "D1-D0", pt <= 70 ~ "Normal",
                          pt <= 90 ~ "W0-W1", pt <= 98 ~ "W2-W3", TRUE ~ "W4"), levels = CLASSES)
}

if (WHAT == "beta") {
  out <- list(); cls <- list()
  for (period in PERIODS) {
    s <- arrow::read_parquet(file.path(D, paste0("hindcast_beta_", period, ".parquet")),
                             col_select = c("method", "wy", "lead", "s0", "z0", "pit")) %>%
      filter(method != "climatology") %>%
      mutate(beta = as.numeric(sub("^b", "", method)), start = start_class(s0, z0),
             in80 = abs(pit - 0.5) <= 0.4, in50 = abs(pit - 0.5) <= 0.25)
    cov <- s %>% group_by(beta, start, lead) %>% summarise(c80 = mean(in80), c50 = mean(in50), .groups = "drop")
    err <- cov %>% group_by(beta) %>% summarise(e = sum((c80 - 0.8)^2 + (c50 - 0.5)^2), .groups = "drop")
    best <- err$beta[which.min(err$e)]
    out[[period]] <- tibble(period = period, beta = best)
    cls[[period]] <- s %>% filter(beta %in% c(0, best)) %>% group_by(period = period, beta, start) %>%
      summarise(n = n(), cov80 = mean(in80), cov50 = mean(in50), .groups = "drop")
  }
  cat("80% and 50% range coverage by starting class, beta = 0 and the chosen beta\n")
  print(as.data.frame(bind_rows(cls) %>% mutate(across(c(cov80, cov50), \(v) round(v, 2))) %>%
    tidyr::pivot_wider(names_from = beta, values_from = c(cov80, cov50), names_prefix = "b")))
  x <- bind_rows(out)
  con <- file(OUT, "w")
  writeLines(c(
    "# Outlook spread dependence on the starting z, by reference period:",
    "# the spread is multiplied by k = exp(-beta * rho * min(z0, 0)), at most 3, so",
    "# stations starting below normal get a wider range while persistence holds.",
    paste0("# Calibrated from the hindcast (normals data-raw/outlook-hindcast.R, data-raw/outlook-calibrate.R; ",
           format(Sys.Date()), ").")), con)
  utils::write.csv(x, con, row.names = FALSE, quote = FALSE)
  close(con)
  print(as.data.frame(x))
  message("Wrote ", OUT)
  quit(save = "no")
}

table <- list(); report <- list()
for (period in PERIODS) {
  s <- arrow::read_parquet(file.path(D, paste0("hindcast_", period, ".parquet")),
                           col_select = c("method", "flavor", "triplet", "wy", "issue", "lead", "crps", "pit")) %>%
    mutate(in80 = abs(pit - 0.5) <= 0.4, in50 = abs(pit - 0.5) <= 0.25)
  clim <- s %>% filter(method == "climatology") %>% select(flavor, triplet, wy, issue, lead, crps_clim = crps)
  s <- s %>% filter(method != "climatology") %>% mutate(alpha = as.numeric(sub("^a", "", method)))
  cnt <- s %>% group_by(flavor, wy, issue, alpha, lead) %>% summarise(n = n(), k80 = sum(in80), k50 = sum(in50), .groups = "drop")
  for (fl in sort(unique(s$flavor))) {
    cf <- filter(cnt, flavor == fl)
    table[[paste(period, fl)]] <- tibble(period = period, flavor = fl, day = 1:365,
                                         alpha = round(fourier(1:365, fit_curve(err_surface(cf))), 3))
    # Leave one year out: the curve and month-by-month picks, each fitted on the other years.
    picks <- purrr::map_dfr(sort(unique(cf$wy)), \(Y) {
      es <- err_surface(filter(cf, wy != Y))
      monthly <- es %>% mutate(m = substr(issue, 1, 2)) %>% group_by(m, alpha) %>% summarise(e = sum(e), .groups = "drop") %>%
        group_by(m) %>% slice_min(e, n = 1, with_ties = FALSE) %>% select(m, monthly = alpha)
      tibble(wy = Y, issue = unique(es$issue)) %>%
        mutate(curve = round(fourier(issue_day(issue), fit_curve(es)), 1), m = substr(issue, 1, 2)) %>%
        left_join(monthly, by = "m") %>% select(-m)
    })
    report[[paste(period, fl)]] <- purrr::map_dfr(c("curve", "monthly"), \(v)
      s %>% filter(flavor == fl) %>% inner_join(transmute(picks, wy, issue, alpha = .data[[v]]), by = c("wy", "issue", "alpha")) %>%
        inner_join(clim, by = c("flavor", "triplet", "wy", "issue", "lead")) %>%
        summarise(period = period, flavor = fl, version = v, cov80 = mean(in80), cov50 = mean(in50),
                  skill = 1 - sum(crps) / sum(crps_clim)))
  }
  rm(s, cnt); gc(verbose = FALSE)
}

cat("Leave one year out: coverage of the 80% and 50% ranges, CRPS skill against the normals alone\n")
print(as.data.frame(bind_rows(report) %>% mutate(across(c(cov80, cov50, skill), \(v) round(v, 3)))))

x <- bind_rows(table) %>% arrange(match(period, PERIODS), flavor, day)
con <- file(OUT, "w")
writeLines(c(
  "# Outlook persistence shape alpha for the smoothed normals (zig_smooth), by",
  "# reference period and day of the 365-day cycle (Oct 1 = 1; Feb 29 is read",
  "# between days 151 and 152).",
  "# rho(h) = exp(-(h / l)^alpha); l is fitted per station in the annual build.",
  "# Calibrated, not calculated: a cyclic Fourier series (2 harmonics) fitted, on",
  "# every test year, to the coverage error of the 50% and 80% outlook ranges in an",
  "# out-of-sample hindcast (checked leaving each year out). Regenerate with",
  "# normals data-raw/outlook-hindcast.R and",
  paste0("# data-raw/outlook-calibrate.R (", format(Sys.Date()), ").")), con)
utils::write.csv(x, con, row.names = FALSE, quote = FALSE)
close(con)
message("Wrote ", OUT, " (", nrow(x), " rows)")
