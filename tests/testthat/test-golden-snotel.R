# mco-snotel calc_station_normals() / calc_station_percentiles() (API.md 6.3),
# reproduced with the 2.0 calls and the Phase 1 renames (6.4).

snotel_legacy_normals <- function(x) {
  q <- grep("^q[0-9]{4}$", names(x), value = TRUE)
  tibble::as_tibble(c(
    list(flavor = flavor_legacy(x$flavor, "legacy"), date = x$date, triplet = x$triplet,
         n_valid = x$n_valid, ref_first_wy = x$ref_first_wy, ref_last_wy = x$ref_last_wy,
         mean = x$mean, zig_shape = x$shape, zig_scale = x$scale, zig_p_zero = x$p_zero),
    as.list(x[q])))
}

snotel_legacy_percentiles <- function(x) {
  tibble::tibble(flavor = flavor_legacy(x$flavor, "legacy"), date = x$date, triplet = x$triplet,
                 swe = x$value, swe_ptile = x$ptile, swe_standardized = x$standardized)
}

# Every period's 2.0 tables, computed once for both tests below.
snotel_results <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    fx  <- golden("snotel")
    rec <- record_matrix(fx$record, "triplet", "date", "wteq")
    # Stations in C-locale order, as today's dplyr::arrange(triplet).
    expect_identical(rec$subjects, sort(fx$stations, method = "radix"))
    subjects <- tibble::tibble(triplet = rec$subjects)
    zig <- family_zig()
    first_wy <- year_of(rec$origin, "wy")
    rules <- list(
      `30yr` = period_years("rolling", fx$wy, length = 30),
      snodas = period_years("rolling", fx$wy, length = 30, earliest = 2004),
      por    = period_years("por", fx$wy, first_year = first_wy))
    out <- lapply(stats::setNames(names(fx$expected), names(fx$expected)), \(period) {
      ex <- fx$expected[[period]]
      years <- rules[[period]]
      years <- years[years >= first_wy]        # mco-snotel normals_years() (P2)
      ref <- \(d) reference_record(rec, d, years, "wy")
      mv  <- min_valid_years(length(years), cap = 10)
      fb  <- c(fitted = "empirical", smoothed = "empirical")
      res <- normals_table(ref = ref, dates = year_dates(fx$wy, "wy"), subjects = subjects,
                           family = zig, year_type = "wy", min_valid = mv, fallback = fb,
                           smoothed_mask = "min_valid")
      obs <- fx$obs
      names(obs)[names(obs) == "wteq"] <- "value"
      pc <- percentiles_table(res$normals, ref = ref, subjects = subjects, obs = obs, family = zig,
                              min_valid = mv, fallback = fb, ptile_digits = 2)
      list(years = years, res = res, subjects = subjects,
           normals = snotel_legacy_normals(res$normals[res$normals$date %in% fx$keep_dates, ]),
           percentiles = snotel_legacy_percentiles(pc))
    })
    cache <<- list(fx = fx, out = out)
    cache
  }
})

test_that("SNOTEL ECDF and ZIG normals and their percentiles match mco-snotel, every period", {
  r <- snotel_results()
  for (period in names(r$out)) {
    ex <- r$fx$expected[[period]]; got <- r$out[[period]]
    expect_identical(got$years, ex$years)
    for (fl in c("ecdf", "zig")) {
      expect_golden(got$normals[got$normals$flavor == fl, ], ex$normals[ex$normals$flavor == fl, ],
                    label = paste(period, fl, "normals"))
      expect_golden(got$percentiles[got$percentiles$flavor == fl, ],
                    ex$percentiles[ex$percentiles$flavor == fl, ], label = paste(period, fl, "percentiles"))
    }
  }
})

test_that("SNOTEL ZIG-smooth normals, fit and percentiles match mco-snotel, every period", {
  skip_smoothed_off_platform()
  r <- snotel_results()
  for (period in names(r$out)) {
    ex <- r$fx$expected[[period]]; got <- r$out[[period]]
    expect_golden(got$normals, ex$normals, label = paste(period, "normals (all flavors)"))
    expect_golden(got$percentiles, ex$percentiles, label = paste(period, "percentiles (all flavors)"))
    tabs <- smooth_tables(got$res$smooth, got$subjects)
    expect_golden(tabs$coef, ex$coef, label = paste(period, "coef"))
    expect_golden(tabs$diagnostics, ex$diagnostics, label = paste(period, "diagnostics"))
  }
})
