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

test_that("SNOTEL normals, smoothed fit and percentiles match mco-snotel, every period", {
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
  for (period in names(fx$expected)) {
    ex <- fx$expected[[period]]
    years <- rules[[period]]
    years <- years[years >= first_wy]        # mco-snotel normals_years() (P2)
    expect_identical(years, ex$years)
    ref <- \(d) reference_record(rec, d, years, "wy")
    mv  <- min_valid_years(length(years), cap = 10)
    res <- normals_table(ref = ref, dates = year_dates(fx$wy, "wy"), subjects = subjects,
                         family = zig, year_type = "wy", min_valid = mv,
                         fallback = c(fitted = "empirical", smoothed = "empirical"),
                         smoothed_mask = "min_valid")
    got <- snotel_legacy_normals(res$normals[res$normals$date %in% fx$keep_dates, ])
    expect_golden(got, ex$normals, label = paste(period, "normals"))
    # Flavor by flavor too, so a platform difference is located.
    for (fl in unique(ex$normals$flavor))
      expect_golden(got[got$flavor == fl, ], ex$normals[ex$normals$flavor == fl, ],
                    label = paste(period, fl, "normals"))

    tabs <- smooth_tables(res$smooth, subjects)
    expect_golden(tabs$coef, ex$coef, label = paste(period, "coef"))
    expect_golden(tabs$diagnostics, ex$diagnostics, label = paste(period, "diagnostics"))

    obs <- fx$obs
    names(obs)[names(obs) == "wteq"] <- "value"
    pc <- percentiles_table(res$normals, ref = ref, subjects = subjects, obs = obs, family = zig,
                            min_valid = mv, fallback = c(fitted = "empirical", smoothed = "empirical"),
                            ptile_digits = 2)
    expect_golden(snotel_legacy_percentiles(pc), ex$percentiles, label = paste(period, "percentiles"))
  }
})
