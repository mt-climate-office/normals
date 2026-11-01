test_that("year_of, year_start, year_dates and dates_through", {
  d <- as.Date(c("2023-09-30", "2023-10-01", "2024-02-29", "2024-12-31"))
  expect_identical(year_of(d, "wy"), c(2023L, 2024L, 2024L, 2025L))
  expect_identical(year_of(d, "cy"), c(2023L, 2023L, 2024L, 2024L))
  expect_identical(year_start(2024, "wy"), as.Date("2023-10-01"))
  expect_identical(year_start(2024, "cy"), as.Date("2024-01-01"))
  expect_length(year_dates(2024, "wy"), 366L)
  expect_length(year_dates(2027, "wy"), 365L)
  expect_equal(range(year_dates(2024, "cy")), as.Date(c("2024-01-01", "2024-12-31")))
  expect_equal(dates_through(as.Date("2023-10-03"), "wy"), as.Date("2023-10-01") + 0:2)
  expect_identical(day_key(as.Date("2024-02-29")), "02-29")
  expect_error(year_of(d, NULL), "year_type")
  expect_error(year_dates(2024), "year_type")
})

test_that("cycle_day: both year types, Feb 29 at the half day", {
  wy <- year_dates(2024, "wy")
  cd <- cycle_day(wy, "wy")
  expect_identical(cd[1], 1); expect_identical(cd[length(cd)], 365)
  expect_identical(cd[day_key(wy) == "02-29"], 151.5)
  expect_identical(cycle_day(year_dates(2027, "wy"), "wy"), as.numeric(1:365))
  cy <- year_dates(2024, "cy")
  expect_identical(cycle_day(cy, "cy")[day_key(cy) %in% c("02-28", "02-29", "03-01")], c(59, 59.5, 60))
  expect_identical(cycle_day(cy, "cy")[366], 365)
  # Today's normals_cycle_day() (mco-snodas R/utils.R), for water years.
  legacy <- function(date) {
    start <- as.Date(sprintf("%d-10-01", year_of(date, "wy") - 1L))
    pos <- as.numeric(date - start) + 1
    feb29 <- as.Date(sprintf("%d-02-29", year_of(date, "wy")), optional = TRUE)
    leap <- !is.na(feb29)
    pos[leap & date > feb29[leap]] <- pos[leap & date > feb29[leap]] - 1
    pos[leap & date == feb29[leap]] <- pos[leap & date == feb29[leap]] - 0.5
    pos
  }
  for (w in 2023:2025) expect_identical(cycle_day(year_dates(w, "wy"), "wy"), legacy(year_dates(w, "wy")))
})

test_that("reference_dates: Feb 29 -> Feb 28, and Oct-Dec dates by water year", {
  expect_identical(reference_dates(as.Date("2024-02-29"), 2021:2024, "wy"),
                   as.Date(c("2021-02-28", "2022-02-28", "2023-02-28", "2024-02-29")))
  expect_identical(reference_dates(as.Date("2023-11-15"), c(2004L, 2023L), "wy"),
                   as.Date(c("2003-11-15", "2022-11-15")))
  expect_identical(reference_dates(as.Date("2023-11-15"), c(2004L, 2023L), "cy"),
                   as.Date(c("2004-11-15", "2023-11-15")))
  # mco-snotel reference_dates() and mco-snodas normals_reference_years():
  # lubridate's %m-% years().
  skip_if_not_installed("lubridate")
  for (d in as.Date(c("2024-02-29", "2023-10-01", "2024-09-30", "2026-12-31"))) {
    d <- as.Date(d, origin = "1970-01-01")
    yrs <- 1990:2025
    expect_identical(reference_dates(d, yrs, "wy"),
                     lubridate::`%m-%`(d, lubridate::years(year_of(d, "wy") - yrs)))
  }
})

test_that("min_valid_years and period rules (spec section 5)", {
  expect_identical(min_valid_years(c(0, 3, 10, 30), cap = 10), c(1L, 3L, 10L, 10L))
  expect_error(min_valid_years(30), "cap")
  expect_identical(period_years("rolling", 2027, length = 30), 1997:2026)
  expect_identical(period_years("rolling", 2027, length = 30, earliest = 2004), 2004:2026)
  # WY2035: snodas = 30yr.
  expect_identical(period_years("rolling", 2035, length = 30, earliest = 2004),
                   period_years("rolling", 2035, length = 30))
  expect_identical(period_years("por", 2027, first_year = 1963), 1963:2026)
  noaa <- list(c(1991, 2020), c(2021, 2050))
  expect_null(period_years("fixed", 2021, periods = noaa, publish_lag = 1))
  expect_identical(period_years("fixed", 2022, periods = noaa, publish_lag = 1), 1991:2020)
  expect_identical(period_years("fixed", 2052, periods = noaa, publish_lag = 1), 2021:2050)
  # mco-snodas's NOAA role: absent before WY2042 with periods from 2011-2040.
  expect_null(period_years("fixed", 2041, periods = list(c(2011, 2040)), publish_lag = 1))
  expect_error(period_years("rolling", 2027, length = 30, first_year = 1963), "first_year")
  expect_identical(period_key(2004, 2026, "wy"), "wy2004-2026")
  expect_identical(period_key(1996, 2025, "cy"), "cy1996-2025")
  expect_identical(period_key(NULL, 2026, "wy", por = TRUE), "wypor-2026")
})

test_that("flavor and quantile vocabularies", {
  expect_identical(flavor_legacy(c("smoothed", "empirical"), "legacy"), c("zig_smooth", "ecdf"))
  expect_identical(flavor_legacy(c("zig", "ecdf"), "spec"), c("fitted", "empirical"))
  expect_error(flavor_legacy("zig", "legacy"), "unknown flavor")
  expect_error(flavor_legacy("fitted"), "to")
  expect_identical(quantile_columns(c(0, 0.025, 1)), c("q0000", "q0025", "q1000"))
  expect_error(quantile_columns(c(0.0101, 0.0102)), "share")
})
