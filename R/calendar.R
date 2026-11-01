#' Years, dates and cycle days
#'
#' Calendar helpers for water years (`"wy"`: Oct 1 - Sep 30, named for the
#' year it ends in) and calendar years (`"cy"`). `year_type` is never
#' defaulted.
#'
#' * `year_of()`: the year a date belongs to.
#' * `year_start()`: the first day of a year (Oct 1 of `year - 1`, or Jan 1).
#' * `year_dates()`: every day of a year, 365 or 366 dates.
#' * `dates_through()`: the start of `date`'s year through `date`.
#' * `cycle_day()`: position on the 365-day cycle used for smoothing. `"wy"`:
#'   Oct 1 = 1 ... Sep 30 = 365, Feb 29 = 151.5. `"cy"`: Jan 1 = 1 ...
#'   Dec 31 = 365, Feb 29 = 59.5.
#' * `day_key()`: `"MM-DD"`.
#'
#' @param date A `Date` vector.
#' @param year Integer year(s).
#' @param year_type `"wy"` or `"cy"`.
#' @return `year_of()`: integer. `year_start()`, `year_dates()`,
#'   `dates_through()`: `Date`. `cycle_day()`: double. `day_key()`: character.
#' @name calendar
NULL

#' @rdname calendar
#' @export
year_of <- function(date, year_type) not_implemented()

#' @rdname calendar
#' @export
year_start <- function(year, year_type) not_implemented()

#' @rdname calendar
#' @export
year_dates <- function(year, year_type) not_implemented()

#' @rdname calendar
#' @export
dates_through <- function(date, year_type) not_implemented()

#' @rdname calendar
#' @export
cycle_day <- function(date, year_type) not_implemented()

#' @rdname calendar
#' @export
day_key <- function(date) not_implemented()

#' Reference dates for one day
#'
#' The same calendar day as `date` in each of `years`, shifted by whole years
#' from `date`'s own year of `year_type`. Feb 29 becomes Feb 28 in a non-leap
#' year.
#'
#' @param date A single `Date`.
#' @param years Integer reference years (water or calendar years).
#' @param year_type `"wy"` or `"cy"`.
#' @return A `Date` vector of `length(years)`.
#' @export
reference_dates <- function(date, years, year_type) not_implemented()

#' Minimum valid years
#'
#' `max(1L, min(cap, n_years))`, elementwise: the rule both pipelines apply
#' today with `cap = 10`.
#'
#' @param n_years Integer, the number of reference years counted (see the
#'   design notes on what each pipeline counts).
#' @param cap Integer cap. No default.
#' @return Integer vector.
#' @export
min_valid_years <- function(n_years, cap) not_implemented()

#' Reference periods
#'
#' `period_years()` evaluates one role rule (spec section 5) for normals year
#' `year`:
#' * `"rolling"`: `length` years ending `year - 1`, starting no earlier than
#'   `earliest` if given;
#' * `"fixed"`: the latest of `periods` (a list of `c(first, last)`) whose
#'   `last + publish_lag < year`;
#' * `"por"`: `first_year` through `year - 1`.
#'
#' Arguments a rule does not use must be `NULL`.
#'
#' `period_key()` formats a period as in spec section 3: `"wy2004-2026"`,
#' `"cy1996-2025"`, or `"wypor-2026"` when `por = TRUE`.
#'
#' @param rule `"rolling"`, `"fixed"` or `"por"`.
#' @param year The normals year.
#' @param length,earliest Rolling rule parameters.
#' @param periods,publish_lag Fixed rule parameters.
#' @param first_year Period-of-record start.
#' @param first,last First and last years of a period.
#' @param year_type `"wy"` or `"cy"`.
#' @param por Format as a period of record.
#' @return `period_years()`: an integer vector, or `NULL` when the rule
#'   yields no period for `year`. `period_key()`: character.
#' @export
period_years <- function(rule, year, length = NULL, earliest = NULL,
                         periods = NULL, publish_lag = NULL,
                         first_year = NULL) not_implemented()

#' @rdname period_years
#' @export
period_key <- function(first, last, year_type, por = FALSE) not_implemented()
