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
year_of <- function(date, year_type) {
  check_date(date, "date")
  check_year_type(year_type)
  lt <- as.POSIXlt(date)
  y <- lt$year + 1900L
  if (year_type == "wy") y <- y + as.integer(lt$mon >= 9L)
  as.integer(y)
}

#' @rdname calendar
#' @export
year_start <- function(year, year_type) {
  check_whole(year, "year")
  check_year_type(year_type)
  year <- as.integer(year)
  if (year_type == "wy") make_date(year - 1L, 10L, 1L) else make_date(year, 1L, 1L)
}

#' @rdname calendar
#' @export
year_dates <- function(year, year_type) {
  check_whole(year, "year", single = TRUE)
  check_year_type(year_type)
  seq(year_start(year, year_type), year_start(year + 1L, year_type) - 1L, by = "day")
}

#' @rdname calendar
#' @export
dates_through <- function(date, year_type) {
  check_date(date, "date", single = TRUE)
  check_year_type(year_type)
  seq(year_start(year_of(date, year_type), year_type), date, by = "day")
}

#' @rdname calendar
#' @export
cycle_day <- function(date, year_type) {
  check_date(date, "date")
  check_year_type(year_type)
  # As mco-snodas normals_cycle_day() (wy): the day of the year, less one
  # after Feb 29 and less a half on it.
  yr    <- year_of(date, year_type)
  pos   <- as.numeric(date - year_start(yr, year_type)) + 1
  feb29 <- make_date(if (year_type == "wy") yr else yr, 2L, 29L)   # NA in non-leap years
  leap  <- !is.na(feb29)
  pos[leap & date >  feb29[leap]] <- pos[leap & date > feb29[leap]] - 1
  pos[leap & date == feb29[leap]] <- pos[leap & date == feb29[leap]] - 0.5
  pos
}

#' @rdname calendar
#' @export
day_key <- function(date) {
  check_date(date, "date")
  format(date, "%m-%d")
}

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
reference_dates <- function(date, years, year_type) {
  check_date(date, "date", single = TRUE)
  check_whole(years, "years")
  check_year_type(year_type)
  # date shifted back by whole years, as lubridate's %m-% years(): Feb 29
  # rolls back to Feb 28 where the target year has none.
  lt  <- as.POSIXlt(date)
  cal <- lt$year + 1900L - (year_of(date, year_type) - as.integer(years))
  out <- make_date(cal, lt$mon + 1L, lt$mday)
  bad <- is.na(out)
  out[bad] <- make_date(cal[bad], lt$mon + 1L, lt$mday - 1L)
  out
}

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
min_valid_years <- function(n_years, cap) {
  check_whole(n_years, "n_years")
  check_whole(cap, "cap", positive = TRUE, single = TRUE)
  as.integer(pmax(1L, pmin(cap, n_years)))
}

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
                         first_year = NULL) {
  check_choice(rule, c("rolling", "fixed", "por"), "rule")
  check_whole(year, "year", single = TRUE)
  year <- as.integer(year)
  given <- c(length = !is.null(length), earliest = !is.null(earliest),
             periods = !is.null(periods), publish_lag = !is.null(publish_lag),
             first_year = !is.null(first_year))
  uses <- switch(rule,
                 rolling = c("length", "earliest"),
                 fixed   = c("periods", "publish_lag"),
                 por     = "first_year")
  extra <- names(given)[given & !names(given) %in% uses]
  if (length(extra))
    stop_arg(extra[1], "must be NULL for rule \"", rule, "\"")
  switch(rule,
    rolling = {
      if (is.null(length)) stop_arg("length", "is required for rule \"rolling\"")
      check_whole(length, "length", positive = TRUE, single = TRUE)
      first <- year - as.integer(length)
      if (!is.null(earliest)) {
        check_whole(earliest, "earliest", single = TRUE)
        first <- max(first, as.integer(earliest))
      }
      if (first > year - 1L) return(NULL)
      first:(year - 1L)
    },
    fixed = {
      if (is.null(periods) || is.null(publish_lag))
        stop_arg(if (is.null(periods)) "periods" else "publish_lag",
                 "is required for rule \"fixed\"")
      if (!is.list(periods) || !all(vapply(periods, \(p) is.numeric(p) && length(p) == 2L, logical(1))))
        stop_arg("periods", "must be a list of c(first, last)")
      check_whole(publish_lag, "publish_lag", single = TRUE)
      last  <- vapply(periods, \(p) as.integer(p[2]), integer(1))
      ok    <- which(last + as.integer(publish_lag) < year)
      if (!length(ok)) return(NULL)
      p <- periods[[ok[which.max(last[ok])]]]
      as.integer(p[1]):as.integer(p[2])
    },
    por = {
      if (is.null(first_year)) stop_arg("first_year", "is required for rule \"por\"")
      check_whole(first_year, "first_year", single = TRUE)
      if (first_year > year - 1L) return(NULL)
      as.integer(first_year):(year - 1L)
    })
}

#' @rdname period_years
#' @export
period_key <- function(first, last, year_type, por = FALSE) {
  check_year_type(year_type)
  check_whole(last, "last")
  if (!is.logical(por) || length(por) != 1L || is.na(por)) stop_arg("por", "must be TRUE or FALSE")
  if (por) return(sprintf("%spor-%d", year_type, as.integer(last)))
  check_whole(first, "first")
  sprintf("%s%d-%d", year_type, as.integer(first), as.integer(last))
}

# A Date from year, month and day vectors; NA where the day does not exist
# (Feb 29 in a non-leap year).
make_date <- function(year, month, day) {
  as.Date(sprintf("%04d-%02d-%02d", as.integer(year), as.integer(month), as.integer(day)),
          optional = TRUE, format = "%Y-%m-%d")
}
