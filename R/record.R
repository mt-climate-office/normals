#' Records and reference matrices
#'
#' * `record_matrix()`: a subject x day matrix from a long table, for fast
#'   lookups by date. Returns a `normals_record`: `list(m, subjects, origin)`.
#'   Duplicate subject x date rows are an error.
#' * `record_at()`: values on the given dates, n_subjects x `length(dates)`,
#'   `NA` outside the record.
#' * `reference_record()`: one day's reference matrix, n_subjects x
#'   `length(years)`, columns named by year.
#' * `reference_span()`: first and last column year with a value, per row
#'   (`NA` for an empty row). Needs year column names.
#'
#' @param x A data frame.
#' @param subject,date,value Column names of `x` (strings). No defaults.
#' @param subjects Row order (default: the sorted unique subjects).
#' @param record A `normals_record`.
#' @param dates A `Date` vector.
#' @param years Integer reference years.
#' @param year_type `"wy"` or `"cy"`.
#' @param ref A reference matrix with year column names.
#' @name record
NULL

#' @rdname record
#' @export
record_matrix <- function(x, subject, date, value, subjects = NULL) {
  if (!is.data.frame(x)) stop_arg("x", "must be a data frame")
  for (a in c("subject", "date", "value")) {
    col <- get(a)
    if (!is.character(col) || length(col) != 1L || !col %in% names(x))
      stop_arg(a, "must name a column of `x`")
  }
  s <- x[[subject]]; d <- x[[date]]; v <- x[[value]]
  check_date(d, "date")
  if (!is.numeric(v)) stop_arg("value", "must name a numeric column")
  if (is.null(subjects)) subjects <- sort(unique(s))
  if (anyDuplicated(subjects)) stop_arg("subjects", "must not repeat a subject")
  if (!nrow(x)) stop_arg("x", "has no rows")
  # As mco-snotel swe_matrix(): the origin is the first date of x (every
  # subject), and values of subjects not in `subjects` are dropped.
  origin <- min(d)
  n_day  <- as.integer(max(d) - origin) + 1L
  m <- matrix(NA_real_, length(subjects), n_day)
  i <- match(s, subjects)
  ok <- !is.na(i)
  j <- as.integer(d[ok] - origin) + 1L
  if (anyDuplicated(cbind(i[ok], j)))
    stop_arg("x", "has more than one row for a subject and date")
  m[cbind(i[ok], j)] <- v[ok]
  structure(list(m = m, subjects = subjects, origin = origin), class = "normals_record")
}

#' @rdname record
#' @export
record_at <- function(record, dates) {
  check_record(record)
  check_date(dates, "dates")
  j <- as.integer(dates - record$origin) + 1L
  j[j < 1L | j > ncol(record$m)] <- NA_integer_
  record$m[, j, drop = FALSE]
}

#' @rdname record
#' @export
reference_record <- function(record, date, years, year_type) {
  h <- record_at(record, reference_dates(date, years, year_type))
  colnames(h) <- years
  h
}

#' @rdname record
#' @export
reference_span <- function(ref) {
  check_matrix(ref, "ref")
  years <- suppressWarnings(as.integer(colnames(ref)))
  if (is.null(colnames(ref)) || anyNA(years) || any(colnames(ref) != as.character(years)))
    stop_arg("ref", "must have years as column names")
  # As mco-snotel record_span().
  ok <- !is.na(ref)
  any_ok <- rowSums(ok) > 0
  first <- ifelse(any_ok, years[max.col(ok, ties.method = "first")], NA_integer_)
  last  <- ifelse(any_ok, years[max.col(ok, ties.method = "last")],  NA_integer_)
  list(first = as.integer(first), last = as.integer(last))
}

check_record <- function(record) {
  if (!inherits(record, "normals_record")) stop_arg("record", "must be a normals_record (record_matrix())")
  invisible(record)
}
