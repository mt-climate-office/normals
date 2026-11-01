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
record_matrix <- function(x, subject, date, value, subjects = NULL) not_implemented()

#' @rdname record
#' @export
record_at <- function(record, dates) not_implemented()

#' @rdname record
#' @export
reference_record <- function(record, date, years, year_type) not_implemented()

#' @rdname record
#' @export
reference_span <- function(ref) not_implemented()
