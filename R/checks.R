# Argument checks shared by every exported function. Messages name the
# argument: "normals: `min_valid` must be ...".

stop_arg <- function(arg, ...) {
  stop("normals: `", arg, "` ", ..., call. = FALSE)
}

check_choice <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices)
    stop_arg(arg, "must be one of ", paste0('"', choices, '"', collapse = ", "))
  invisible(x)
}

check_year_type <- function(year_type) check_choice(year_type, c("wy", "cy"), "year_type")

check_date <- function(x, arg, single = FALSE) {
  if (!inherits(x, "Date")) stop_arg(arg, "must be a Date")
  if (single && length(x) != 1L) stop_arg(arg, "must be a single Date")
  invisible(x)
}

check_whole <- function(x, arg, positive = FALSE, single = FALSE, na_ok = FALSE) {
  ok <- is.numeric(x) && (na_ok || !anyNA(x)) &&
    all(is.na(x) | x == round(x)) && (!positive || all(is.na(x) | x >= 1))
  if (!ok || (single && length(x) != 1L))
    stop_arg(arg, "must be ", if (single) "a single " else "",
             if (positive) "positive " else "", "whole number", if (single) "" else "s")
  invisible(x)
}

check_matrix <- function(x, arg) {
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)))
    stop_arg(arg, "must be a numeric matrix")
  invisible(x)
}

# A numeric matrix in double storage, for the C++ (integer NA becomes NA_real_).
as_double_matrix <- function(x) {
  if (!is.double(x)) storage.mode(x) <- "double"
  x
}

# min_valid as the C++ takes it: integer, length 1 or n.
as_min_valid <- function(min_valid, n, arg = "min_valid") {
  if (!is.numeric(min_valid) || anyNA(min_valid) || any(min_valid != round(min_valid)) ||
      !(length(min_valid) %in% c(1L, n)) || any(min_valid < 1))
    stop_arg(arg, "must be a positive whole number, length 1 or one per row")
  as.integer(min_valid)
}
