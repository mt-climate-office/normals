# Golden fixtures (data-raw/make-fixtures.R): today's pipeline outputs, which
# the 2.0 calls must reproduce exactly.

golden <- function(name) {
  path <- test_path("fixtures", paste0(name, ".rds"))
  if (!file.exists(path)) skip(paste("fixture", name, "not found"))
  readRDS(path)
}

# TRUE when this R session is the fixtures' platform (R version, OS
# architecture, BLAS and LAPACK).
golden_platform <- local({
  same <- NULL
  function() {
    if (is.null(same)) {
      fp <- readRDS(test_path("fixtures", "platform.rds"))
      si <- utils::sessionInfo()
      same <<- identical(fp, list(R = R.version.string, platform = R.version$platform,
                                  blas = si$BLAS, lapack = si$LAPACK))
    }
    same
  }
})

# Largest relative difference between the doubles of two objects that are
# otherwise alike (same structure, types, NA and NaN positions); NA if they
# are not alike.
golden_rel_diff <- function(new, old) {
  if (is.list(old)) {
    if (!is.list(new) || !identical(names(new), names(old)) || length(new) != length(old) ||
        !identical(class(new), class(old))) return(NA_real_)
    d <- vapply(seq_along(old), \(i) golden_rel_diff(new[[i]], old[[i]]), numeric(1))
    return(if (!length(d)) 0 else max(d))
  }
  if (!identical(attributes(new), attributes(old)) || typeof(new) != typeof(old) ||
      length(new) != length(old)) return(NA_real_)
  if (!is.double(old)) return(if (identical(new, old)) 0 else NA_real_)
  if (!identical(is.na(new), is.na(old)) || !identical(is.nan(new), is.nan(old))) return(NA_real_)
  ok <- !is.na(old)
  if (!any(ok)) return(0)
  a <- unclass(new)[ok]; b <- unclass(old)[ok]
  if (!identical(is.infinite(a), is.infinite(b)) || any(a[is.infinite(a)] != b[is.infinite(b)]))
    return(NA_real_)
  f <- is.finite(b)
  max(0, abs(a[f] - b[f]) / pmax(abs(b[f]), .Machine$double.xmin))
}

# Exact on the fixtures' platform: identical(), which tells NA from NaN and
# integer from double. Elsewhere (API.md section 7) a result may differ in its
# last bits through linear algebra or libm: identical() is tried first, then
# the largest relative difference is reported and must not exceed 1e-12.
expect_golden <- function(new, old, label = deparse(substitute(new))) {
  if (identical(new, old)) return(expect_identical(new, old))
  if (golden_platform()) return(expect_identical(new, old, label = label))
  d <- golden_rel_diff(new, old)
  if (!is.na(d)) message(sprintf("golden: %s differs by %.3g (relative) off the fixture platform", label, d))
  expect(!is.na(d) && d <= 1e-12,
         sprintf("%s: not identical to the fixture, largest relative difference %s", label,
                 if (is.na(d)) "n/a (structure, type or NA pattern differs)" else format(d)))
}
