#' normals: climate normals from empirical, fitted and smoothed distributions
#'
#' Daily climate normals for any set of subjects (stations, zones or pixels).
#' Three flavors: `empirical` (the zero-inflated type-7 ECDF of the reference
#' years), `fitted` (a distribution family fitted for each day) and `smoothed`
#' (the same family with its parameters smoothed through the year). Families are
#' plug-ins ([new_family()]); the first is the zero-inflated gamma
#' ([family_zig()]).
#'
#' The package reads and writes no storage and builds no paths. Its inputs and
#' outputs are matrices, vectors and tibbles. The API is described in
#' `docs/API.md` in the source repository.
#'
#' @keywords internal
#' @importFrom stats approx pnorm qnorm qgamma pgamma plogis median
#' @importFrom utils read.csv
#' @importFrom grDevices dev.size rgb
#' @importFrom parallel mclapply
#' @importFrom mgcv smoothCon
#' @importFrom tibble tibble
#' @importFrom grid makeContent
"_PACKAGE"

# Every exported function is a stub until the API (docs/API.md) is approved.
not_implemented <- function() stop("not implemented", call. = FALSE)
