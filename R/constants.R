#' Quantile probabilities of every normal
#'
#' The 19 probabilities at which normals report quantiles: 0, .01, .02, .025,
#' .05, .1, .2, .25, .3, .5, .7, .75, .8, .9, .95, .975, .98, .99 and 1. Names
#' are `sprintf("q%04d", round(p * 1000))` (`q0000` ... `q1000`), which sort
#' naturally. `unname()` gives the probabilities, `names()` the column names.
#'
#' This is the one definition. It replaces `SNODAS_QUANTILES` (mco-snodas
#' `R/utils.R`, vendored in mco-snotel).
#'
#' @format A named numeric vector of length 19.
#' @export
NORMALS_QUANTILES <- local({
  p <- c(0.0, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5,
         0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.98, 0.99, 1.0)
  stats::setNames(p, sprintf("q%04d", round(p * 1000)))
})

#' Flavor names
#'
#' `NORMALS_FLAVORS` maps the spec's flavor names (the names) to the names
#' published until the layout v1 switch (the values). `NORMALS_FLAVOR_CAPTIONS`
#' gives each flavor's plot caption.
#'
#' @format Named character vectors of length 3.
#' @seealso [flavor_legacy()]
#' @export
NORMALS_FLAVORS <- c(empirical = "ecdf", fitted = "zig", smoothed = "zig_smooth")

#' @rdname NORMALS_FLAVORS
#' @export
NORMALS_FLAVOR_CAPTIONS <- c(empirical = "Empirical History",
                             fitted    = "Fitted History",
                             smoothed  = "Smoothed History")

#' Translate flavor names between the spec and the legacy names
#'
#' @param x Character vector of flavor names.
#' @param to `"legacy"` (`empirical` -> `ecdf`, `fitted` -> `zig`,
#'   `smoothed` -> `zig_smooth`) or `"spec"` (the reverse). No default.
#' @return Character vector, the same length as `x`. Unknown names are an
#'   error.
#' @export
flavor_legacy <- function(x, to) not_implemented()

#' Quantile column names
#'
#' @param probs Probabilities in \[0, 1\].
#' @return `sprintf("q%04d", round(probs * 1000))`. An error if two
#'   probabilities give the same name.
#' @export
quantile_columns <- function(probs = NORMALS_QUANTILES) not_implemented()
