#' Normals and percentile tables
#'
#' `normals_table()` builds every flavor for every date of one normals year,
#' for subjects held as rows of reference matrices (stations, zones, bands).
#' The smoother is fitted at every date except Feb 29 and evaluated at every
#' date.
#'
#' It returns `list(normals, smooth)`. `smooth` is the `normals_smooth` fit
#' (`NULL` without the smoothed flavor). `normals` is a tibble, one row per
#' flavor x date x subject, ordered by flavor (in `flavors` order), date and
#' subject (input order), with columns `flavor`, `date`, `day`, the subject
#' columns, `distribution`, `n_valid`, `ref_first_wy`, `ref_last_wy`, `mean`,
#' the family's parameters, `fit_status` and the quantile columns.
#'
#' `percentiles_table()` places observed values in those normals. It returns a
#' tibble ordered by flavor, date and subject: `flavor`, `date`, the subject
#' columns, `value`, `ptile` and `standardized` (Acklam, to 0.01). Observations
#' of unknown subjects, or on dates without normals, are dropped.
#'
#' @param ref `function(date)` returning that date's reference matrix (rows in
#'   `subjects` order, year column names), or a list of matrices named by
#'   `format(date)`.
#' @param dates The normals year's dates ([year_dates()]).
#' @param subjects Data frame of subject key columns, one row per `ref` row.
#' @param family The family for the fitted and smoothed flavors.
#' @param year_type `"wy"` or `"cy"`.
#' @param min_valid Integer, or `function(ref)` returning one. No default.
#' @param fallback Named character `c(fitted = ., smoothed = .)`, each
#'   `"empirical"` or `"none"`. No default.
#' @param smoothed_mask `"none"` or `"min_valid"`. No default.
#' @param flavors Flavors to build.
#' @param probs Quantile probabilities.
#' @param cores Forked workers for the smoother.
#' @param normals A `normals_table()$normals`.
#' @param obs Data frame of subject key columns, `date` and `value`.
#' @param ptile_digits Decimal places of `ptile`. No default.
#' @export
normals_table <- function(ref, dates, subjects, family, year_type, min_valid,
                          fallback, smoothed_mask,
                          flavors = names(NORMALS_FLAVORS),
                          probs = NORMALS_QUANTILES, cores = 1L) {
  check_year_type(year_type)
  check_date(dates, "dates")
  if (!length(dates) || !identical(as.numeric(dates), as.numeric(year_dates(year_of(dates[1], year_type), year_type))))
    stop_arg("dates", "must be every day of one year of `year_type` (year_dates())")
  if (!is.data.frame(subjects) || !ncol(subjects)) stop_arg("subjects", "must be a data frame of key columns")
  n <- nrow(subjects)
  flavors <- check_flavors(flavors)
  q_names <- quantile_columns(probs)
  probs <- unname(probs)
  model <- setdiff(flavors, "empirical")
  if (length(model)) {
    check_family(family)
    if (is.null(family$fit)) stop_arg("family", "must have a fit for the fitted and smoothed flavors")
    if ("smoothed" %in% model && is.null(family$smooth)) stop_arg("family", "has no smoother")
  }
  fallback <- check_fallback(fallback, model)
  check_choice(smoothed_mask, c("none", "min_valid"), "smoothed_mask")
  check_whole(cores, "cores", positive = TRUE, single = TRUE)
  check_min_valid_arg(min_valid)
  ref_at <- ref_getter(ref, n)
  mv_of  <- \(h) as_min_valid(if (is.function(min_valid)) min_valid(h) else min_valid, n)

  # Per date, per flavor: a block of n rows.
  block <- function(flavor, d, h, dist, params, status = NULL) {
    mean_v <- dist_mean(dist)
    q <- dist_quantile(dist, probs)
    # The family where its parameters are complete; "empirical" for the
    # empirical flavor and for rows that fell back; NA without a normal.
    has <- !is.na(mean_v)
    distribution <- rep(NA_character_, n)
    if (is.null(params)) {
      distribution[has] <- "empirical"
    } else {
      complete <- stats::complete.cases(params)
      distribution[has & complete]  <- family$name
      distribution[has & !complete] <- if (dist$fallback == "empirical") "empirical" else family$name
    }
    span <- ref_span_or_na(h)
    list(flavor = flavor, date = d, distribution = distribution,
         n_valid = as.integer(rowSums(!is.na(h))), ref_first_wy = span$first,
         ref_last_wy = span$last, mean = mean_v, params = params,
         fit_status = if (is.null(status)) rep(NA_character_, n) else fit_status_label(status, family),
         q = q)
  }

  blocks <- list(empirical = NULL, fitted = NULL, smoothed = NULL)
  daily <- intersect(c("empirical", "fitted"), flavors)
  if (length(daily)) {
    per <- lapply(seq_along(dates), \(j) {
      d <- dates[j]; h <- ref_at(d); mv <- mv_of(h)
      out <- list()
      if ("empirical" %in% daily)
        out$empirical <- block("empirical", d, h,
                               normals_dist(family_empirical(), NULL, h, mv, NULL), NULL)
      if ("fitted" %in% daily) {
        params <- fit_day(day_stats(h, family), family, mv)
        out$fitted <- block("fitted", d, h,
                            normals_dist(family, params, h, mv, fallback[["fitted"]]), params)
      }
      out
    })
    for (fl in daily) blocks[[fl]] <- lapply(per, `[[`, fl)
  }

  smooth <- NULL
  if ("smoothed" %in% flavors) {
    # Feb 29 is left out of the fit and evaluated between Feb 28 and Mar 1.
    fit_dates <- dates[format(dates, "%m-%d") != "02-29"]
    stats <- stack_stats(lapply(fit_dates, \(d) day_stats(ref_at(d), family)))
    smooth <- smooth_fit(stats, cycle_day(fit_dates, year_type), family, cores)
    ev <- smooth_eval(smooth, cycle_day(dates, year_type))
    blocks$smoothed <- lapply(seq_along(dates), \(j) {
      d <- dates[j]; h <- ref_at(d); mv <- mv_of(h)
      params <- smooth_params(ev, j)
      if (smoothed_mask == "min_valid") params <- mask_sparse(params, h, mv)
      block("smoothed", d, h, normals_dist(family, params, h, mv, fallback[["smoothed"]]),
            params, smooth$status)
    })
  }

  list(normals = bind_blocks(unname(unlist(blocks[flavors], recursive = FALSE)), subjects,
                             if (length(model)) family$params else character(), q_names),
       smooth = smooth)
}

# One tibble from per-date blocks, in block order.
bind_blocks <- function(blocks, subjects, param_names, q_names) {
  n <- nrow(subjects)
  col <- \(name) do.call(c, lapply(blocks, `[[`, name))
  keys <- lapply(as.list(subjects), \(v) rep(v, length(blocks)))
  params <- lapply(seq_along(param_names), \(i)
    unlist(lapply(blocks, \(b) if (is.null(b$params)) rep(NA_real_, n) else unname(b$params[, i]))))
  q <- do.call(rbind, lapply(blocks, `[[`, "q"))
  tibble::as_tibble(c(
    list(flavor = rep(vapply(blocks, `[[`, "", "flavor"), each = n),
         date   = rep(do.call(c, lapply(blocks, `[[`, "date")), each = n),
         day    = rep(format(do.call(c, lapply(blocks, `[[`, "date")), "%m-%d"), each = n)),
    keys,
    list(distribution = col("distribution"), n_valid = col("n_valid"),
         ref_first_wy = col("ref_first_wy"), ref_last_wy = col("ref_last_wy"),
         mean = col("mean")),
    stats::setNames(params, param_names),
    list(fit_status = col("fit_status")),
    stats::setNames(lapply(seq_along(q_names), \(i) unname(q[, i])), q_names)
  ))
}

#' @rdname normals_table
#' @export
percentiles_table <- function(normals, ref, subjects, obs, family, min_valid,
                              fallback, ptile_digits,
                              flavors = names(NORMALS_FLAVORS)) {
  if (!is.data.frame(subjects) || !ncol(subjects)) stop_arg("subjects", "must be a data frame of key columns")
  keys <- names(subjects)
  n <- nrow(subjects)
  flavors <- check_flavors(flavors)
  model <- setdiff(flavors, "empirical")
  if (length(model)) check_family(family)
  fallback <- check_fallback(fallback, model)
  check_whole(ptile_digits, "ptile_digits", single = TRUE)
  check_min_valid_arg(min_valid)
  if (!is.data.frame(normals) || !all(c("flavor", "date", keys) %in% names(normals)))
    stop_arg("normals", "must be a normals_table()$normals with the subject columns")
  if (length(model) && !all(family$params %in% names(normals)))
    stop_arg("normals", "must hold the family's parameter columns")
  if (!is.data.frame(obs) || !all(c(keys, "date", "value") %in% names(obs)))
    stop_arg("obs", "must have the subject columns, `date` and `value`")
  check_date(obs$date, "obs$date")
  ref_at <- ref_getter(ref, n)

  # As mco-snotel calc_station_percentiles(): observations of unknown
  # subjects, or on dates without normals, are dropped.
  subj_key <- key_of(subjects)
  obs <- obs[key_of(obs[keys]) %in% subj_key, , drop = FALSE]
  norm_key <- key_of(normals[keys])
  by_date <- split(seq_len(nrow(normals)), format(normals$date))
  out <- list()
  for (d in sort(unique(obs$date))) {
    d <- as.Date(d, origin = "1970-01-01")
    rows_d <- by_date[[format(d)]]
    if (is.null(rows_d)) next
    o <- obs[obs$date == d, , drop = FALSE]
    o_key <- key_of(o[keys])
    i <- match(o_key, subj_key)
    h_all <- ref_at(d)
    mv_all <- as_min_valid(if (is.function(min_valid)) min_valid(h_all) else min_valid, n)
    h  <- h_all[i, , drop = FALSE]
    mv <- if (length(mv_all) == 1L) mv_all else mv_all[i]
    for (fl in flavors) {
      nf <- rows_d[normals$flavor[rows_d] == fl]
      nf <- nf[match(o_key, norm_key[nf])]
      p <- if (fl == "empirical") {
        dist_cdf(normals_dist(family_empirical(), NULL, h, mv, NULL), o$value)
      } else {
        params <- as.matrix(as.data.frame(normals[nf, family$params, drop = FALSE]))
        dist_cdf(normals_dist(family, params, h, mv, fallback[[fl]]), o$value)
      }
      out[[length(out) + 1L]] <- list(flavor = fl, date = d, o = o, i = i, p = p)
    }
  }
  empty <- tibble::as_tibble(c(list(flavor = character(), date = as.Date(character())),
                               lapply(as.list(subjects), \(v) v[0]),
                               list(value = numeric(), ptile = numeric(), standardized = numeric())))
  if (!length(out)) return(empty)
  ord <- order(match(vapply(out, `[[`, "", "flavor"), flavors), vapply(out, \(b) as.numeric(b$date), 0))
  res <- tibble::as_tibble(c(
    list(flavor = unlist(lapply(out, \(b) rep(b$flavor, nrow(b$o)))),
         date   = do.call(c, lapply(out, \(b) rep(b$date, nrow(b$o))))),
    lapply(stats::setNames(keys, keys), \(k) do.call(c, lapply(out, \(b) b$o[[k]]))),
    list(value        = unlist(lapply(out, \(b) b$o$value)),
         ptile        = unlist(lapply(out, \(b) as_ptile(b$p, ptile_digits))),
         standardized = unlist(lapply(out, \(b) standardize(b$p, "acklam", "z"))),
         .i           = unlist(lapply(out, `[[`, "i")),
         .b           = rep(order(ord), vapply(out, \(b) nrow(b$o), 0L)))
  ))
  # Ordered by flavor, date and subject (input order); ties keep obs order.
  res <- res[order(res$.b, res$.i), , drop = FALSE]
  res$.i <- NULL
  res$.b <- NULL
  res
}

check_flavors <- function(flavors) {
  if (!is.character(flavors) || !length(flavors) || anyDuplicated(flavors) ||
      !all(flavors %in% names(NORMALS_FLAVORS)))
    stop_arg("flavors", "must be distinct flavors among ",
             paste0('"', names(NORMALS_FLAVORS), '"', collapse = ", "))
  names(NORMALS_FLAVORS)[names(NORMALS_FLAVORS) %in% flavors][order(match(
    names(NORMALS_FLAVORS)[names(NORMALS_FLAVORS) %in% flavors], flavors))]
}

check_fallback <- function(fallback, model) {
  if (!length(model)) return(fallback)
  if (!is.character(fallback) || is.null(names(fallback)) || !all(model %in% names(fallback)) ||
      !all(fallback[model] %in% c("empirical", "none")))
    stop_arg("fallback", "must be c(", paste0(model, ' = "empirical" or "none"', collapse = ", "), ")")
  fallback
}

check_min_valid_arg <- function(min_valid) {
  if (!is.function(min_valid) && !is.numeric(min_valid))
    stop_arg("min_valid", "must be a positive integer or a function")
  invisible(min_valid)
}

# ref as a function of the date: a function, or a list named by format(date).
ref_getter <- function(ref, n) {
  get <- if (is.function(ref)) ref
         else if (is.list(ref) && !is.null(names(ref))) \(d) ref[[format(d)]]
         else stop_arg("ref", "must be a function(date) or a list of matrices named by date")
  function(d) {
    h <- get(d)
    if (is.null(h)) stop_arg("ref", "has no reference matrix for ", format(d))
    if (!is.matrix(h) || nrow(h) != n)
      stop_arg("ref", "must give one row per subject (", n, ") for ", format(d))
    h
  }
}

# reference_span() where the columns are years; NA otherwise (e.g. dates).
ref_span_or_na <- function(h) {
  cn <- colnames(h)
  if (!is.null(cn) && all(grepl("^[0-9]+$", cn))) return(reference_span(h))
  na <- rep(NA_integer_, nrow(h))
  list(first = na, last = na)
}

# One string per row of a key data frame.
key_of <- function(x) {
  do.call(paste, c(lapply(unname(as.list(x)), as.character), sep = "\r"))
}
