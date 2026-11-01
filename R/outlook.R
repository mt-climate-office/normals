#' The outlook: the likely value for the rest of the year
#'
#' Today's value is placed in its day's distribution (z0, clamped), and that
#' percentile is carried forward on the standard-normal scale, fading toward
#' the normals:
#'
#' * rho(h) = exp(-(h / length)^alpha), h days ahead (`outlook_persistence()`);
#' * k = min(exp(-beta rho min(z0, 0)), k_max) widens the spread below normal
#'   (`outlook_spread()`);
#' * q_p = Q_target(pnorm(rho z0 + qnorm(p) k sqrt(1 - rho^2))) x anchor^rho
#'   (`outlook_quantiles()`), with the anchor scaling a value beyond z = +/-4
#'   (`outlook_start()`).
#'
#' For zero-mass families, a fitted p_zero below `p_zero_floor` / n_valid is
#' taken as zero for subjects with a positive value today, and
#' `outlook_p_zero()` gives the chance of exactly zero on each target date.
#'
#' The alpha table (persistence shape by cycle day) and beta (spread) are per
#' dataset and are passed in: `outlook_alpha_at()` interpolates an alpha table
#' the caller has already filtered to one period and flavor. `outlook_lengths()`
#' fits the persistence length per subject and date from the reference years.
#' `outlook()` composes the pieces for one issue date. `outlook_crps()` and
#' `outlook_pit()` score hindcasts.
#'
#' @param k_max,fit_lags,pool_days,min_pairs,min_lags,z_max,l_grid,p_zero_floor,anchor_range,rho_floor
#'   The method's constants (today's values).
#' @param quantile_digits,p_zero_digits,length_digits,alpha_digits Roundings
#'   of the published columns.
#' @param control An `outlook_control()` list.
#' @param table Data frame with `day` (cycle day) and `alpha`.
#' @param dates,target_dates `Date` vectors.
#' @param year_type `"wy"` or `"cy"`.
#' @param lead Days ahead.
#' @param length Persistence length (days), per subject or scalar.
#' @param alpha Persistence shape, per subject or scalar (`outlook()`), or one
#'   per date (`outlook_lengths()`).
#' @param z0 Starting z.
#' @param rho Persistence.
#' @param beta Spread factor.
#' @param issue The issue date's `normals_dist`.
#' @param target A target date's `normals_dist`.
#' @param targets List of `normals_dist`, one per target date.
#' @param x Values on the issue date, length n.
#' @param start An `outlook_start()` result.
#' @param probs Quantile probabilities.
#' @param issue_date The issue date.
#' @param subjects Data frame of subject key columns.
#' @param dists List of `normals_dist`, one per date of the normals year.
#' @param q n x m quantiles at `probs`.
#' @param y Outcomes, length n.
#' @param u Uniform draws for zero outcomes, length n.
#' @name outlook
NULL

#' @rdname outlook
#' @export
outlook_control <- function(k_max = 3,
                            fit_lags = c(3, 7, 14, 21, 30, 45, 60, 90, 120),
                            pool_days = 15L, min_pairs = 8L, min_lags = 3L,
                            z_max = 4,
                            l_grid = exp(seq(log(2), log(3000), length.out = 400)),
                            p_zero_floor = 0.5, anchor_range = c(0.1, 10),
                            rho_floor = 1e-12, quantile_digits = 1,
                            p_zero_digits = 3, length_digits = 2,
                            alpha_digits = 3) not_implemented()

#' @rdname outlook
#' @export
outlook_alpha_at <- function(table, dates, year_type) not_implemented()

#' @rdname outlook
#' @export
outlook_persistence <- function(lead, length, alpha) not_implemented()

#' @rdname outlook
#' @export
outlook_spread <- function(z0, rho, beta, control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook_start <- function(issue, x, control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook_quantiles <- function(target, start, rho, beta,
                              probs = NORMALS_QUANTILES,
                              control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook_p_zero <- function(target, start, rho, beta,
                           control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook <- function(issue, targets, x, issue_date, target_dates, subjects,
                    length, alpha, beta, probs = NORMALS_QUANTILES,
                    control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook_lengths <- function(dists, dates, subjects, alpha,
                            control = outlook_control()) not_implemented()

#' @rdname outlook
#' @export
outlook_crps <- function(q, y, probs) not_implemented()

#' @rdname outlook
#' @export
outlook_pit <- function(target, start, rho, beta, y, u,
                        control = outlook_control()) not_implemented()
