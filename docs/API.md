# normals 2.0: API design (step 1A-0)

**Status:** proposal for gate G1a, 2026-11-01. Nothing here is implemented; every function in
`R/` is a stub that stops with "not implemented".

**Governs:** steps 1A-i to 1A-v (the port), 1B (mco-snodas) and 1C (mco-snotel) of the migration
plan. The spec is mco-snotel `outlook:docs/bucket-layout.md` ("the spec"); §-numbers below are
its sections.

**Baselines read:** mco-snodas `main` at 9b45d74 (`src/` and the vendored R unchanged since
62fe19a), mco-snotel `origin/main` at f39fe56, and mco-snotel `outlook` at bb8d06a for the
outlook. File:line references are to those commits.

## Contents

1. [Principles](#1-principles)
2. [Vocabulary and conventions](#2-vocabulary-and-conventions)
3. [Exported API](#3-exported-api)
4. [The family plug-in interface](#4-the-family-plug-in-interface)
5. [Differences between the pipelines, and the arguments that encode them](#5-differences-between-the-pipelines-and-the-arguments-that-encode-them)
6. [Mapping from today's code to 2.0 calls](#6-mapping-from-todays-code-to-20-calls)
7. [Golden-fixture test plan](#7-golden-fixture-test-plan)
8. [Package layout, dependencies and CI](#8-package-layout-dependencies-and-ci)
9. [Open questions](#9-open-questions)
10. [Surprises found in the code](#10-surprises-found-in-the-code)

---

## 1. Principles

1. **Phase 1 changes no output.** Every behavior on which mco-snodas and mco-snotel differ is an
   argument with **no default**. Each pipeline passes the values that reproduce it (§5), and the
   golden fixtures (§7) prove it bit for bit. Shared constants (the L-moment minimum of 4
   positive values, the spline's k = 24 and λ, …) keep today's values as defaults, because
   no pipeline differs on them.
2. **No storage, no paths, no AWS** (spec §2 rule 3). Inputs and outputs are matrices, vectors
   and tibbles. The package never reads a file except its own `inst/` assets through
   `system.file()`.
3. **Stations, zones and pixels are all "subjects".** Two levels of API:
   - *matrix level* (`day_stats()`, `fit_day()`, `smooth_fit()`, `normals_dist()`, `dist_*()`):
     works on any block of rows. This is what the gridded SNODAS build calls, block by block.
   - *table level* (`normals_table()`, `percentiles_table()`, `outlook()`,
     `outlook_lengths()`): long tibbles in the §8 schemas, for stations and zones.
4. **Spec names inside, legacy names at the edge.** The package says `empirical` / `fitted` /
   `smoothed` and names family parameters plainly (`shape`, `scale`, `p_zero`). In Phase 1 the
   pipelines rename at their edge with `flavor_legacy()` and a column rename (§6.4), so published
   names stay `ecdf` / `zig` / `zig_smooth`, `zig_shape`, `swe_ptile`, … until the v1 switch.
5. **The C++ moves, it is not rewritten.** 1A-i moves `ecdf_impl.h`, `ecdf.cpp`, `zig.cpp`,
   `zig_spline.cpp` and `qnorm_grid.cpp` into `src/` with only these changes: Rcpp exports
   become internal (`_normals_*`), and `MIN_VALID_YEARS` / `min_valid_years(n_time)` is replaced
   by an explicit `min_valid` argument (§5, D3). Arithmetic and evaluation order stay as they are.

## 2. Vocabulary and conventions

| Term | Meaning |
|---|---|
| subject | One row: a station, a zone (or zone × elevation band), or a pixel. |
| `ref` | A **reference matrix**: numeric, n_subjects × n_years, the values of one calendar day in each reference year. `NA` (or `NaN`) is a missing year. Column names are the years (`"2004"`, …) where known; `reference_record()` sets them. |
| normals year | The year whose normals are built: water year W (`year_type = "wy"`, Oct 1–Sep 30, named for the year it ends in) or calendar year (`"cy"`). Always explicit; there is no default `year_type`. |
| `dates` | The normals year's days, `year_dates(year, year_type)`: 365 or 366 `Date`s. |
| `day` | `MM-DD`, `02-29` included (spec §3). |
| cycle day | Position on the 365-day cycle used by smoothing: `wy` Oct 1 = 1 … Sep 30 = 365, Feb 29 = 151.5; `cy` Jan 1 = 1 … Dec 31 = 365, Feb 29 = 59.5. |
| family | A distribution family object (§4), e.g. `family_zig()` (`"zero_inflated_gamma"`) or `family_empirical()` (`"empirical"`). |
| flavor | How a distribution is estimated (spec §3): `empirical`, `fitted`, `smoothed`. |
| `params` | n_subjects × n_params numeric matrix, columns named by `family$params`. |
| `min_valid` | Minimum number of non-missing reference years for a subject to get a normal on a day. Always passed; see §5 D3. |
| `fallback` | What a fitted or smoothed distribution does on rows where the family's parameters are missing: `"empirical"` (use the zero-inflated type-7 ECDF of `ref`) or `"none"` (return `NA`). Always passed; see §5 D1. |

**Errors.** Every exported function validates its inputs and stops with a message naming the
argument (`"normals: `min_valid` must be a positive integer or a function"`), never with an
internal Rcpp error. Shapes that must agree (rows of `params` and `ref`, `length(x)` and rows,
`length(at)` and columns) are checked first. Missing required arguments fail at the call
(they have no defaults).

**Types.** Numeric results are `double`. Counts (`n_valid`, years) are `integer`. Tables are
tibbles, with no grouping and no row names. Integer storage encodings (INT2S COGs, integer
percentiles) are the pipeline's job; the package returns unrounded doubles except where a
rounding is part of a documented output (`as_ptile()`, `standardize()`, `outlook()`).

## 3. Exported API

Signatures are final for 2.0.0 unless the review changes them. "Today" names the code each
function absorbs; §6 has the full mapping.

### 3.1 Constants and vocabulary (`R/constants.R`)

| Object | Value / signature | Contract |
|---|---|---|
| `NORMALS_QUANTILES` | Named numeric, 19 values: 0, .01, .02, .025, .05, .1, .2, .25, .3, .5, .7, .75, .8, .9, .95, .975, .98, .99, 1, named `sprintf("q%04d", round(p * 1000))` (`q0000` … `q1000`). | **The** definition of the quantile probabilities (today `SNODAS_QUANTILES`, mco-snodas `R/utils.R:33-37`, vendored in mco-snotel). `unname()` gives the probabilities, `names()` the column names. Defined, not a stub. |
| `NORMALS_FLAVORS` | `c(empirical = "ecdf", fitted = "zig", smoothed = "zig_smooth")` | Spec name → Phase 1 published name. Defined. |
| `NORMALS_FLAVOR_CAPTIONS` | `c(empirical = "Empirical History", fitted = "Fitted History", smoothed = "Smoothed History")` | Plot captions (spec §3; today `SNODAS_PLOT_FLAVORS` / `SNOTEL_PLOT_FLAVORS`). Defined. |
| `flavor_legacy(x, to)` | `to = "legacy"` or `"spec"`, no default | Maps flavor names either way; unknown names are an error. |
| `quantile_columns(probs = NORMALS_QUANTILES)` | → character | `q%04d` column names for any probabilities (error if two probabilities share a name). |

### 3.2 Calendar and periods (`R/calendar.R`)

`year_type` is `"wy"` or `"cy"` everywhere, never defaulted.

| Function | Returns | Contract / today |
|---|---|---|
| `year_of(date, year_type)` | integer | Water year (`month >= 10` → year + 1) or calendar year. Today `water_year()`. |
| `year_start(year, year_type)` | Date | Oct 1 of `year - 1`, or Jan 1. Today `water_year_start()` (which takes a date; use `year_start(year_of(d, t), t)`). |
| `year_dates(year, year_type)` | Date vector | Every day of the year, 365 or 366. Today `normals_dates()`. |
| `dates_through(date, year_type)` | Date vector | Start of `date`'s year through `date`. Today `water_year_dates()`. |
| `cycle_day(date, year_type)` | double | §2 cycle position. Today `normals_cycle_day()` (wy only). The `cy` form is new. |
| `day_key(date)` | character | `format(date, "%m-%d")`. |
| `reference_dates(date, years, year_type)` | Date vector, `length(years)` | The same calendar day as `date` in each year of `years`, shifted by whole years from `date`'s own year; Feb 29 → Feb 28 in non-leap years. Today mco-snotel `reference_dates()` (`R/normals.R:50-53`), equivalent to mco-snodas `normals_reference_years()` (`R/pipeline_sync.R:231-234`) before its availability filter. Base R, no lubridate. |
| `min_valid_years(n_years, cap)` | integer | `max(1L, min(cap, n_years))`, elementwise. Today `min_valid_years()` in `src/ecdf_impl.h:16-18` with `cap = 10`. `cap` has no default. |
| `period_years(rule, year, length = NULL, earliest = NULL, periods = NULL, publish_lag = NULL, first_year = NULL)` | integer vector, or `NULL` when the role is absent for `year` | Evaluates one spec §5 rule for normals year `year`: `"rolling"` (`length` years ending `year - 1`, from `earliest` if set), `"fixed"` (the latest of `periods`, a list of `c(first, last)`, whose `last + publish_lag < year`), `"por"` (`first_year`:`year - 1`). Arguments a rule does not use must be `NULL`. Today mco-snotel `period_years()` (`R/normals.R:34-40`). See open question Q3. |
| `period_key(first, last, year_type, por = FALSE)` | character | `"wy2004-2026"`, `"cy1996-2025"`, `"wypor-2026"` (spec §3). |

### 3.3 Records (`R/record.R`)

For subjects stored as a long table (stations). Zones and pixels arrive as matrices already.

| Function | Returns | Contract / today |
|---|---|---|
| `record_matrix(x, subject, date, value, subjects = NULL)` | `normals_record`: `list(m, subjects, origin)` | `x` a data frame; `subject`, `date`, `value` are its column names (strings, no defaults). `m` is n_subjects × n_days, `NA` where absent. `subjects` sets the row order (default: sorted unique). Error on duplicate subject × date. Today `swe_matrix()` (`R/normals.R:59-67`, which hard-codes `triplet` / `date` / `wteq`). |
| `record_at(record, dates)` | matrix n_subjects × `length(dates)` | `NA` outside the record. Today `swe_at()`. |
| `reference_record(record, date, years, year_type)` | `ref` matrix, colnames = `years` | `record_at(record, reference_dates(date, years, year_type))`. Today `reference_record()`. |
| `reference_span(ref)` | `list(first, last)` integer vectors | First and last column year with a value, per row; `NA` for an empty row. Needs year colnames. Today `record_span()`. |

### 3.4 Families (`R/family.R`; interface in §4)

| Function | Returns | Contract |
|---|---|---|
| `family_zig(control = zig_control())` | `normals_family` | Zero-inflated gamma, `name = "zero_inflated_gamma"`, params `shape`, `scale`, `p_zero`; the ZIG fit, smoother, CDF, quantile and mean of today's `src/zig.cpp`, `src/zig_spline.cpp` and `R/calc_snodas_zig_spline_normals.R`. |
| `zig_control(min_pos_fit = 4L, k = 24L, lambda = c(nu = 100, l1 = 1000, lcv = 1000), lcv_range = c(0.001, 0.999), min_pos = 10L, min_days = 30L, maxit = 100L, tol = 1e-6)` | list | Today's constants: the L-moment minimum (`zig.cpp:32`), `ZIG_SMOOTH_K`, `ZIG_SMOOTH_LAMBDA`, `ZIG_SMOOTH_LCV_RANGE`, and `zig_spline_fit()`'s `min_pos`, `min_days`, `maxit`, `tol`. Identical in both pipelines, so defaulted. Recorded in provenance by the pipeline. |
| `family_empirical()` | `normals_family` | `name = "empirical"`, no parameters, no smoother; CDF, quantile and mean from `ref` by the zero-inflated type-7 ECDF (`src/ecdf.cpp`). |
| `new_family(name, params, zero_mass, stats_names, stats, fit, smooth, cdf, quantile, mean, p_zero, control = list())` | `normals_family` | Constructor for plug-ins; runs `validate_family()`. |
| `validate_family(family)` | `family`, invisibly | Checks §4's contract (slot types, function formals, `params` non-empty iff `fit` is given, `p_zero` given iff `zero_mass`). |
| `is_family(x)` | logical | |
| `register_family(family)` / `get_family(name)` | invisible family / family | A session registry so provenance's `distribution` names resolve to objects. `zero_inflated_gamma` and `empirical` are registered on load. Error on an unknown name or on re-registering a name with a different object. |

### 3.5 Fitting and smoothing, matrix level (`R/fit.R`, `R/smooth.R`)

| Function | Returns | Contract / today |
|---|---|---|
| `day_stats(ref, family)` | n × s matrix, colnames `family$stats_names` | The per-day statistics the fit and the smoother read. ZIG: `n` (valid years), `z` (zero years), `l1`, `lcv` (NA with < `min_pos_fit` positive values or a degenerate sample). Today `zig_lmom()` / `zig_spline_suffstats()`. Error for a family without `stats`. |
| `fit_day(stats, family, min_valid)` | `params` matrix | One day's parameters from `day_stats()` output. All `NA` where `n < min_valid`; for ZIG, `shape` and `scale` `NA` without L-moments. Today `zig_params_from_lmom(lmom, n_time)` and `zig_fit(mat)` (= `fit_day(day_stats(mat, f), f, min_valid_years(ncol(mat), 10))`). |
| `smooth_basis(at, k)` | `list(X, S, k)` | Cyclic cubic P-spline basis (mgcv `bs = "cp"`, knots 0.5–365.5, `absorb.cons = FALSE`) at cycle positions `at`, and its penalty. Today `zig_spline_basis()`. |
| `stack_stats(stats_by_day)` | named list of n × d matrices | Turns a list of `day_stats()` matrices (one per day) into the stat-major form `smooth_fit()` takes. |
| `smooth_fit(stats, at, family, cores = 1L)` | `normals_smooth` (below) | Fits the family's smoother to stat-major daily statistics (`stats$n`, `stats$z`, … each n × `length(at)`) at cycle positions `at`. Rows are split into `cores` contiguous chunks and fitted with `parallel::mclapply()` (fork); the result does not depend on `cores`. A worker that dies is an error naming the count. Today `zig_spline_fit_stats()` + `zig_spline_chunk_stats()` + `zig_spline_fit_parallel()`. |
| `smooth_eval(fit, at)` | named list (by `family$params`) of n × `length(at)` matrices | Parameters on cycle positions `at`. ZIG: `p_zero = plogis(·)`, L-CV clamped to `lcv_range`, Hosking shape, `scale = l1 / shape`, status 1 → `p_zero = 1`. Today `zig_spline_params()`. |
| `smooth_params(evaluated, j)` | `params` matrix | Convenience: the n × p matrix for the `j`-th position of a `smooth_eval()` result (today's `cbind(p$shape[, 1], p$scale[, 1], p$p_zero[, 1])`). |
| `fit_status_label(status, family)` | character | Family status codes → spec §8 `fit_status`: ZIG 0 `ok`, 1 `no_snow`, 2 `p_zero_only`, 3 `no_data`, 4 `p_zero_not_converged`, 5 `constant`. |
| `smooth_tables(fit, subjects)` | `list(coef, diagnostics)` tibbles | Today's layout-v0 side tables (mco-snotel `R/normals.R:165-176`): `coef` = subject columns + `nu_01`…`nu_24`, `l1_01`…, `lcv_01`…; `diagnostics` = subject columns + `status`, `iter`, `lambda_nu`, `lambda_l1`, `lambda_lcv`, `edf_nu`, `edf_l1`, `edf_lcv`. Dropped in layout v1 (spec §13). |

`normals_smooth` is a list: `coef` (n × 3k for ZIG), `edf` (n × 3), `iter` (n), `status`
(integer n), `family` (the family object), `at` (fit positions), `k`, `lambda`.

### 3.6 Distributions and percentiles, matrix level (`R/dist.R`, `R/percentiles.R`)

A `normals_dist` is one day's distribution for n subjects: a family, its parameters, the
reference matrix (for the empirical family and for fallback), `min_valid` and `fallback`. It is
the one object the CDF, quantile, mean, p_zero, percentile and outlook functions take.

| Function | Returns | Contract / today |
|---|---|---|
| `normals_dist(family, params, ref, min_valid, fallback)` | `normals_dist` | No defaults. `params`: n × p, or `NULL` for `family_empirical()`. `ref`: n × y, or `NULL` when `fallback = "none"` and the family is not empirical. `min_valid`: integer (length 1 or n), or `NULL` when `ref` is `NULL`. `fallback`: `"empirical"` or `"none"` (ignored, and may be `NULL`, for the empirical family). Errors if rows disagree, if `fallback = "empirical"` without `ref`, or if params columns differ from `family$params`. |
| `dist_cdf(dist, x)` | double n, in [0, 1] | P(X ≤ x), with the zero-mass convention: x = 0 → p_zero / 2 (spec §11). NA where x is NA or the row has no normal. ZIG today `zig_p()`; empirical `ecdf_p()`. |
| `dist_quantile(dist, probs)` | n × m double | `probs` a vector (shared by all rows; colnames `quantile_columns(probs)`) or an n × m matrix of row-specific probabilities (the outlook). ZIG with `fallback = "empirical"` is today's `zig_q()`; with `"none"` it is `zig_smooth_q()` (including its p_zero ≥ 1 → 0 rule); empirical is `ecdf_q()`. p = 1 is `NA` for ZIG rows with a positive part (today's behavior, kept). |
| `dist_mean(dist)` | double n | ZIG `"empirical"` = today's `zig_mean()`; `"none"` = `zig_smooth_mean()`; empirical = `ecdf_mean()`. |
| `dist_p_zero(dist)` | double n | The chance of exactly zero: ZIG `p_zero` (NA without a normal); empirical: zeros / valid years where valid ≥ `min_valid`. Error for a family without a zero mass. |
| `dist_n_valid(dist)` | integer n | Non-missing years per row of `ref`. |
| `mask_sparse(params, ref, min_valid)` | `params` | Sets every parameter `NA` on rows with fewer than `min_valid` valid years. Today mco-snotel `R/normals.R:148-154` (smoothed only). |
| `as_ptile(p, digits)` | double | `round(p * 100, digits)`: today's `round(p * 100, 2)` (tables) and, wrapped in `as.integer()`, `round(p * 100)` (grids). `digits` has no default. |
| `standardize(p, method, as, z_max = 4)` | double or integer | Standard normal quantile of `p`, clamped to ±`z_max`. `method = "acklam"`: today's `qnorm_grid(p * 100)` (Acklam's approximation, ×100, rounded) — `as = "hundredths"` returns that integer (grids), `as = "z"` returns it / 100 (tables). `method = "qnorm"`: `stats::qnorm`, unrounded, `as = "z"` only (the outlook's z0). No defaults for `method` or `as`. |

### 3.7 Normals and percentiles, table level (`R/tables.R`)

```r
normals_table(ref, dates, subjects, family, year_type, min_valid, fallback, smoothed_mask,
              flavors = names(NORMALS_FLAVORS), probs = NORMALS_QUANTILES, cores = 1L)
```

- `ref`: `function(date)` returning the reference matrix for that date (rows in `subjects`
  order), or a list of matrices named by `format(date)`. Called once per date per use; the
  caller may memoise.
- `dates`: the normals year's days (`year_dates()`). Must be one whole year of `year_type`.
- `subjects`: a data frame of key columns, one row per `ref` row (`triplet`; or `zone`; or
  `zone`, `swe_elevation`, `area`). Bound to every output row as-is.
- `family`: the family for `fitted` and `smoothed`.
- `min_valid`: integer, or `function(ref)` returning one (see D3).
- `fallback`: named character, `c(fitted = …, smoothed = …)`, each `"empirical"` or `"none"`.
- `smoothed_mask`: `"none"` or `"min_valid"` (D2).
- The smoother is fitted at every date except Feb 29 and evaluated at every date (today's rule
  in both repos).

Returns `list(normals, smooth)`: `smooth` is the `normals_smooth` fit (NULL without `smoothed`),
and `normals` a tibble, one row per flavor × date × subject, ordered by flavor (in `flavors`
order), date, then subject (input order):

| Column | Type | Contents |
|---|---|---|
| `flavor` | chr | `empirical`, `fitted`, `smoothed` |
| `date` | Date | The normals date (Phase 1 pipelines key on it) |
| `day` | chr | `MM-DD` (spec §8) |
| subject columns | as given | |
| `distribution` | chr | `family$name` where the row's parameters are complete; `"empirical"` for the empirical flavor and for rows that fell back; `NA` where the row has no normal |
| `n_valid` | int | Reference years with a value |
| `ref_first_wy`, `ref_last_wy` | int | `reference_span()` (needs year colnames on `ref`) |
| `mean` | dbl | `dist_mean()` |
| family parameters | dbl | `shape`, `scale`, `p_zero` for ZIG; `NA` on empirical rows |
| `fit_status` | chr | `smoothed` rows: `fit_status_label()`; `NA` on other rows |
| `q0000` … `q1000` | dbl | `dist_quantile(probs)` |

`outlook_length` and `outlook_alpha` are not added here; `outlook_lengths()` returns them for the
pipeline to join (as mco-snotel `outlook` does in `R/stages.R`).

```r
percentiles_table(normals, ref, subjects, obs, family, min_valid, fallback,
                  ptile_digits, flavors = names(NORMALS_FLAVORS))
```

- `normals`: a `normals_table()$normals` (spec names; parameters are read from it).
- `ref`, `subjects`, `min_valid`: as above, for the same normals year.
- `obs`: data frame with the subject key columns, `date` and `value`. Rows whose subject is not
  in `subjects`, or whose date has no normals, are dropped (today's behavior).
- `fallback`: `c(fitted = …, smoothed = …)` for the CDF (D1).
- `ptile_digits`: no default (both pipelines' tables use 2).

Returns a tibble ordered by flavor, date, subject: `flavor`, `date`, subject columns, `value`,
`ptile` (`as_ptile(p, ptile_digits)`), `standardized` (`standardize(p, "acklam", "z")`), the
spec §8 `percentiles` schema.

### 3.8 Outlook (`R/outlook.R`)

The outlook of mco-snotel `outlook:R/outlook.R`, made independent of stations and of the
normals table: it takes `normals_dist` objects, so it serves zones and other families. The α and
β tables are arguments (spec §11), so a pipeline always names the table it uses. The calibrated
tables ship in `inst/extdata/` (`outlook_alpha.csv`, `outlook_spread.csv`, with `PROVENANCE`);
pipelines read them with `system.file("extdata", …, package = "normals")`.

```r
outlook_control(k_max = 3, fit_lags = c(3, 7, 14, 21, 30, 45, 60, 90, 120),
                pool_days = 15L, min_pairs = 8L, min_lags = 3L, z_max = 4,
                l_grid = exp(seq(log(2), log(3000), length.out = 400)),
                p_zero_floor = 0.5, anchor_range = c(0.1, 10), rho_floor = 1e-12,
                quantile_digits = 1, p_zero_digits = 3, length_digits = 2, alpha_digits = 3)
```
Today's constants (`OUTLOOK_K_MAX`, `OUTLOOK_FIT_LAGS`, `OUTLOOK_POOL_DAYS`,
`OUTLOOK_MIN_PAIRS`, `OUTLOOK_Z_MAX`, `OUTLOOK_L_GRID`, the 0.5-year floor, the anchor clamp, the
`1e-12` in `sqrt(pmax(1e-12, 1 - rho^2))`, and the roundings at `outlook.R:238, 280, 285`).

| Function | Returns | Contract / today |
|---|---|---|
| `outlook_alpha_at(table, dates, year_type)` | double, `length(dates)` | `table`: data frame with `day` (cycle day 1–365) and `alpha`, already filtered to one period and flavor by the caller. Linear interpolation on `cycle_day()`, constant beyond the ends. Today `outlook_alpha()`. |
| `outlook_persistence(lead, length, alpha)` | double | ρ = exp(−(lead / length)^alpha). |
| `outlook_spread(z0, rho, beta, control = outlook_control())` | double | k = min(exp(−β ρ min(z0, 0)), `k_max`). Today `outlook_spread_k()`. |
| `outlook_start(issue, x, control = outlook_control())` | `list(z0, anchor, floor)` | `issue`: the issue date's `normals_dist`; `x`: observed values (n). `floor` = `x > 0` (where the half-year p_zero floor applies, zero-mass families only); `z0` = clamped `standardize(dist_cdf(floored issue, x), "qnorm", "z")`; `anchor` = today's `outlook_anchor()`. |
| `outlook_quantiles(target, start, rho, beta, probs = NORMALS_QUANTILES, control = outlook_control())` | n × m double, unrounded | The target date's quantiles at u = Φ(ρ z0 + Φ⁻¹(clamped p) σ k), σ = √max(`rho_floor`, 1 − ρ²), times `anchor^rho`. |
| `outlook_p_zero(target, start, rho, beta, control = outlook_control())` | double n | Chance of exactly zero on the target date (spec §8 `p_zero`; today `p_snow_free`). Zero-mass families only. |
| `outlook(issue, targets, x, issue_date, target_dates, subjects, length, alpha, beta, probs = NORMALS_QUANTILES, control = outlook_control())` | tibble | Composes the above for one issue date and every target date (`targets`: list of `normals_dist`, one per `target_dates`, all with `issue`'s rows). `length`, `alpha`: per subject (or scalar) at the issue date. Columns: `date`, subject columns, `q0000` … `q1000` (rounded to `quantile_digits`), `p_zero` (rounded to `p_zero_digits`). Rows without z0 or length are dropped; sorted by subject then date. `NULL` when `target_dates` is empty. Today `calc_station_outlook()`. |
| `outlook_lengths(dists, dates, subjects, alpha, control = outlook_control())` | tibble | `dists`: a `normals_dist` per normals date (each with its `ref`), `dates` those dates, `alpha` one per date (`outlook_alpha_at()`). Correlates each reference year's clamped z between dates `fit_lags` apart, pools ± `pool_days`, fits ρ(h) by least squares over `l_grid`, fills gaps with the subject's median. Columns `date`, subject columns, `outlook_length` (rounded to `length_digits`), `outlook_alpha` (`alpha_digits`). Today `outlook_correlations()` + `outlook_fit_lengths()` + `outlook_lengths()`. |
| `outlook_crps(q, y, probs)` | double n | Twice the mean pinball loss of quantiles `q` (n × m at `probs`) against outcomes `y`. Today `crps_q()` in `scripts/outlook-hindcast.R`. |
| `outlook_pit(target, start, rho, beta, y, u, control = outlook_control())` | double n | Probability integral transform of outcome `y` on the anchored scale; a zero outcome is placed at `u ×` the chance of zero (`u` supplied by the caller, so the package draws no random numbers). Today inline in `scripts/outlook-hindcast.R`. |

The hindcast driver and the α / β calibration (mco-snotel `outlook:scripts/outlook-hindcast.R`,
`outlook-calibrate.R`) move to `data-raw/` as scripts built on these functions. They take a
dataset's record and write the CSVs in `inst/extdata/`.

### 3.9 Plot helpers (`R/plot.R`)

The shared frame, palette, USDM classes, legends and logos, from mco-snodas `R/plot.R:1-278`
(everything before `huc_display_name`, which is the part mco-snotel vendors). Per-dataset plots
(`plot_snodas()`, `plot_hypsomeswe()`, `plot_snotel()`) stay in the pipelines. ggplot2 and the
other plotting packages are in Suggests (Q4); each helper checks for them and stops with the
package to install.

| Function / object | Today |
|---|---|
| `plot_fix_frame(g, axis_width = grid::unit(0.45, "in"), xlab_height = grid::unit(0.2, "in"))` | `fix_plot_frame()` |
| `plot_bands(usdm)` | `snodas_bands()` (no default for `usdm`) |
| `plot_month_breaks(limits)`, `plot_month_labels(x)` | `snodas_month_breaks()`, `snodas_month_labels()` |
| `plot_legend_theme()` | `snodas_legend_theme()` |
| `plot_logo(which)` | `snodas_logo()` (`"mco"`, resampled to 400 px tall) and `snodas_umt_logo()` (`"umt"`), read once per session from `system.file("plot", "mco_logo.png" / "umt_logo.png", package = "normals")` |
| `plot_logo_grob(img, height, ...)` | `snodas_logo_grob()`; class `normals_logo`, whose `makeContent()` method is registered (`S3method(grid::makeContent, normals_logo)`) |
| `plot_current_legend(obs_lwd, likely_range = FALSE)` | `snodas_current_legend()`; with `likely_range = TRUE`, mco-snotel `snotel_outlook_legend()` |
| `plot_legend_row(g)` | `snodas_legend_row()` |
| `plot_palette()` | `snodas_palette()`, from `system.file("plot", "plot_palette.png", package = "normals")` |
| `plot_quantize(nr)` | `snodas_quantize()` |
| `plot_save_png(p, file, width = 8, height = 16 / 3, dpi = 187.5)` | `snodas_save_png()` |
| `plot_conditions(value, median, ptile)` | `snodas_conditions()` |
| `register_plot_fonts()` | `register_snodas_fonts()` (`R/utils.R:41-72`) |
| `PLOT_TITLE_SIZE` (11.5), `PLOT_UMT_LOGO_HEIGHT` (0.3) | `SNODAS_TITLE_SIZE`, `UMT_LOGO_HEIGHT` |

`snodas_resample()` becomes internal. `snodas_ref_years()` is not moved: it encodes the SNODAS
period rule, which belongs to the pipeline (or, after Phase 2, the resolver).
`scripts/build-plot-palette.R` moves to `data-raw/` and writes `inst/plot/plot_palette.png`.

## 4. The family plug-in interface

A family is an S3 object of class `normals_family`: a named list. `new_family()` builds one and
`validate_family()` checks it. Adding a family means writing these slots, in any package or
script; the core (`normals_dist()`, `dist_*()`, `normals_table()`, `percentiles_table()`,
`outlook*()`) dispatches through them and never tests a family's name.

| Slot | Type | Contract |
|---|---|---|
| `name` | string | The `distribution` value (spec §8): `"zero_inflated_gamma"`. |
| `params` | character | Parameter names, in column order: `c("shape", "scale", "p_zero")`. Empty for `empirical`. |
| `zero_mass` | logical | TRUE if P(X = 0) is part of the family. Turns on the zero conventions: x = 0 → p_zero / 2 in the CDF, p ≤ p_zero → 0 in the quantile function, `dist_p_zero()`, and the outlook's half-year floor and `p_zero` column. |
| `stats_names` | character | Column names of `stats()`'s output: `c("n", "z", "l1", "lcv")`. |
| `stats` | `function(ref, control)` → n × s matrix | One day's statistics from a reference matrix. The fit and the smoother read only these, so the gridded build can store them per day and fit later (today's `suffstats` row blocks). |
| `fit` | `function(stats, min_valid, control)` → n × p | One day's parameters. Rows with `stats[, "n"] < min_valid` are all `NA`. Rows the family cannot fit have `NA` in the non-zero-mass parameters (the fallback signal). |
| `smooth` | list or `NULL` | `NULL`: the family has no `smoothed` flavor. Otherwise: |
| `smooth$terms` | data frame | Documentation of the smoothing transform: one row per smoothed quantity with `term`, `scale` (`"logit"`, `"log"`), `method` (`"binomial"`, `"wls"`), `weight`. ZIG: `p_zero` on logit by penalised binomial on (z, n); `l1` and `lcv` on log by penalised WLS weighted by positive years. Written into provenance by the pipeline. |
| `smooth$fit` | `function(stats, basis, cores, control)` → list(coef, edf, iter, status) | Fits all rows. `stats` stat-major; `basis` from `smooth_basis()`. Must be row-independent (so chunking cannot change results). |
| `smooth$eval` | `function(coef, X, status, control)` → named list of n × nrow(X) matrices | Back-transform to `params` at basis rows `X`. |
| `smooth$status` | named integer | Status codes and their spec §8 labels. |
| `cdf` | `function(params, x, ref, min_valid, fallback, control)` → double n | P(X ≤ x). Fallback rows (see `fit`) use the ECDF of `ref` when `fallback = "empirical"`, `NA` when `"none"`. |
| `quantile` | `function(params, probs, ref, min_valid, fallback, control)` → n × m | `probs` a vector or an n × m matrix. Same fallback rule. |
| `mean` | `function(params, ref, min_valid, fallback, control)` → double n | Same fallback rule. |
| `p_zero` | `function(params, ref, min_valid, control)` → double n, or `NULL` | Required iff `zero_mass`. |
| `control` | list | The family's constants (`zig_control()`). |

Notes:
- **Fallback is the family's job,** not the core's, because it mixes the family and the ECDF.
  ZIG's fallback CDF puts x = 0 at the *fitted or smoothed* p_zero / 2 but x > 0 on the
  *empirical* zero fraction (`zig.cpp:199-212`); its fallback quantile tests p ≤ the smoothed
  p_zero before reading the ECDF (`zig.cpp:259-270`). A generic "use the empirical family"
  would change smoothed output.
- **`family_empirical()`** fills the same slots with `params = character()`, `fit = NULL`,
  `smooth = NULL`, `zero_mass = TRUE` (its zero fraction is the empirical p_zero).
- **Smoothing engines.** ZIG's `smooth$fit` wraps today's single C++ routine (`zig_spline_fit`),
  so its output is unchanged. A future family supplies its own `smooth$fit`; whether the core
  should also export generic engines (penalised binomial and WLS curves on the shared basis) for
  such families is Q9.
- **Example, a normal family for temperature (gridMET).** `params = c("mean", "sd")`,
  `zero_mass = FALSE`, `stats` = n, mean, log sd; `fit` = moments; `smooth$terms` = mean on the
  identity scale and sd on log, both WLS weighted by n; `cdf` / `quantile` = `pnorm` / `qnorm`;
  `p_zero = NULL`. Nothing in the core changes to add it.

## 5. Differences between the pipelines, and the arguments that encode them

"Grid" is mco-snodas gridded (pixels), "zonal" is mco-snodas zonal (HUCs and bands), "SNOTEL" is
mco-snotel. Each row gives today's behavior, where it lives, and the 2.0 argument that each
pipeline sets to reproduce it. Rows D1–D7 are method differences the package must express;
rows P1–P7 are differences that stay in pipeline code (inputs or storage) but that the golden
fixtures must reproduce.

| # | Difference | mco-snodas | mco-snotel | 2.0 argument (Phase 1 values) |
|---|---|---|---|---|
| D1a | **Smoothed quantiles and mean: empirical fallback.** Rows whose smoothed `shape`/`scale` are NA (fit status 2, or 1) | Grid: **none** — `zig_smooth_q()` / `zig_smooth_mean()` (`R/calc_snodas_zig_spline_normals.R:135-151`, called at `R/normals_annual.R:473-476`), NA above p_zero, p_zero ≥ 1 → 0. Zonal: **empirical** — `zig_q(params, …, h)` / `zig_mean(params, h)` (`R/normals_annual.R:312, 317`) | **Empirical** — `zig_q(params, q_probs, h)` / `zig_mean(params, h)` (`R/normals.R:155, 157`) | `normals_dist(fallback = …)` for the quantile/mean dist; `normals_table(fallback = c(…, smoothed = …))`. Grid `"none"`, zonal `"empirical"`, SNOTEL `"empirical"` |
| D1b | **Smoothed percentiles (CDF): empirical fallback** | Grid and zonal: **empirical** — `zig_p(smooth params, swe, hist_mat)` (`R/calc_snodas_normals.R:160-162`; `R/calc_zonal_normals.R:322-331`). So the grid's smoothed flavor falls back in its percentiles but not in its quantiles (D1a). | **Empirical** — `zig_p(params_of(nf), o$wteq, h)` (`R/normals.R:204`) | The percentile dist's `fallback`; `percentiles_table(fallback = c(…, smoothed = "empirical"))`. All `"empirical"` |
| D1c | **Fitted flavor: empirical fallback** (rows with < 4 positive values) | **Empirical** everywhere: `zig_q`/`zig_mean`/`zig_p` with the record (`R/calc_snodas_normals.R:108, 117, 157`; `R/calc_zonal_normals.R:248-258, 318-320`) | **Empirical** (`R/normals.R:122, 126, 204`) | `fallback = c(fitted = "empirical", …)`. Same in all; explicit for symmetry, no default |
| D2 | **Sparse-record masking of the smoothed flavor.** The spline fit has no minimum record | **None** (grid `R/normals_annual.R:259-266` only blanks status 3, which is already all NA; zonal `R/normals_annual.R:310-311`) | **Masked**: parameters set NA where valid years < `min_valid` (`R/normals.R:144-154`), before quantiles, mean and (via the stored table) percentiles | `normals_table(smoothed_mask = …)` / `mask_sparse()`. mco-snodas `"none"`, SNOTEL `"min_valid"` |
| D3 | **`min_valid`: the threshold and what "years in the period" counts.** The rule is min(10, n) in both, but n differs | n = `ncol(hist_mat)`, the reference **layers found**: missing COGs are dropped before stacking (`R/calc_snodas_normals.R:56-60`; `R/calc_zonal_normals.R:97-101`), so a missing year lowers n. Hard-coded in C++ (`src/ecdf_impl.h:14-18`, used at `ecdf_impl.h:37`, `ecdf.cpp:80`, `zig.cpp:75, 314`). Not applied to the smoothed flavor (D2) | n = `length(years)`, the **period's years within the archive**, all-NA columns kept (`R/normals.R:43-47, 76-81`); the same C++ rule, plus `min(10L, length(years))` in R for the mask (`R/normals.R:148`) | `min_valid` on `fit_day()`, `normals_dist()`, `mask_sparse()`, `normals_table()`, `percentiles_table()` — required. mco-snodas: `function(ref) min_valid_years(ncol(ref), cap = 10)`. SNOTEL: `min_valid_years(length(years), cap = 10)` (equal to `ncol(ref)` there). The C++ takes `min_valid` instead of deriving it |
| D4 | **Percentile rounding** | Grid: integer percent, `as.integer(round(p * 100))` (`R/calc_snodas_normals.R:173, 179, 200`). Zonal: `round(p * 100, 2)` (`R/calc_zonal_normals.R:306`) | `round(p * 100, 2)` (`R/normals.R:206`) | `as_ptile(p, digits)` / `percentiles_table(ptile_digits = …)`. Grid `0` (then `as.integer`), zonal and SNOTEL `2` |
| D5 | **Standardized anomaly encoding** | Grid: `qnorm_grid(p * 100)`, integer hundredths (`R/calc_snodas_normals.R:186, 192, 205`). Zonal: `qnorm_grid(ptile) / 100` (`R/calc_zonal_normals.R:307`) | `qnorm_grid(p * 100) / 100` (`R/normals.R:206`) | `standardize(p, method = "acklam", as = …)`. Grid `"hundredths"`, zonal and SNOTEL `"z"` |
| D6 | **z for the outlook vs published standardized** | — | Published: Acklam, rounded to 0.01 (D5). Outlook z0 and the persistence correlations: `stats::qnorm`, unrounded (`outlook:R/outlook.R:183, 272`) | `standardize(method = "qnorm")` inside `outlook_start()` / `outlook_lengths()`; `"acklam"` for published values. Fixed in the outlook functions (they are new, not published yet) |
| D7 | **Empirical p_zero's minimum record in the outlook** | — | ECDF flavor only (hindcast option, not in production): fixed `n >= 10L` (`outlook:R/outlook.R:110`), not min(10, n) | `dist_p_zero()` uses the dist's `min_valid`. Differs from today only for the ECDF flavor with < 10 period years, which no published output uses; recorded, no argument |
| P1 | **Reference dates** | `d − 0:30` calendar years, kept if ≥ 2003-09-30 and before the date's water year (`R/pipeline_sync.R:231-234`, `R/normals_annual.R:124-128`). For **Sep 30** this keeps **2003-09-30 (WY2003)**: 24 reference years, not 23 | `period_years()` then `reference_dates()` (`R/normals.R:34-53`); the `snodas` period starts at WY2004 for every day | Pipeline: mco-snodas calls `reference_dates(d, (W-30):(W-1), "wy")` and keeps dates ≥ 2003-09-30, as today. See Q2 |
| P2 | **Period of record (`por`)** | — | First year of the **whole archive** for every station (`R/normals.R:44`), not each station's first year (spec §5). Each station's earlier years are NA columns | Pipeline (or resolver) passes `years`; `min_valid` follows D3. See Q1 |
| P3 | **Float32 parameters on the grid** | Grid `zig` percentiles (nightly) use the frozen `{d}_zig_params.tif`, written FLT4S (`R/calc_snodas_normals.R:113-115, 145`). Grid smoothed parameters go through FLT4S tiles (`R/normals_annual.R:80-86, 263-267`) before `zig_smooth_q`/`zig_smooth_mean` and the smoothed percentiles. Zonal and the grid's annual `zig` quantiles use doubles | Doubles from the parquet table | Pipeline: passes the float32-rounded matrix it reads, as today. The package takes `params` as given |
| P4 | **INT2S grids** | Grid means and quantiles are written as INT2S COGs (`R/snodas_raw.R:151`, `R/normals_annual.R:473-476`) | — | Pipeline storage; the package returns doubles |
| P5 | **Frozen record for percentiles** | Zonal: the saved raw means (`R/calc_zonal_normals.R:188-196`). Grid: the COGs, re-read each night | The frozen `normals/reference/wy{W}.parquet` (`R/stages.R:103-106, 165-170`) | Pipeline: supplies `ref` |
| P6 | **NaN for an all-missing zone mean** | Zone and band means are `sum(x·a, na.rm) / sum(a·!is.na(x))`, so 0/0 = **NaN**, not NA (`R/calc_zonal_normals.R:139-141, 160-161`) | — | The package treats NaN as missing (`ISNAN`), as today. Fixtures compare with `identical()`, which tells NaN from NA |
| P7 | **Output columns and order** | Zonal per-date tables: `date`, `zone` (+ `swe_elevation`, `area`), [`swe`, `swe_ptile`, `swe_standardized`], `mean`, `zig_*`, `q*` | `flavor`, `date`, `triplet`, `n_valid`, `ref_first_wy`, `ref_last_wy`, `mean`, `zig_*`, `q*`, sorted flavor, date, triplet | Pipeline: selects and renames from `normals_table()` (§6.4) |

Constants that are the **same** in both pipelines (so defaulted in `zig_control()` /
`outlook_control()`): L-moment minimum 4 positive values; k = 24; λ = (100, 1000, 1000); L-CV
clamp [0.001, 0.999]; `min_pos = 10`, `min_days = 30`, `maxit = 100`, `tol = 1e-6`; Feb 29 left out
of the spline fit and evaluated at the half day; the 19 probabilities; ECDF type 7 with
x = 0 → p_zero / 2; the ±4 clamp.

## 6. Mapping from today's code to 2.0 calls

`zig <- family_zig()`, `emp <- family_empirical()`, `P <- NORMALS_QUANTILES`. "mv" is the
pipeline's `min_valid` (D3).

### 6.1 Shared (vendored) functions

| Today (both repos) | 2.0 |
|---|---|
| `SNODAS_QUANTILES` | `NORMALS_QUANTILES` (a pipeline may keep `SNODAS_QUANTILES <- normals::NORMALS_QUANTILES` in Phase 1) |
| `water_year(d)` | `year_of(d, "wy")` |
| `water_year_start(d)` | `year_start(year_of(d, "wy"), "wy")` |
| `water_year_dates(d)` | `dates_through(d, "wy")` |
| `normals_dates(wy)` | `year_dates(wy, "wy")` |
| `normals_cycle_day(d)` | `cycle_day(d, "wy")` |
| `ecdf_p(raw, obs)` | `dist_cdf(normals_dist(emp, NULL, raw, mv, NULL), obs)` |
| `ecdf_q(raw, probs)` | `dist_quantile(normals_dist(emp, NULL, raw, mv, NULL), probs)` |
| `ecdf_mean(raw)` | `dist_mean(normals_dist(emp, NULL, raw, mv, NULL))` |
| `ecdf_stats(raw)` | `list(mean = dist_mean(d), quantiles = dist_quantile(d, P))` with `d` the empirical dist |
| `zig_lmom(mat)` / `zig_spline_suffstats(mat)` | `day_stats(mat, zig)` |
| `zig_params_from_lmom(lm, n_time)` | `fit_day(lm, zig, min_valid_years(n_time, 10))` |
| `zig_fit(mat)` | `fit_day(day_stats(mat, zig), zig, mv)` |
| `zig_shape_from_lcv(lcv)` | internal (`smooth$eval` and `fit` use it) |
| `zig_p(params, obs, raw)` | `dist_cdf(normals_dist(zig, params, raw, mv, "empirical"), obs)` |
| `zig_q(params, probs, raw)` | `dist_quantile(normals_dist(zig, params, raw, mv, "empirical"), probs)` |
| `zig_mean(params, raw)` | `dist_mean(normals_dist(zig, params, raw, mv, "empirical"))` |
| `zig_smooth_q(params, probs)` | `dist_quantile(normals_dist(zig, params, NULL, NULL, "none"), probs)` |
| `zig_smooth_mean(params)` | `dist_mean(normals_dist(zig, params, NULL, NULL, "none"))` |
| `zig_spline_basis(k, at)` | `smooth_basis(at, k)` |
| `zig_spline_fit_stats(stats, basis, …)` | `smooth_fit(stats, at, zig)` (basis built inside from `at` and `k`) |
| `zig_spline_chunk_stats()` + `zig_spline_fit_parallel(chunks, basis, cores)` | `smooth_fit(stats, at, zig, cores)` |
| `zig_spline_params(coef, X, status)` | `smooth_eval(fit, at)`; one date: `smooth_params(smooth_eval(fit, at[j]), 1)` |
| `qnorm_grid(p * 100)` | `standardize(p, "acklam", "hundredths")`; `/ 100` → `as = "z"` |
| `register_snodas_fonts()` | `register_plot_fonts()` |
| `ZIG_SMOOTH_K`, `ZIG_SMOOTH_LAMBDA`, `ZIG_SMOOTH_LCV_RANGE` | `zig_control()$k`, `$lambda`, `$lcv_range` |
| `write_provenance()` (snodas `R/utils.R:115-124`) | stays in the pipeline (it writes a file) |
| plot helpers (`snodas_*`, `fix_plot_frame`) | `plot_*` (§3.9) |
| `MIN_VALID_YEARS` (C++) | `min_valid_years(n, cap = 10)`, passed explicitly |

### 6.2 mco-snodas (1B)

| Today | 2.0 |
|---|---|
| `calc_snodas_normals()` ECDF block (`:92-101`) | `d_e <- normals_dist(emp, NULL, hist_mat, mv, NULL)`; `dist_mean(d_e)`, `dist_quantile(d_e, P)` |
| ZIG block (`:106-117`) | `lm <- day_stats(hist_mat, zig)`; `pz <- fit_day(lm, zig, mv)`; `d_z <- normals_dist(zig, pz, hist_mat, mv, "empirical")`; `dist_quantile(d_z, P)`, `dist_mean(d_z)`; `lm` is also the `suffstats` block written per row block |
| nightly refit (`:148`) | `fit_day(day_stats(hist_mat, zig), zig, mv)` |
| percentiles (`:156-162`) | `dist_cdf()` on the ECDF dist, on `normals_dist(zig, frozen_f32_params, hist_mat, mv, "empirical")`, and on the smoothed one built the same way |
| ptile / standardized COGs (`:173-209`) | `as.integer(as_ptile(p, 0))`; `standardize(p, "acklam", "hundredths")` |
| `calc_zonal_normals()` normals (`:223-276`) | `normals_table(ref = <hist columns of the raw means, per date>, dates = d, subjects = zone keys, family = zig, year_type = "wy", min_valid = function(r) min_valid_years(ncol(r), 10), fallback = c(fitted = "empirical", smoothed = "empirical"), smoothed_mask = "none", flavors = c("empirical", "fitted"))` for one date, then per-flavor legacy tables (§6.4) |
| `calc_zonal_normals()` percentiles (`:293-332`) | `dist_cdf()` per flavor with `fallback = "empirical"`, `as_ptile(p, 2)`, `standardize(p, "acklam", "z")`; or `percentiles_table(…, ptile_digits = 2)` |
| `annual_fit_block()` (`R/normals_annual.R:216-271`) | `fit <- smooth_fit(chunked stats, cycle_day(fit_dates, "wy"), zig, cores)`; `smooth_tables()`-style coef/diagnostic bands; per date `smooth_params(smooth_eval(fit, cycle_day(dates[j], "wy")), 1)` |
| `annual_fit_zonal()` (`:276-334`) | `smooth_fit(stack_stats(lapply(fit_dates, \(d) day_stats(hist(d), zig))), cycle_day(fit_dates, "wy"), zig, cores)`; per date `normals_dist(zig, params, h, mv, "empirical")` → `dist_mean`, `dist_quantile` (D1a zonal). Equivalently `normals_table(…, flavors = "smoothed", fallback = c(fitted = "empirical", smoothed = "empirical"), smoothed_mask = "none")` over the year |
| `annual_materialize()` (`:457-487`) | `terra::app(params_r, \(v) dist_mean(normals_dist(zig, v, NULL, NULL, "none")))`, likewise `dist_quantile(…, P)` (D1a grid) |
| `normals_reference_years(d)` | stays; may become `reference_dates(d, (W-30):(W-1), "wy")` filtered to ≥ `SNODAS_START` (P1) |
| `snodas_ref_years()` | stays in the pipeline |

### 6.3 mco-snotel (1C)

| Today | 2.0 |
|---|---|
| `swe_matrix(x)` | `record_matrix(x, "triplet", "date", "wteq")` |
| `swe_at()` / `reference_record()` / `record_span()` | `record_at()` / `reference_record(rec, d, years, "wy")` / `reference_span()` |
| `period_years(wy, period, first_wy)` | `period_years("rolling", wy, length = 30)` (30yr); `period_years("rolling", wy, length = 30, earliest = 2004)` (snodas); `period_years("por", wy, first_year = first_wy)` (por) — or kept in the pipeline (Q3) |
| `normals_years(sm, wy, period)` | stays in the pipeline (it clips to the archive's first year, P2) |
| `calc_station_normals(sm, wy, period, cores)` | `res <- normals_table(ref = \(d) reference_record(rec, d, years, "wy"), dates = year_dates(wy, "wy"), subjects = tibble(triplet = rec$subjects), family = zig, year_type = "wy", min_valid = min_valid_years(length(years), 10), fallback = c(fitted = "empirical", smoothed = "empirical"), smoothed_mask = "min_valid", cores = cores)`; `smooth_tables(res$smooth, subjects)` for `zig_smooth_coef` / `_diagnostics` |
| `calc_station_percentiles(sm, obs, normals, years)` | `percentiles_table(normals, ref, subjects, obs = rename(obs, value = wteq), family = zig, min_valid = mv, fallback = c(fitted = "empirical", smoothed = "empirical"), ptile_digits = 2)` |
| `outlook_alpha(period, dates)` | `outlook_alpha_at(alpha_tbl[alpha_tbl$period == period & alpha_tbl$flavor == "zig_smooth", ], dates, "wy")` (`alpha_tbl` read from the package's `inst/extdata/outlook_alpha.csv`) |
| `outlook_spread_beta(period)` | stays in the pipeline (reads `config/outlook_spread.csv`); passes `beta` |
| `outlook_spread_k()`, `outlook_params()`, `outlook_cdf()`, `outlook_p_zero()`, `outlook_quantiles()`, `outlook_anchor()` | `outlook_spread()`, `outlook_start()`, `dist_cdf()`, `outlook_p_zero()`, `outlook_quantiles()` (internal pieces otherwise) |
| `outlook_correlations()` + `outlook_fit_lengths()` + `outlook_lengths(sm, normals, years, period)` | `outlook_lengths(dists = <smoothed normals_dist per date>, dates, subjects, alpha = outlook_alpha_at(…))`, then `flavor = "zig_smooth"` added at the edge |
| `calc_station_outlook(sm, normals, years, obs, issue, beta)` | `outlook(issue = <dist at issue>, targets = <dists after issue>, x = obs$wteq, issue_date, target_dates, subjects, length = <outlook_length at issue>, alpha = <outlook_alpha at issue>, beta)`, then rename `p_zero` → `p_snow_free` (Phase 1 edge) |
| `scripts/test_outlook.R`, `scripts/test_distributions.R` | `tests/testthat/test-outlook.R`, `test-ecdf.R`, `test-zig.R`, `test-smooth.R` (1A-v) |
| `scripts/outlook-hindcast.R`, `scripts/outlook-calibrate.R` | `data-raw/outlook-hindcast.R`, `data-raw/outlook-calibrate.R` (1A-iii) |
| `scripts/check-vendored.sh`, `src/`, `R/zig_spline_normals.R`, vendored `R/utils.R`, `R/plot.R`, logos, palette | removed in 1C |

### 6.4 Legacy names at the pipeline edge (Phase 1)

| 2.0 (spec) | Phase 1 published |
|---|---|
| `flavor` `empirical` / `fitted` / `smoothed` | `ecdf` / `zig` / `zig_smooth` — `flavor_legacy(x, "legacy")` |
| `shape`, `scale`, `p_zero` | `zig_shape`, `zig_scale`, `zig_p_zero` |
| percentiles `value`, `ptile`, `standardized` | `swe`, `swe_ptile`, `swe_standardized` |
| outlook `p_zero` | `p_snow_free` |
| `day`, `distribution`, `fit_status` | dropped (not in today's tables) |
| `ref_first_wy`, `ref_last_wy`, `n_valid` | kept by mco-snotel, dropped by mco-snodas |

## 7. Golden-fixture test plan

**Purpose:** prove that the 2.0 calls of §6 reproduce the baseline bit for bit, in the package's
CI and in each pipeline's CI (spec §11 "Parity").

**Baseline.** The Phase 0 control image (production `mco-snodas:latest`
`sha256:ab56aa88…`, rebuilt with pinned dependencies, plan 0C/0F) for mco-snodas, and
mco-snotel `main` at the commit Phase 0 records. The fixtures are produced once, by a script run
inside each baseline image, and committed with that image digest and commit in a
`fixtures/MANIFEST.json`. They are never regenerated by the 2.0 code.

**Bit-exactness is defined on the pinned platform** (decision of 2026-11-01, after 1A-v). The
fixtures are made, and CI checks them, inside `rocker/r-ver:4.5.3` by digest
(`sha256:f5a72ffe…`, the mco-snotel base image) with CRAN from its P3M snapshot (2026-04-23) and
the image's OpenBLAS 0.3.26 single-threaded with its kernel pinned (`OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, `OPENBLAS_CORETYPE=Nehalem`; `data-raw/fixtures.Dockerfile`). On that platform
every golden comparison is `identical()`. OpenBLAS otherwise picks its kernels by CPU, and the
penalised spline solve then differs between machines in the 7th to 11th significant digit, so off
that platform (a different R, BLAS/LAPACK or OpenBLAS setting; `tests/testthat/fixtures/platform.rds`)
the smoothed comparisons skip with a message, and the empirical, fitted and percentile comparisons
still run (`identical()` first, else within 1e-12 relative: libm differs in the last bits). **The G1
gate is not these fixtures:** it is shadow runs of each pipeline's own image (Phase 0's control
image for mco-snodas, the mco-snotel image), comparing 2.0 outputs with the baseline's in that same
image.

**Inputs and expected outputs** (all `.rds`, `version = 3`, xz):

| Set | Subjects | Inputs | Expected outputs |
|---|---|---|---|
| Stations | 5 SNOTEL stations, chosen for coverage: (1) a long full record; (2) a station with fewer than 10 years in `30yr` (D2 mask active); (3) a low station where some days have < 4 positive years (D1 fallback, status 2 or 5); (4) a station with missing years inside the period (NA columns); (5) a station whose record starts after WY2004 (`por` vs `snodas`, P2) | The frozen reference record rows for those stations (WY2027 build) and WY2027 observations through a fixed date | For each period (`30yr`, `snodas`, `por`): `calc_station_normals()$normals`, `$coef`, `$diagnostics`; `calc_station_percentiles()`; `outlook_lengths()`; `calc_station_outlook()` for 3 issue dates (Dec 1, Feb 1, Apr 15) |
| HUCs | 5 zones: a HUC2 and a HUC8 (both levels of aggregation); a HUC with no snow on record (status 1); a HUC with sparse snow (fallback rows); and one HUC8's `hypsomeswe` bands | Raw zone and band means for every date of WY2027 (the `raw/{d}` tables, 366 dates, those rows only) | ECDF, ZIG and ZIG-smooth per-date tables (`swe` and `hypsomeswe`) for the 8 G1 dates (10-01, 11-15, 01-01, 02-29, 04-01, 05-15, 07-01, 09-30) and the smoothed fit (`coef`, `status`) |
| Pixels | One 50 × 50 block: mixed mountain and plain, including pixels that are NA in early years (the SNODAS mask grew) and pixels with sparse snow | SWE (integer mm) for the 50 × 50 window on every reference date of the 8 G1 dates, plus the 365-day `suffstats` (n, z, l1, lcv) for the block | For the 8 dates: ECDF and ZIG mean and quantiles, ZIG params; the smoothed coefficients, status and per-date params (float32, as the tiles hold them); smoothed mean and quantiles from those params (no fallback); ptile and standardized grids for one observed date per flavor |

The pixel fixture is the largest: about 2,500 × 24 × 8 integers for the dates plus
2,500 × 365 × 4 doubles for the statistics (~29 MB uncompressed). Q5 asks where it lives.

**Tests (1A-v):**

1. **Parity** (`tests/testthat/test-golden-*.R`): each expected output is recomputed with the
   §6 2.0 calls and compared with `identical()` (which distinguishes NA from NaN, P6, and
   integer from double) on the pinned platform above. Off it, results that pass through the
   spline's linear algebra (smoothed coefficients, EDF, and what follows from them) are
   skipped, and the others are compared within 1e-12 relative after `identical()` fails. No
   tolerance is widened to pass (plan rule). In the pipeline images (G1) the target is
   `identical()` throughout.
2. **Ported unit tests:** `scripts/test_distributions.R` (Parts 1–4) and
   `scripts/test_outlook.R` (cases 1–7) as testthat files, with their tolerances unchanged.
3. **Interface tests:** every required argument missing → error; shape mismatches → the
   documented error; a toy plug-in family (normal) passes `validate_family()` and runs through
   `normals_table()` and `percentiles_table()` without touching the core; `cores = 1` and
   `cores = 2` give identical smoothed fits.
4. **Calendar tests:** `cycle_day()` for both year types including Feb 29; `reference_dates()` on
   Feb 29 and on Oct–Dec dates; `period_years()` against every role example in spec §5 (WY2035
   `snodas` = `30yr`; `noaa` absent before WY2042 for mco-snodas's configuration).
5. **Pipeline parity in CI:** each pipeline runs the same fixtures through its own call sites
   (1B, 1C), so a pipeline cannot drift from the package it pins.

## 8. Package layout, dependencies and CI

```
DESCRIPTION  NAMESPACE  LICENSE  LICENSE.md  README.md
R/          constants.R calendar.R record.R family.R fit.R smooth.R dist.R
            percentiles.R tables.R outlook.R plot.R normals-package.R
src/        (1A-i) ecdf_impl.h ecdf.cpp zig.cpp zig_spline.cpp qnorm_grid.cpp RcppExports.cpp
inst/plot/  (1A-iv) mco_logo.png umt_logo.png plot_palette.png
data-raw/   (1A-iii, 1A-iv) outlook-hindcast.R outlook-calibrate.R build-plot-palette.R
tests/testthat/
docs/API.md (this file; build-ignored)
```

**Dependencies.** Imports, all CRAN or base: `mgcv` (the spline basis; a recommended package),
`tibble`, `stats`, `utils`, `grid`, `grDevices`, `parallel`. `Rcpp` joins Imports and `Rcpp`,
`RcppArmadillo` join LinkingTo when 1A-i adds `src/` (R CMD check rejects LinkingTo without `src/`, so the skeleton leaves them out). Suggests:
`ggplot2`, `gtable`, `patchwork`, `png`, `ragg`, `fastpng`, `scales`, `systemfonts` (plot
helpers) and `testthat (>= 3.0.0)`. No lubridate, dplyr or purrr: the package does its date and
table work in base R and tibble. No GitHub remotes. Pipelines install the package as
`github::mt-climate-office/normals@<sha>` (spec §2 rule 6); its own CRAN dependencies are pinned
by each pipeline's Posit Package Manager snapshot (plan 0D).

**CI** (`.github/workflows/R-CMD-check.yaml`): one job on `ubuntu-24.04`, R 4.5, actions pinned
to release tags (`actions/checkout@v7.0.1`, `r-lib/actions/*@v2.14.0`), CRAN from a dated Posit
Package Manager snapshot, `R CMD check --no-manual` with warnings as errors. Runs on
pushes and pull requests to `main` and `v2`.

## 9. Open questions

**Q1. `por` start: whole archive or each subject?** Spec §5 says `por` runs from "each
subject's first year". mco-snotel uses the archive's first year for every station
(`R/normals.R:44`), so a station's earlier years are NA columns and `min_valid` = 10 for all.
Under the spec's reading a station with 6 years would get min_valid = 6 and a normal. Phase 1
reproduces the code. **For the owner, before 2.1:** which is intended? If the code is, §5's
wording should change; if the spec is, `por` normals change at the v1 switch.

**Q2. SNODAS Sep 30 uses WY2003.** mco-snodas's reference dates for Sep 30 include 2003-09-30
(P1), so that day has 24 reference years while the spec's `snodas` role (`earliest: 2004`), the
plots' "WY2004–" label and mco-snotel's `snodas` period have 23. Phase 1 reproduces it (and
whether the 2003-09-30 COG exists in production should be confirmed). **For the owner:** fix at
the v1 switch (drop it), or keep?

**Q3. Who evaluates role rules?** Spec §1 gives roles to the resolver; the plan lists
`period_years` under 1A-ii. Proposal: the package exports the pure `period_years()` and
`period_key()`; the resolver calls them, so the rule logic exists once. The alternative is to
leave both in the resolver and drop them here. Orchestrator's call, with 2A.

**Q4. Plot packages in Suggests, not Imports.** Keeps the core light for the gridded workers and
other users; both pipeline images install them anyway. The cost is a runtime check in each plot
helper. Orchestrator's call.

**Q5. Where the pixel fixture lives.** ~29 MB of statistics plus the date stacks. Options:
(a) commit it under `tests/fixtures-large/`, build-ignored, run by a CI step; (b) a GitHub
release asset downloaded by CI and checksummed; (c) shrink to a 50 × 50 block of 8 dates plus a
10 × 10 block for the smoothed fit. The brief asks for 50 × 50; (b) is the proposal.

**Q6. Spec vocabulary that is SWE- or water-year-specific.** `fit_status = "no_snow"` is a ZIG
status in a generic package (a gridMET precipitation family would mean "no precipitation"), and
`ref_first_wy` / `ref_last_wy` name water years even for a `cy` period. Proposal: keep the spec
names in 2.0 (they are what pipelines publish) and revisit when gridMET adopts the spec (§14).

**Q7. `date` and `day` in `normals_table()`.** The spec's table has `day` only; Phase 1 pipelines
key on `date`. 2.0 returns both. Should 2.1 drop `date`?

**Q8. Unifying at the v1 switch (2.1).** D1 (fallback everywhere?), D2 (mask everywhere?), D3
(what n counts), D6 (Acklam vs exact `qnorm`; the published z is rounded to 0.01 anyway), P3
(float32 parameters on the grid) and Q2 are the decisions the plan's "methods unify at the v1
switch" needs. Listed here so 3B starts from them.

**Q9. Generic smoothing engines.** Should 2.0 export the penalised binomial and WLS curve fitters
(today inside `zig_spline.cpp`) so the next family can reuse them? Proposal: not in 2.0.0; add
them with the first second family, in 2.x, without changing ZIG's path.

**Q10. DESCRIPTION authors.** The skeleton lists Kyle Bocinsky (aut, cre) and Colin Brust (aut,
for 1.x). The owner should confirm.

## 10. Surprises found in the code

1. **The grid's smoothed flavor has two fallback rules** (D1a vs D1b): its percentiles fall
   back to the ECDF (`zig_p(…, hist_mat)`) but its quantiles and mean do not (`zig_smooth_q`).
   Zonal and SNOTEL fall back in both.
2. **Grid percentiles use float32 parameters** (P3) for both `zig` (nightly) and `zig_smooth`,
   while the grid's own annual `zig` quantiles use doubles. A percentile and a quantile COG for
   the same pixel and day can therefore disagree in the last float32 bit.
3. **Sep 30 in SNODAS includes WY2003** (Q2).
4. **SNOTEL `por` starts at the archive's first year** (Q1).
5. **`min_valid` counts columns, not years in the period** (D3): a missing SNODAS COG lowers the
   threshold for that day; in SNOTEL an all-NA year still counts.
6. **Two standard normal quantile functions** (D6): published z is Acklam's approximation
   rounded to 0.01; the outlook uses exact `qnorm`.
7. **All-missing zone means are NaN, not NA** (P6), and are written as NaN.
8. **`q1000` is always NA for fitted and smoothed rows with a positive part** (`zig.cpp:262`,
   by design: the gamma has no maximum). Consumers reading `q1000` get NA except for the
   empirical flavor and all-zero rows.
9. **The outlook's ECDF p_zero hard-codes 10 years** (D7), unlike everything else.
