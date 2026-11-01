#!/usr/bin/env Rscript
# Golden fixtures for normals 2.0 (API.md section 7), produced by TODAY's
# code: the pipelines' own R functions, sourced, and their C++, compiled with
# Rcpp::sourceCpp(). The 2.0 package is never used here. The tests in
# tests/testthat/test-golden-*.R recompute every expected output with the 2.0
# calls and compare them exactly.
#
# Usage (from the package root):
#   MCO_SNODAS=/path/to/mco-snodas MCO_SNOTEL=/path/to/mco-snotel \
#     Rscript data-raw/make-fixtures.R [--large]
#
#   MCO_SNODAS  an mco-snodas checkout at 9b45d74 (main; src/ and the R used
#               here unchanged since 62fe19a)
#   MCO_SNOTEL  an mco-snotel checkout whose R/normals.R, R/utils.R,
#               R/zig_spline_normals.R and src/ equal main at f39fe56
#   --large     also write the 50 x 50 pixel fixture to
#               data-raw/fixtures-large/grid50.rds (a release asset; see
#               tests/testthat/test-golden-grid-large.R), not committed
#
# Inputs (downloaded read-only to data-raw/cache/, which is git-ignored):
#   SNOTEL  s3://mco-snotel/normals/reference/wy2027.parquet (the frozen record,
#           WY2027 build) and s3://mco-snotel/archive/snotel.parquet
#           (observations), both fetched 2026-11-01
#   Zonal   s3://mco-snodas/normals/zonal/wbd/raw/{d}/{d}_raw_swe_means.parquet
#           and _raw_hypsomeswe_means.parquet for every date d of WY2026
#           (2025-10-01 .. 2026-09-30), fetched 2026-11-01; the hypsomeswe
#           tables are reduced to HUC 10020001's bands on download
#           (data-raw/cache/zonal_hyps/{d}.rds)
#   Pixels  synthetic: tests/testthat/helper-grid.R (deterministic)
#
# Outputs: tests/testthat/fixtures/{snotel,zonal,grid}.rds (version 3, xz)
# and tests/testthat/fixtures/MANIFEST.json.

suppressPackageStartupMessages({
  library(magrittr)    # today's code uses %>% unqualified
  library(lubridate)   # ... and %m-%
})
options(warn = 1)

SNODAS_SHA <- "9b45d74"
SNOTEL_SHA <- "f39fe56"
SNODAS <- Sys.getenv("MCO_SNODAS")
SNOTEL <- Sys.getenv("MCO_SNOTEL")
LARGE  <- "--large" %in% commandArgs(TRUE)
# Sections to (re)build; the others' fixture files are kept as they are.
ONLY   <- strsplit(Sys.getenv("MCO_FIXTURES", "snotel,zonal,grid"), ",")[[1]]
stopifnot(dir.exists(SNODAS), dir.exists(SNOTEL), file.exists("DESCRIPTION"))

CACHE <- "data-raw/cache"
large <- "data-raw/fixtures-large/grid50.rds"
FIX   <- "tests/testthat/fixtures"
dir.create(FIX, showWarnings = FALSE, recursive = TRUE)
WORK  <- file.path(CACHE, "work")   # scratch for the pipelines' file I/O
unlink(WORK, recursive = TRUE)
dir.create(WORK, recursive = TRUE)

git <- function(repo, ...) system2("git", c("-C", repo, ...), stdout = TRUE)
snodas_head <- git(SNODAS, "rev-parse", "HEAD")
snotel_head <- git(SNOTEL, "rev-parse", "HEAD")
# The files used must equal mco-snodas 9b45d74 and mco-snotel main f39fe56.
SNODAS_FILES <- c("src", "R/utils.R", "R/calc_snodas_zig_spline_normals.R", "R/calc_zonal_normals.R",
                  "R/normals_annual.R", "R/pipeline_sync.R", "R/calc_snodas_normals.R")
SNOTEL_FILES <- c("R/normals.R", "R/utils.R", "R/zig_spline_normals.R", "src")
same_as <- function(repo, sha, files)
  length(git(repo, "diff", "--stat", sha, "HEAD", "--", files)) == 0L &&
    length(git(repo, "status", "--porcelain", "--", files)) == 0L
stopifnot(same_as(SNODAS, SNODAS_SHA, SNODAS_FILES), same_as(SNOTEL, SNOTEL_SHA, SNOTEL_FILES))

save_fixture <- function(x, name) {
  path <- file.path(FIX, paste0(name, ".rds"))
  saveRDS(x, path, version = 3, compress = "xz")
  message("  wrote ", path, " (", format(file.size(path) / 1024, digits = 3), " KiB)")
  path
}

# Evaluate selected top-level assignments of a file into env (for files whose
# other top-level code needs the whole pipeline).
source_defs <- function(file, names, env) {
  exprs <- parse(file, keep.source = FALSE)
  for (e in exprs) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        as.character(e[[2]]) %in% names) eval(e, env)
  }
  missing <- setdiff(names, ls(env))
  if (length(missing)) stop("not found in ", file, ": ", paste(missing, collapse = ", "))
}

source_cpp <- function(repo, env) {
  for (f in c("ecdf", "zig", "zig_spline", "qnorm_grid"))
    Rcpp::sourceCpp(file.path(repo, "src", paste0(f, ".cpp")), env = env, rebuild = TRUE)
}

G1_DAYS <- c("10-01", "11-15", "01-01", "02-29", "04-01", "05-15", "07-01", "09-30")


# ══ SNOTEL (mco-snotel R/normals.R) ═════════════════════════════════════════
# Five stations, normals year WY2024 (a leap year, so Feb 29 is covered),
# every period. The record is the WY2027 frozen record before 2023-10-01
# (what a WY2024 build reads from it); the observations are WY2024's.
#   578:MT:SNTL   long full record (from 1963-09-30, the archive's start)
#   1287:MT:SNTL  6 years in the 30yr period (the D2 mask is active)
#   1084:OR:SNTL  low station: few positive years on many days (D1 fallback,
#                 status 5)
#   932:MT:SNTL   whole years missing inside the record (NA columns)
#   1069:WA:SNTL  record starts WY2008, after WY2004 (snodas vs por)
# Observations of 1000:OR:SNTL are included and must be dropped.

if ("snotel" %in% ONLY) {
message("SNOTEL")
snotel <- new.env()
for (f in c("R/utils.R", "R/zig_spline_normals.R", "R/normals.R"))
  sys.source(file.path(SNOTEL, f), envir = snotel)
source_cpp(SNOTEL, snotel)

STATIONS <- c("1069:WA:SNTL", "1084:OR:SNTL", "1287:MT:SNTL", "578:MT:SNTL", "932:MT:SNTL")
SNOTEL_WY <- 2024L
ref_all <- arrow::read_parquet(file.path(CACHE, "snotel_reference_wy2027.parquet"))
record <- ref_all[ref_all$triplet %in% STATIONS & ref_all$date < as.Date("2023-10-01"), ]
record <- as.data.frame(record[order(record$triplet, record$date), c("triplet", "date", "wteq")])
rownames(record) <- NULL
arch <- arrow::read_parquet(file.path(CACHE, "snotel_archive.parquet"),
                            col_select = c("triplet", "date", "wteq"))
OBS_DATES <- as.Date(c("2023-10-15", "2023-12-01", "2024-01-15", "2024-02-01", "2024-02-29",
                       "2024-03-15", "2024-04-15", "2024-05-15", "2024-06-15", "2024-09-30"))
obs <- arch[arch$triplet %in% c(STATIONS, "1000:OR:SNTL") & arch$date %in% OBS_DATES, ]
obs <- as.data.frame(obs[order(obs$date, obs$triplet), ])
rownames(obs) <- NULL
rm(ref_all, arch); gc()

sm <- snotel$swe_matrix(record)
stopifnot(identical(sm$stations, sort(STATIONS, method = "radix")))
# Every date's normals rows are kept.
keep_dates <- snotel$normals_dates(SNOTEL_WY)
snotel_out <- lapply(stats::setNames(snotel$SNOTEL_PERIODS, snotel$SNOTEL_PERIODS), \(period) {
  res <- snotel$calc_station_normals(sm, SNOTEL_WY, period, cores = 1L)
  pc  <- snotel$calc_station_percentiles(sm, obs, res$normals, res$years)
  list(years = res$years,
       normals = res$normals[res$normals$date %in% keep_dates, ],
       coef = res$coef, diagnostics = res$diagnostics, percentiles = pc)
})
snotel_fixture <- list(
  wy = SNOTEL_WY, stations = STATIONS, record = record, obs = obs,
  keep_dates = keep_dates, expected = snotel_out)
snotel_path <- save_fixture(snotel_fixture, "snotel")
}


# ══ SNODAS zonal (mco-snodas R/calc_zonal_normals.R, R/normals_annual.R) ════
# Normals year WY2026 (its raw tables hold the WY2026 value, so the
# percentiles see a whole winter). WY2026 has no Feb 29; 02-28 stands in for
# the G1 date 02-29 (Feb 29 is covered by the SNOTEL and pixel fixtures).
#   10        a HUC2 (Missouri)
#   10020001  a HUC8 (Red Rock, MT); its hypsomeswe bands too
#   03090102  no snow on record (status 1)
#   03050202  sparse snow: too few positive values (status 2, fallback rows)
#   02080108  few days with L-moments (status 5, fallback rows)
# One band of 10020001 is set to NaN in its first three reference years
# (2004-2006) on every date, as an all-missing band mean is in production
# (0 / 0, API.md P6), so NaN handling is exercised.

message("SNODAS zonal")
snodas <- new.env()
sys.source(file.path(SNODAS, "R/utils.R"), envir = snodas)
sys.source(file.path(SNODAS, "R/calc_snodas_zig_spline_normals.R"), envir = snodas)
sys.source(file.path(SNODAS, "R/calc_zonal_normals.R"), envir = snodas)
source_defs(file.path(SNODAS, "R/normals_annual.R"),
            c("annual_fit_zonal", "spline_fit_dates", "annual_reference_dates"), snodas)
source_defs(file.path(SNODAS, "R/pipeline_sync.R"), "normals_reference_years", snodas)
source_cpp(SNODAS, snodas)
ZONES <- c("02080108", "03050202", "03090102", "10", "10020001")
HYPS_ZONE <- "10020001"
# The pipeline's setup and I/O, stubbed: files stay under WORK.
with(snodas, {
  SNODAS_START <- as.Date("2003-09-30")
  S3_PREFIX_NORMALS_ZONAL_WBD <- "normals/zonal/wbd"
  KEEP_LOCAL <- TRUE
  s3_sync_down <- function(...) invisible(NULL)
  s3_sync_up   <- function(...) invisible(NULL)
})
snodas$OUTPUT_DIR   <- WORK
snodas$wbd_zone_ids <- function() ZONES

if ("zonal" %in% ONLY) {
ZONAL_WY <- 2026L
zonal_dates <- snodas$normals_dates(ZONAL_WY)
g1_zonal <- zonal_dates[format(zonal_dates, "%m-%d") %in% c(G1_DAYS, "02-28")]
raw_dir <- file.path(WORK, "normals/zonal/wbd/raw")
NAN_BAND <- 7000   # swe_elevation (ft) of the band set to NaN in 2004-2006
nan_cols <- function(cols) {
  dts <- as.Date(cols, format = "%Y-%m-%d", optional = TRUE)
  !is.na(dts) & as.integer(format(dts, "%Y")) + (as.integer(format(dts, "%m")) >= 10) %in% 2004:2006
}
raw_swe  <- list()
raw_hyps <- list()
for (d in format(zonal_dates)) {
  z <- arrow::read_parquet(file.path(CACHE, "zonal_raw", paste0(d, "_raw_swe_means.parquet")))
  z <- z[as.character(z$zone) %in% ZONES, ]
  z$zone <- droplevels(z$zone)
  h <- readRDS(file.path(CACHE, "zonal_hyps", paste0(d, ".rds")))
  h$zone <- droplevels(h$zone)
  h[h$swe_elevation == NAN_BAND, nan_cols(names(h))] <- NaN
  dir.create(file.path(raw_dir, d), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(z, file.path(raw_dir, d, paste0(d, "_raw_swe_means.parquet")))
  arrow::write_parquet(h, file.path(raw_dir, d, paste0(d, "_raw_hypsomeswe_means.parquet")))
  # Read back, as the pipeline reads them.
  raw_swe[[d]]  <- arrow::read_parquet(file.path(raw_dir, d, paste0(d, "_raw_swe_means.parquet")))
  raw_hyps[[d]] <- arrow::read_parquet(file.path(raw_dir, d, paste0(d, "_raw_hypsomeswe_means.parquet")))
}
stopifnot(any(is.nan(as.matrix(raw_hyps[[1]][, -(1:3)]))))

# The smoothed fit (annual_fit_zonal(), which writes the zig_smooth tables).
snodas$annual_fit_zonal(ZONAL_WY, cores = 1L)

# The fit itself, as annual_fit_zonal() computes it (it keeps only the
# tables), and the daily statistics it fits.
wy_start <- snodas$water_year_start(zonal_dates[1])
hist_of <- function(x) {
  cols <- names(x)[!is.na(suppressWarnings(lubridate::as_date(names(x))))]
  as.matrix(x[, cols[lubridate::as_date(cols) < wy_start], drop = FALSE])
}
zonal_fit <- function(raw) {
  fit_dates <- snodas$spline_fit_dates(zonal_dates)
  vars  <- c(n = "n", z = "z", l1 = "l1", lcv = "lcv")
  stats <- lapply(vars, \(v) matrix(NA_real_, nrow(raw[[1]]), length(fit_dates)))
  for (j in seq_along(fit_dates)) {
    ss <- snodas$zig_spline_suffstats(hist_of(raw[[format(fit_dates[j])]]))
    for (v in vars) stats[[v]][, j] <- ss[, v]
  }
  fit <- snodas$zig_spline_fit_parallel(snodas$zig_spline_chunk_stats(stats, 1L),
                                        snodas$zig_spline_basis(at = snodas$normals_cycle_day(fit_dates)), 1L)
  list(stats = stats, fit = fit)
}
fit_swe  <- zonal_fit(raw_swe)
fit_hyps <- zonal_fit(raw_hyps)

# ECDF and ZIG normals and every flavor's percentiles, per G1 date
# (calc_zonal_normals(), annual build: write_normals = TRUE; it reads the
# zig_smooth tables written above and adds their percentiles).
zonal_tables <- lapply(stats::setNames(g1_zonal, format(g1_zonal)), \(d) {
  ds <- format(d)
  res <- snodas$calc_zonal_normals(d, zonal_dir = file.path(WORK, "normals/zonal/wbd"),
                                   force = TRUE, write_normals = TRUE,
                                   means = list(hypsomeswe = raw_hyps[[ds]], zone = raw_swe[[ds]]))
  res[c("ecdf", "zig", "zig_smooth")]
})

zonal_fixture <- list(
  wy = ZONAL_WY, zones = ZONES, hyps_zone = HYPS_ZONE, g1_dates = g1_zonal,
  raw_swe = raw_swe,                                 # every date, 5 zones
  raw_hyps = raw_hyps[format(g1_zonal)],             # G1 dates, the bands
  stats_hyps = fit_hyps$stats,                       # the bands' daily statistics
  expected = list(
    tables = zonal_tables,
    fit = list(swe = fit_swe$fit, hypsomeswe = fit_hyps$fit),
    stats_swe = fit_swe$stats))
zonal_path <- save_fixture(zonal_fixture, "zonal")
}


# ══ SNODAS pixels (mco-snodas R/calc_snodas_normals.R, R/normals_annual.R) ═
# A synthetic 20 x 20 block (tests/testthat/helper-grid.R), normals year
# WY2024 (a leap year), run through the gridded build's own functions in
# its order, with its raster I/O replaced by in-memory matrices: the COG
# stack of a date's reference years is the block on normals_reference_years()
# (newest first, as terra stacks the files); the FLT4S parameter files are
# written and read back with terra (float32, API.md P3); the INT2S encoding
# of means and quantiles (P4) is storage and is not applied. One reference
# COG is missing for every date in WY2024's April (the 2015 layer), so those
# dates have one reference year fewer (D3).

source("tests/testthat/helper-grid.R")

grid_fixture_for <- function(nr, nc) {
  cube <- grid_cube(nr, nc)
  GRID_WY <- 2024L
  dates <- snodas$normals_dates(GRID_WY)
  g1 <- dates[format(dates, "%m-%d") %in% G1_DAYS]
  wy_start <- snodas$water_year_start(dates[1])
  missing_cog <- \(d, rd) format(d, "%Y-%m") == "2024-04" & format(rd, "%Y") == "2015"
  layers <- function(d) {
    rd <- snodas$normals_reference_years(d)
    rd[!missing_cog(d, rd)]
  }
  hist_mat <- function(d) {
    rd <- layers(d)
    h <- grid_at(cube, rd[rd < wy_start])
    colnames(h) <- format(rd[rd < wy_start])
    h
  }
  current <- function(d) grid_at(cube, d)[, 1]
  f32 <- function(m) {        # a FLT4S GeoTIFF round trip, as the pipeline's params files
    r <- terra::rast(nrows = nr, ncols = nc, nlyrs = ncol(m), xmin = 0, xmax = nc, ymin = 0, ymax = nr)
    r <- terra::setValues(r, m)
    names(r) <- c("shape", "scale", "p_zero")
    f <- tempfile(tmpdir = WORK, fileext = ".tif")
    terra::writeRaster(r, f, datatype = "FLT4S", overwrite = TRUE)
    terra::rast(f)[]
  }
  P <- unname(snodas$SNODAS_QUANTILES)

  # The smoothed fit (annual_fit_block(), one block, two workers: its
  # cut(, cores) needs cores > 1).
  fit_dates <- snodas$spline_fit_dates(dates)
  vars <- c(n = "n", z = "z", l1 = "l1", lcv = "lcv")
  n_cell <- nr * nc
  cores <- 2L
  rows <- split(seq_len(n_cell), cut(seq_len(n_cell), max(1L, cores), labels = FALSE))
  chunks <- lapply(rows, \(r) lapply(vars, \(v) matrix(NA_real_, length(r), length(fit_dates))))
  for (j in seq_along(fit_dates)) {
    t <- tibble::as_tibble(snodas$zig_lmom(hist_mat(fit_dates[j])))   # the suffstats parquet
    for (i in seq_along(chunks)) for (v in vars) chunks[[i]][[v]][, j] <- t[[v]][rows[[i]]]
  }
  fit <- snodas$zig_spline_fit_parallel(chunks, snodas$zig_spline_basis(at = snodas$normals_cycle_day(fit_dates)), cores)
  status <- cbind(fit$status)[, 1]
  X <- snodas$zig_spline_basis(at = snodas$normals_cycle_day(dates))$X

  per_date <- lapply(stats::setNames(g1, format(g1)), \(d) {
    j <- match(d, dates)
    h <- hist_mat(d)
    cur <- current(d)
    # calc_snodas_normals(), annual build.
    es <- snodas$ecdf_stats(h)
    zig_lm <- snodas$zig_lmom(h)
    zig_params <- snodas$zig_params_from_lmom(zig_lm, ncol(h))
    zig_quants <- snodas$zig_q(zig_params, P, h)
    zig_mean <- snodas$zig_mean(zig_params, h)
    zig_f32 <- f32(zig_params)
    # annual_fit_block(): the date's parameter tile; annual_materialize().
    p <- snodas$zig_spline_params(fit$coef, X[j, , drop = FALSE], status)
    v <- cbind(p$shape[, 1], p$scale[, 1], p$p_zero[, 1])
    v[status == 3, ] <- NA
    smooth_f32 <- f32(v)
    smooth_mean <- snodas$zig_smooth_mean(smooth_f32)
    smooth_q <- snodas$zig_smooth_q(smooth_f32)
    # calc_snodas_normals(), nightly: percentiles against the frozen params.
    ptile <- list(ecdf = snodas$ecdf_p(h, cur),
                  zig = snodas$zig_p(zig_f32, cur, h),
                  zig_smooth = snodas$zig_p(smooth_f32, cur, h))
    list(n_layers = ncol(h),
         ecdf = es, zig = list(lmom = zig_lm, params = zig_params, quantiles = zig_quants,
                               mean = zig_mean),
         zig_f32 = zig_f32,
         smooth = list(params = v, mean = smooth_mean, quantiles = smooth_q),
         smooth_f32 = smooth_f32,
         ptile = lapply(ptile, \(x) as.integer(round(x * 100))),
         standardized = lapply(ptile, \(x) snodas$qnorm_grid(x * 100)))
  })
  list(wy = GRID_WY, nr = nr, nc = nc, checksum = grid_checksum(cube),
       missing_cog = "2015 layer missing for 2024-04-01 .. 2024-04-30",
       expected = list(fit = fit, per_date = per_date))
}

if ("grid" %in% ONLY) {
message("SNODAS pixels")
grid_path <- save_fixture(grid_fixture_for(20L, 20L), "grid")
if (LARGE) {
  dir.create("data-raw/fixtures-large", showWarnings = FALSE)
  saveRDS(grid_fixture_for(50L, 50L), large, version = 3, compress = "xz")
  message("  wrote ", large, " (", format(file.size(large) / 2^20, digits = 3), " MiB)")
  saveRDS(sub(" .*", "", system2("shasum", c("-a", "256", shQuote(large)), stdout = TRUE)),
          file.path(FIX, "grid50-sha256.rds"), version = 3)
}
}


# ══ Manifest ═════════════════════════════════════════════════════════════
snotel_path <- file.path(FIX, "snotel.rds")
zonal_path  <- file.path(FIX, "zonal.rds")
grid_path   <- file.path(FIX, "grid.rds")
sha256 <- function(f) sub(" .*", "", system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE))
si <- sessionInfo()
manifest <- list(
  created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  generator = "data-raw/make-fixtures.R",
  sources = list(
    `mco-snodas` = list(commit = SNODAS_SHA, checkout = snodas_head,
                        files = c("src/*.cpp", "src/ecdf_impl.h", "R/utils.R",
                                  "R/calc_snodas_zig_spline_normals.R", "R/calc_zonal_normals.R",
                                  "R/normals_annual.R (annual_fit_zonal, spline_fit_dates, annual_reference_dates)",
                                  "R/pipeline_sync.R (normals_reference_years)")),
    `mco-snotel` = list(commit = SNOTEL_SHA, checkout = snotel_head,
                        files = c("src/*.cpp", "src/ecdf_impl.h", "R/utils.R", "R/zig_spline_normals.R",
                                  "R/normals.R"))),
  inputs = list(
    snotel = c("s3://mco-snotel/normals/reference/wy2027.parquet",
               "s3://mco-snotel/archive/snotel.parquet"),
    zonal = "s3://mco-snodas/normals/zonal/wbd/raw/{d}/{d}_raw_{swe,hypsomeswe}_means.parquet, d in WY2026",
    pixels = "synthetic, tests/testthat/helper-grid.R",
    input_md5 = as.list(tools::md5sum(file.path(CACHE, c("snotel_reference_wy2027.parquet",
                                                         "snotel_archive.parquet"))))),
  platform = list(
    R = R.version.string, platform = R.version$platform, os = si$running,
    blas = si$BLAS, lapack = si$LAPACK, la_version = La_version(),
    compiler = system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE),
    packages = lapply(c(Rcpp = "Rcpp", RcppArmadillo = "RcppArmadillo", mgcv = "mgcv",
                        arrow = "arrow", terra = "terra", dplyr = "dplyr", lubridate = "lubridate"),
                      \(p) as.character(utils::packageVersion(p)))),
  comparison = paste(
    "identical() on the fixture's platform. On another platform, results that pass",
    "through linear algebra or libm may differ in the last bits: the tests report the",
    "largest relative difference and fail above 1e-12 (API.md section 7)."),
  fixtures = lapply(stats::setNames(c(snotel_path, zonal_path, grid_path),
                                    basename(c(snotel_path, zonal_path, grid_path)))[
                                      file.exists(c(snotel_path, zonal_path, grid_path))],
                    \(f) list(bytes = file.size(f), sha256 = sha256(f))),
  release_assets = if (file.exists(large)) list(grid50.rds = list(
    bytes = file.size(large), sha256 = sha256(large),
    url = "https://github.com/mt-climate-office/normals/releases/download/fixtures-v1/grid50.rds"))
)
jsonlite::write_json(manifest, file.path(FIX, "MANIFEST.json"), auto_unbox = TRUE, pretty = TRUE)
# The platform fingerprint the tests compare against (no JSON reader needed).
saveRDS(list(R = R.version.string, platform = R.version$platform, blas = si$BLAS,
             lapack = si$LAPACK),
        file.path(FIX, "platform.rds"), version = 3)
unlink(WORK, recursive = TRUE)
message("done")
