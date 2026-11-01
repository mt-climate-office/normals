# A synthetic SNODAS-like pixel block for the gridded golden fixture.
#
# Integer SWE (mm) for every day from 2003-09-30 (SNODAS's first day) through
# 2024-09-30, for nr x nc pixels. Built with exact double arithmetic only
# (+, *, /, floor, %% on integers below 2^53): no random number generator and
# no log/exp, so every platform and R version regenerates the same block.
# data-raw/make-fixtures.R sources this file; the tests check the block's
# checksum before comparing anything.
#
# Pixel classes, by pixel index p (row-major, 1-based):
#   p %% 10 == 0  outside the domain: NA every day (fit status 3)
#   p %% 10 == 1  the mask grew: NA before WY2011, a snowpack after
#   p %% 10 == 2  never snow: 0 every day (fit status 1)
#   p %% 10 == 3  very sparse snow: a few small values in a few years
#                 (ECDF fallback rows, fit status 2)
#   p %% 10 == 4  sparse snow: small values on a few days of some years
#                 (few days with L-moments, fit status 5)
#   otherwise     a seasonal snowpack whose size and timing follow an
#                 elevation-like field, varying by year and day
grid_dates <- function() seq(as.Date("2003-09-30"), as.Date("2024-09-30"), by = "day")

# Uniform [0, 1) from integers, exactly: a 32-bit LCG step.
grid_unif <- function(a) ((a * 69069 + 1) %% 4294967296) / 4294967296

grid_cube <- function(nr = 20L, nc = 20L) {
  dates <- grid_dates()
  n <- nr * nc
  p <- seq_len(n)
  r <- (p - 1) %/% nc
  cc <- (p - 1) %% nc
  elev <- floor(100 * (r / max(1, nr - 1) + cc / max(1, nc - 1)) / 2)   # 0 .. 100
  lt <- as.POSIXlt(dates)
  wy <- lt$year + 1900 + (lt$mon >= 9)
  t  <- as.numeric(dates - as.Date(sprintf("%d-10-01", wy - 1))) + 1    # day of water year
  di <- as.numeric(dates - dates[1]) + 1                                 # day index
  cls <- p %% 10

  out <- matrix(0, n, length(dates))
  for (i in p) {
    u_y <- grid_unif((i * 7919 + wy * 104729 + 12345) %% 4294967296)   # per year
    u_d <- grid_unif((i * 31337 + di * 2711 + 777) %% 4294967296)      # per day
    if (cls[i] == 0) {
      out[i, ] <- NA_real_
    } else if (cls[i] == 2) {
      out[i, ] <- 0
    } else if (cls[i] == 3) {
      on <- u_y < 0.15 & t == 100
      out[i, ] <- ifelse(on, floor(1 + 6 * u_d), 0)
    } else if (cls[i] == 4) {
      on <- u_y < 0.3 & t >= 100 & t <= 103
      out[i, ] <- ifelse(on, floor(1 + 6 * u_d), 0)
    } else {
      e  <- elev[i]
      s0 <- 60 - floor(e / 4)          # first snow
      tp <- 160 + floor(e / 2)         # peak
      te <- 210 + e                    # melt-out
      shape <- ifelse(t <= s0 | t >= te, 0,
               ifelse(t <= tp, (t - s0) / (tp - s0), (te - t) / (te - tp)))
      amp <- 40 + 8 * e
      fy  <- 0.2 + 1.6 * u_y           # a dry or a big year
      v   <- floor(amp * fy * shape * (0.85 + 0.3 * u_d))
      # Early-season dustings that melt: some years have snow before s0.
      v[t > s0 - 20 & t <= s0 & u_d < 0.15] <- 1
      out[i, ] <- v
      if (cls[i] == 1) out[i, wy < 2011] <- NA_real_
    }
  }
  list(m = out, dates = dates, nr = nr, nc = nc)
}

# An exact checksum of the block (integers below 2^53).
grid_checksum <- function(cube) {
  w <- (seq_along(cube$m) %% 9973) + 1
  c(sum = sum(cube$m * w, na.rm = TRUE), na = sum(is.na(cube$m)))
}

# The block's values on the given dates: n_pixels x length(dates).
grid_at <- function(cube, dates) {
  cube$m[, match(dates, cube$dates), drop = FALSE]
}
