# Reference for the pixel tests (test-plot-pixels.R): the code the plot
# helpers were ported from, copied verbatim, so the tests render each plot
# with today's pipeline code and with the package and compare the pixels.
# Do not edit, other than to re-copy from upstream.
#
# 1. mco-snodas R/plot.R lines 1-278 (the shared helpers, through
#    snodas_conditions()), at 9b45d74 (file last changed in 62fe19a).
#    mco-snotel vendors the same lines (outlook:R/plot.R).
# Pin a plot's frame so panels sit in the same place on every plot and HUC,
# and flipping between plots moves as little as possible. Takes a
# ggplotGrob() or a patchworkGrob() (whose first panel's columns and rows
# carry a "-1" suffix).
#   axis_width:  the left axis column. 0.45 in fits the widest 8.8 pt tick
#                labels (0.34 in for "14000" or "1.25M") plus tick and margin.
#   xlab_height: the x-axis title row. A bold 11 pt title needs 0.18 in. The
#                seasonal plot has no title but two-line date labels, and 5 pt
#                puts its panel's bottom edge where the hypsomeSWE panels' is.
fix_plot_frame <- function(g, axis_width  = grid::unit(0.45, "in"),
                           xlab_height = grid::unit(0.2, "in")) {
  col <- unique(g$layout$l[g$layout$name %in% c("axis-l", "axis-l-1")])
  stopifnot(length(col) == 1L)
  g$widths[col] <- axis_width
  row <- unique(g$layout$t[grepl("^xlab-b(-\\d+)?$", g$layout$name)])
  g$heights[row] <- xlab_height
  g
}

# Normal-condition bands shared by the seasonal and hypsomeSWE plots. The
# outer limits are the 1st and 99th percentiles for every normals flavor:
# ECDF cannot resolve beyond them with ~23 reference years, and the fitted
# flavors have no record min/max.
snodas_bands <- function(usdm = FALSE) {
  if (usdm) {
    list(
      bands  = list(D4 = c("q0010", "q0020"), D3 = c("q0020", "q0050"), D2 = c("q0050", "q0100"),
                    D1 = c("q0100", "q0200"), D0 = c("q0200", "q0300"), Neutral = c("q0300", "q0700"),
                    W0 = c("q0700", "q0800"), W1 = c("q0800", "q0900"), W2 = c("q0900", "q0950"),
                    W3 = c("q0950", "q0980"), W4 = c("q0980", "q0990")),
      breaks = c(0, 2, 5, 10, 20, 30, 70, 80, 90, 95, 98, 100),   # percentile classes
      colors = c(D4 = "#730000", D3 = "#E60000", D2 = "#FFAA00", D1 = "#FCD37F", D0 = "#FFFF00",
                 Neutral = "white", W0 = "#c6f0ff", W1 = "#73b2ff", W2 = "#0070ff",
                 W3 = "#0000c8", W4 = "#00007d"),
      legend = c("D4", "D3", "D2", "D1", "D0", "W0", "W1", "W2", "W3", "W4"),
      ncol   = 2
    )
  } else {
    list(
      bands  = list(`1st–5th` = c("q0010", "q0050"), `5th–30th` = c("q0050", "q0300"),
                    `30th–70th` = c("q0300", "q0700"), `70th–95th` = c("q0700", "q0950"),
                    `95th–99th` = c("q0950", "q0990")),
      breaks = c(0, 5, 30, 70, 95, 100),
      colors = c(`1st–5th` = "#730000", `5th–30th` = "orange", `30th–70th` = "green",
                 `70th–95th` = "cyan", `95th–99th` = "#0A0099"),
      legend = c("1st–5th", "5th–30th", "30th–70th", "70th–95th", "95th–99th"),
      ncol   = 1
    )
  }
}

# Reference water years for the normals of the water year containing date.
snodas_ref_years <- function(date) {
  wy <- water_year(max(date))
  c(max(wy - 30L, 2004L), wy - 1L)
}

# Seasonal plot date axis: a tick and gridline on the first of every month.
snodas_month_breaks <- function(limits) {
  b <- seq(lubridate::floor_date(limits[1], "month"), limits[2], by = "month")
  b[b >= limits[1]]
}

# Month labels with the year under each, so any month reads on its own.
snodas_month_labels <- function(x) format(x, "%b\n%Y")

# Legend sizing shared by the seasonal and hypsomeSWE plots so they read as a set.
snodas_legend_theme <- function()
  ggplot2::theme(legend.margin     = ggplot2::margin(0, 0, 0, 0),
                 legend.key.size   = grid::unit(14, "pt"),
                 legend.text       = ggplot2::element_text(size = 9),
                 # 11 pt: at 10 pt bold and 150 dpi, Helvetica's "i" loses its dot.
                 legend.title      = ggplot2::element_text(size = 11, face = "bold"),
                 legend.box.spacing = grid::unit(6, "pt"))

# Shrink an image (rows x cols x RGB[A] array, 0-1) to out_h rows by averaging
# the source pixels each output pixel covers, weighted by alpha. A device's
# raster scaling samples instead, which leaves logo edges jagged.
snodas_resample <- function(img, out_h) {
  if (out_h >= dim(img)[1]) return(img)
  out_w <- max(1L, round(dim(img)[2] * out_h / dim(img)[1]))
  # n_out x n_in weights: the share of each output pixel each input pixel covers.
  overlap <- function(n_in, n_out) {
    e_in <- 0:n_in / n_in
    e_out <- 0:n_out / n_out
    m <- outer(seq_len(n_out), seq_len(n_in), \(o, i)
      pmax(0, pmin(e_out[o + 1], e_in[i + 1]) - pmax(e_out[o], e_in[i])))
    m / rowSums(m)
  }
  wr <- overlap(dim(img)[1], out_h)
  wc <- t(overlap(dim(img)[2], out_w))
  a  <- if (dim(img)[3] == 4) img[, , 4] else matrix(1, dim(img)[1], dim(img)[2])
  sa <- wr %*% a %*% wc
  out <- array(0, c(out_h, out_w, 4))
  for (k in 1:3) out[, , k] <- ifelse(sa > 0, (wr %*% (img[, , k] * a) %*% wc) / sa, 0)
  out[, , 4] <- sa
  pmin(pmax(out, 0), 1)
}

# A logo drawn `height` inches tall, resampled with snodas_resample() to the
# device's pixel size when it is drawn, so it is clean at any resolution (the
# plots' and the animations' no-snow frames). `...` goes to rasterGrob().
snodas_logo_grob <- function(img, height, ...)
  grid::gTree(img = img, height = height, raster_args = list(...), cl = "snodas_logo")

makeContent.snodas_logo <- function(x) {
  dpi <- grDevices::dev.size("px")[2] / grDevices::dev.size("in")[2]
  r   <- snodas_resample(x$img, max(1L, round(x$height * dpi)))
  grid::setChildren(x, grid::gList(do.call(grid::rasterGrob, c(
    list(r, height = grid::unit(x$height, "in"), interpolate = TRUE), x$raster_args))))
}

# The MCO logo image, read once per session. Read with png, not terra, because
# some terra versions flip an ungeoreferenced PNG. The 3600 px original is
# shrunk to 400 px tall on reading, which keeps forked plotting workers light
# and is still twice what a 300 dpi plot draws (0.55 in).
snodas_logo <- local({
  img <- NULL
  function(path = "mco_logo.png") {
    if (is.null(img)) img <<- snodas_resample(png::readPNG(path), 400L)
    img
  }
})

# The University of Montana wordmark, read once per session, for the top-right
# corner of every plot (snodas_legend_row()), centred on the title and
# subtitle. UMT_LOGO_HEIGHT is in inches.
UMT_LOGO_HEIGHT <- 0.3
snodas_umt_logo <- local({
  img <- NULL
  function(path = "umt_logo.png") {
    if (is.null(img)) img <<- png::readPNG(path)
    img
  }
})

# Title size for every plot, set so the longest title in the domain
# (Upper Cumberland-Cordell Hull Reservoir, 05130106) clears the UM logo by
# ~0.3 in at 8 in wide (12 pt leaves ~0.07 in).
SNODAS_TITLE_SIZE <- 11.5

# The "Current Conditions" key shared by both plots: this year's line, then
# the median's gray dashed line. Both keys show on every plot, even one with no
# observed SWE, so the legend row is the same from day to day.
snodas_current_legend <- function(obs_lwd)
  ggplot2::scale_linetype_manual(
    name = "Current Conditions",
    values = c(Observed = "solid", Median = "dashed"),
    limits = c("Observed", "Median"),
    drop = FALSE,
    # One key per row, Median under Observed; keys wide enough to show dashes.
    guide = ggplot2::guide_legend(order = 3, title.position = "top", ncol = 1,
                                  keywidth = grid::unit(24, "pt"),
                                  override.aes = list(colour = c("black", "gray50"),
                                                      linewidth = c(obs_lwd, 0.5)))
  )

# Rebuild a plot's bottom legend row: the first key (normal conditions)
# centred on the figure, the MCO logo centred in the space to its left, and
# the remaining keys (present conditions) centred in the space to its right.
# Takes a ggplotGrob() or patchworkGrob() and returns it wrapped so it prints
# and ggsave()s like a ggplot.
snodas_legend_row <- function(g) {
  i   <- which(g$layout$name %in% c("guide-box-bottom", "guide-box"))
  box <- g$grobs[[i]]
  keys <- box$grobs[box$layout$name == "guides"]

  centred <- function(grob, width, height) {
    cell <- gtable::gtable(widths  = grid::unit.c(grid::unit(1, "null"), width,  grid::unit(1, "null")),
                           heights = grid::unit.c(grid::unit(1, "null"), height, grid::unit(1, "null")))
    gtable::gtable_add_grob(cell, grob, t = 2, l = 2, clip = "off")
  }
  key_cell <- function(k) centred(k, sum(k$widths), sum(k$heights))

  normal <- keys[[1]]
  h      <- sum(normal$heights)
  logo   <- snodas_logo_grob(snodas_logo(), height = 0.55)   # same on every plot
  right  <- if (length(keys) > 1)
    do.call(gtable::gtable_row, c(list(name = "present", grobs = lapply(keys[-1], key_cell)),
                                  list(widths = do.call(grid::unit.c, lapply(keys[-1], \(k) sum(k$widths))))))

  row <- gtable::gtable(widths  = grid::unit.c(grid::unit(1, "null"), sum(normal$widths), grid::unit(1, "null")),
                        heights = h)
  row <- gtable::gtable_add_grob(row, logo, t = 1, l = 1, clip = "off")
  row <- gtable::gtable_add_grob(row, normal, t = 1, l = 2, clip = "off")
  if (!is.null(right))
    row <- gtable::gtable_add_grob(row, centred(right, sum(right$widths), h), t = 1, l = 3, clip = "off")

  # The row spans the full figure width so "centred" means on the figure.
  g$grobs[[i]]  <- row
  g$layout$l[i] <- 1
  g$layout$r[i] <- ncol(g)
  # ggplot and patchwork size the legend box differently (0.68 vs 0.83 in for
  # the same keys); sizing the row to its content keeps it at one height.
  g$heights[g$layout$t[i]] <- h

  # The UM logo, right-aligned at the plot margin and centred vertically
  # between the top of the image and the top of the panels.
  rows <- seq_len(min(g$layout$t[grepl("^panel", g$layout$name)]) - 1L)
  g <- gtable::gtable_add_grob(
    g, snodas_logo_grob(snodas_umt_logo(), height = UMT_LOGO_HEIGHT,
                        x = grid::unit(1, "npc") - grid::unit(5.5, "pt"), hjust = 1),
    t = min(rows), b = max(rows), l = 1, r = ncol(g), clip = "off", name = "umt-logo")
  patchwork::wrap_elements(full = g)
}

# Plots are written as 8-bit palette PNGs on one fixed 256-colour palette
# (plot_palette.png; scripts/build-plot-palette.R), so a colour is the same in
# every plot: legends, logos and axes are pixel-identical across plots and
# animation frames. A per-image palette (pngquant) shifted them by up to
# ~20/255 from plot to plot. The fixed palette is also more accurate (mean
# error 0.12/255 vs 0.14) and its files ~26% smaller.
snodas_palette <- local({
  pal <- NULL
  function(path = "plot_palette.png") {
    if (is.null(pal)) {
      x   <- png::readPNG(path)
      rgb <- unique(round(cbind(c(t(x[, , 1])), c(t(x[, , 2])), c(t(x[, , 3]))) * 255))
      pal <<- list(rgb = rgb, hex = grDevices::rgb(rgb, maxColorValue = 255))
    }
    pal
  }
})

# Palette indices (0-based, rows x cols) for a nativeRaster: each distinct
# colour goes to its nearest palette colour (RGB distance, no dithering).
# Matches are remembered across plots, since nearly all plots share their
# colours (~5,000 in a plot, ~92% of pixels in 20 of them).
snodas_quantize <- local({
  keys <- integer()
  idx  <- integer()
  function(nr) {
    pal <- snodas_palette()$rgb
    v   <- as.vector(nr)   # nativeRaster: row-major, 0xAABBGGRR
    k   <- bitwAnd(v, 0xFFFFFFL)                  # RGB, alpha dropped (opaque on white)
    u   <- unique(k)
    new <- u[is.na(match(u, keys))]
    if (length(new)) {
      r <- bitwAnd(new, 255L)
      g <- bitwAnd(bitwShiftR(new, 8L), 255L)
      b <- bitwShiftR(new, 16L)
      d <- outer(r, pal[, 1], "-")^2 + outer(g, pal[, 2], "-")^2 + outer(b, pal[, 3], "-")^2
      keys <<- c(keys, new)
      idx  <<- c(idx, max.col(-d, ties.method = "first") - 1L)
    }
    matrix(idx[match(k, keys)], nrow = dim(nr)[1], byrow = TRUE)
  }
})

# Save a plot as an 8-bit palette PNG at the shared 3:2 size: drawn in memory
# with ragg, quantized to the fixed palette, written with fastpng. 1500 x 1000
# (187.5 dpi at 8 x 5.33 in): as sharp as 2400 x 1600 at the 600-1200 CSS px
# consumers show them (CONSUMERS.md).
snodas_save_png <- function(p, file, width = 8, height = 16 / 3, dpi = 187.5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  cap <- ragg::agg_capture(width = width, height = height, units = "in", res = dpi,
                           background = "white")
  tryCatch({
    print(p)
    img <- cap(native = TRUE)
  }, finally = grDevices::dev.off())
  fastpng::write_png(snodas_quantize(img), file, palette = snodas_palette()$hex,
                     use_filter = FALSE, compression_level = 9L)
  invisible(file)
}

# Subtitle conditions from the basin-wide record on one date: percent of the
# basin median and the basin percentile. Percent of normal is undefined when
# the median is zero; with no snow as well, it says so instead.
snodas_conditions <- function(swe, median, ptile) {
  # No snow where the median is no snow: the percentile is the middle of the
  # tie at zero (~50th), which reads as more than it says.
  if (isTRUE(swe == 0) && isTRUE(median == 0)) return("No snow (typical for this date)")
  pct  <- if (isTRUE(median > 0)) paste0(round(100 * swe / median), "% of Median")
  rank <- if (!is.na(ptile)) paste0("~", scales::ordinal(round(ptile)), " Percentile")
  if (is.null(pct) && is.null(rank)) return("No normal for this date")
  if (is.null(pct)) rank else if (is.null(rank)) pct else paste0(pct, " (", rank, ")")
}

# 2. mco-snotel R/plot_snotel.R snotel_outlook_legend(), at outlook (bb8d06a).
# snodas_current_legend() (R/plot.R) with the outlook's likely range as a third
# key: Observed and Median in one column, Likely range beside them, so the
# legend row keeps its two-row height.
snotel_outlook_legend <- function(obs_lwd)
  ggplot2::scale_linetype_manual(
    name = "Current Conditions",
    values = c(Observed = "solid", Median = "dashed", `Likely range` = "solid"),
    limits = c("Observed", "Median", "Likely range"),
    drop = FALSE,
    guide = ggplot2::guide_legend(order = 3, title.position = "top", nrow = 2, byrow = FALSE,
                                  keywidth = grid::unit(24, "pt"),
                                  override.aes = list(colour = c("black", "gray50", NA),
                                                      linewidth = c(obs_lwd, 0.5, 0),
                                                      fill = c(NA, NA, "gray15"),
                                                      alpha = c(1, 1, 0.22)))
  )
