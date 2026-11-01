#' Shared plot helpers
#'
#' The frame, bands, legends, logos and palette shared by the SNODAS and
#' SNOTEL plots, so every plot reads as one set. Per-dataset plots stay in the
#' pipelines and are built on these. Ported from mco-snodas `R/plot.R` at
#' 62fe19a (through `snodas_conditions()`); a plot drawn with these helpers is
#' pixel-identical to one drawn with that code. The plotting packages
#' (ggplot2, gtable, patchwork, png, ragg, fastpng, scales, systemfonts) are
#' suggested, not imported; each helper stops with the package to install.
#'
#' * `plot_fix_frame()`: pins the left axis column and x-axis title row of a
#'   ggplot or patchwork grob, so panels sit in the same place on every plot.
#'   `axis_width`: 0.45 in fits the widest 8.8 pt tick labels plus tick and
#'   margin. `xlab_height`: a bold 11 pt title needs 0.18 in.
#' * `plot_bands()`: the normal-condition bands, USDM classes (`usdm = TRUE`)
#'   or five percentile classes: a list of `bands` (quantile column pairs),
#'   `breaks` (percentile classes), `colors`, `legend` (keys) and `ncol`. The
#'   outer limits are the 1st and 99th percentiles for every flavor.
#' * `plot_month_breaks()`, `plot_month_labels()`: a monthly date axis, a tick
#'   on the first of every month, labelled with the month and the year under it.
#' * `plot_legend_theme()`: legend sizing shared by every plot.
#' * `plot_logo()`: the MCO (`"mco"`, shrunk to 400 px tall) or University of
#'   Montana (`"umt"`) logo, read once per session from the package's
#'   `inst/plot/`.
#' * `plot_logo_grob()`: a logo `height` inches tall that resamples itself
#'   (area-averaged, alpha-weighted) to the device resolution when it is drawn,
#'   so its edges are clean at any resolution. Its `makeContent()` method is
#'   registered.
#' * `plot_current_legend()`: the "Current Conditions" key: Observed (solid,
#'   black) over Median (dashed, gray). Both keys show on every plot, so the
#'   legend row is the same from day to day. With `likely_range = TRUE`, a
#'   third key, "Likely range" (a gray15 fill at alpha 0.22, as mco-snotel draws
#'   the outlook), sits beside them and the row keeps its two-row height.
#' * `plot_legend_row()`: rebuilds the bottom legend row: the first key (normal
#'   conditions) centred on the figure, the MCO logo centred to its left and
#'   the other keys centred to its right; adds the University of Montana logo
#'   at the top right. Returns the grob wrapped by
#'   [patchwork::wrap_elements()], so it prints like a ggplot.
#' * `plot_palette()`: the fixed 256-color palette every plot is written with
#'   (`rgb`, a 256 x 3 matrix, and `hex`). One palette makes the shared parts
#'   pixel-identical across plots and animation frames.
#' * `plot_quantize()`: 0-based palette indices (rows x cols) for a
#'   `nativeRaster`, each color to its nearest palette color (RGB distance, no
#'   dithering). Matches are remembered for the session.
#' * `plot_save_png()`: draws a plot in memory with ragg, quantizes it to the
#'   palette and writes an 8-bit palette PNG with fastpng (1500 x 1000 at the
#'   defaults). Returns `file`, invisibly.
#' * `plot_conditions()`: the subtitle's "% of Median (~Nth Percentile)" text,
#'   or "No snow (typical for this date)" / "No normal for this date".
#' * `register_plot_fonts()`: makes "Helvetica" available to systemfonts (and
#'   so ragg). macOS has it; on Linux, Nimbus Sans (fonts-urw-base35) is
#'   registered under that name, with a warning if neither is installed.
#'
#' @param g A `ggplotGrob()` or `patchworkGrob()`.
#' @param axis_width,xlab_height Grid units.
#' @param usdm `TRUE` for USDM classes, `FALSE` for five percentile classes.
#'   No default.
#' @param limits Date axis limits.
#' @param x Dates to label; for the `makeContent()` method, a `normals_logo`
#'   grob.
#' @param which `"mco"` or `"umt"`.
#' @param img An image array (rows x cols x RGB\[A\], 0-1).
#' @param height Logo height in inches.
#' @param ... Passed to [grid::rasterGrob()].
#' @param obs_lwd Line width of the observed series.
#' @param likely_range `TRUE` to add the outlook's "Likely range" key.
#' @param nr A `nativeRaster`.
#' @param p A plot.
#' @param file Output path.
#' @param width,dpi Size (in) and resolution.
#' @param value,median,ptile The value, its normal's median and its percentile
#'   on the plotted date.
#' @name plot_helpers
NULL

# Stop with the packages to install when a plot helper's suggested packages
# are missing.
plot_needs <- function(...) {
  pkgs <- c(...)
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("The normals plot helpers need ", paste(pkgs[!ok], collapse = ", "),
         ": install.packages(", paste(deparse(pkgs[!ok]), collapse = ""), ")",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @rdname plot_helpers
#' @export
plot_fix_frame <- function(g, axis_width = grid::unit(0.45, "in"),
                           xlab_height = grid::unit(0.2, "in")) {
  col <- unique(g$layout$l[g$layout$name %in% c("axis-l", "axis-l-1")])
  stopifnot(length(col) == 1L)
  g$widths[col] <- axis_width
  row <- unique(g$layout$t[grepl("^xlab-b(-\\d+)?$", g$layout$name)])
  g$heights[row] <- xlab_height
  g
}

#' @rdname plot_helpers
#' @export
plot_bands <- function(usdm) {
  if (usdm) {
    list(
      bands  = list(D4 = c("q0010", "q0020"), D3 = c("q0020", "q0050"), D2 = c("q0050", "q0100"),
                    D1 = c("q0100", "q0200"), D0 = c("q0200", "q0300"), Neutral = c("q0300", "q0700"),
                    W0 = c("q0700", "q0800"), W1 = c("q0800", "q0900"), W2 = c("q0900", "q0950"),
                    W3 = c("q0950", "q0980"), W4 = c("q0980", "q0990")),
      breaks = c(0, 2, 5, 10, 20, 30, 70, 80, 90, 95, 98, 100),
      colors = c(D4 = "#730000", D3 = "#E60000", D2 = "#FFAA00", D1 = "#FCD37F", D0 = "#FFFF00",
                 Neutral = "white", W0 = "#c6f0ff", W1 = "#73b2ff", W2 = "#0070ff",
                 W3 = "#0000c8", W4 = "#00007d"),
      legend = c("D4", "D3", "D2", "D1", "D0", "W0", "W1", "W2", "W3", "W4"),
      ncol   = 2
    )
  } else {
    # En dashes, written as escapes to keep the source ASCII.
    nm <- c("1st\u20135th", "5th\u201330th", "30th\u201370th", "70th\u201395th", "95th\u201399th")
    list(
      bands  = stats::setNames(list(c("q0010", "q0050"), c("q0050", "q0300"), c("q0300", "q0700"),
                                    c("q0700", "q0950"), c("q0950", "q0990")), nm),
      breaks = c(0, 5, 30, 70, 95, 100),
      colors = stats::setNames(c("#730000", "orange", "green", "cyan", "#0A0099"), nm),
      legend = nm,
      ncol   = 1
    )
  }
}

#' @rdname plot_helpers
#' @export
plot_month_breaks <- function(limits) {
  # The first of limits[1]'s month (lubridate::floor_date(, "month") in the
  # original).
  first <- as.Date(format(limits[1], "%Y-%m-01"))
  b <- seq(first, limits[2], by = "month")
  b[b >= limits[1]]
}

#' @rdname plot_helpers
#' @export
plot_month_labels <- function(x) format(x, "%b\n%Y")

#' @rdname plot_helpers
#' @export
plot_legend_theme <- function() {
  plot_needs("ggplot2")
  ggplot2::theme(legend.margin     = ggplot2::margin(0, 0, 0, 0),
                 legend.key.size   = grid::unit(14, "pt"),
                 legend.text       = ggplot2::element_text(size = 9),
                 # 11 pt: at 10 pt bold and 150 dpi, Helvetica's "i" loses its dot.
                 legend.title      = ggplot2::element_text(size = 11, face = "bold"),
                 legend.box.spacing = grid::unit(6, "pt"))
}

# Shrink an image (rows x cols x RGB[A] array, 0-1) to out_h rows by averaging
# the source pixels each output pixel covers, weighted by alpha. A device's
# raster scaling samples instead, which leaves logo edges jagged.
plot_resample <- function(img, out_h) {
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

#' @rdname plot_helpers
#' @export
plot_logo_grob <- function(img, height, ...)
  grid::gTree(img = img, height = height, raster_args = list(...), cl = "normals_logo")

#' @rdname plot_helpers
#' @exportS3Method grid::makeContent
makeContent.normals_logo <- function(x) {
  dpi <- grDevices::dev.size("px")[2] / grDevices::dev.size("in")[2]
  r   <- plot_resample(x$img, max(1L, round(x$height * dpi)))
  grid::setChildren(x, grid::gList(do.call(grid::rasterGrob, c(
    list(r, height = grid::unit(x$height, "in"), interpolate = TRUE), x$raster_args))))
}

#' @rdname plot_helpers
#' @export
plot_logo <- local({
  imgs <- list()
  function(which) {
    which <- match.arg(which, c("mco", "umt"))
    if (is.null(imgs[[which]])) {
      plot_needs("png")
      path <- system.file("plot", paste0(which, "_logo.png"), package = "normals",
                          mustWork = TRUE)
      img <- png::readPNG(path)
      # Read with png, not terra (some terra versions flip an ungeoreferenced
      # PNG). The 3600 px MCO original is shrunk to 400 px tall, which keeps
      # forked plotting workers light and is still twice what a 300 dpi plot
      # draws (0.55 in).
      if (which == "mco") img <- plot_resample(img, 400L)
      imgs[[which]] <<- img
    }
    imgs[[which]]
  }
})

#' @rdname plot_helpers
#' @export
plot_current_legend <- function(obs_lwd, likely_range = FALSE) {
  plot_needs("ggplot2")
  if (!likely_range) {
    return(ggplot2::scale_linetype_manual(
      name = "Current Conditions",
      values = c(Observed = "solid", Median = "dashed"),
      limits = c("Observed", "Median"),
      drop = FALSE,
      # One key per row, Median under Observed; keys wide enough to show dashes.
      guide = ggplot2::guide_legend(order = 3, title.position = "top", ncol = 1,
                                    keywidth = grid::unit(24, "pt"),
                                    override.aes = list(colour = c("black", "gray50"),
                                                        linewidth = c(obs_lwd, 0.5)))
    ))
  }
  # Observed and Median in one column, Likely range beside them, so the legend
  # row keeps its two-row height (mco-snotel's snotel_outlook_legend()).
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
}

#' @rdname plot_helpers
#' @export
plot_legend_row <- function(g) {
  plot_needs("gtable", "patchwork", "png")
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
  logo   <- plot_logo_grob(plot_logo("mco"), height = 0.55)   # same on every plot
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
    g, plot_logo_grob(plot_logo("umt"), height = PLOT_UMT_LOGO_HEIGHT,
                      x = grid::unit(1, "npc") - grid::unit(5.5, "pt"), hjust = 1),
    t = min(rows), b = max(rows), l = 1, r = ncol(g), clip = "off", name = "umt-logo")
  patchwork::wrap_elements(full = g)
}

#' @rdname plot_helpers
#' @export
plot_palette <- local({
  pal <- NULL
  function() {
    if (is.null(pal)) {
      plot_needs("png")
      x   <- png::readPNG(system.file("plot", "plot_palette.png", package = "normals",
                                      mustWork = TRUE))
      rgb <- unique(round(cbind(c(t(x[, , 1])), c(t(x[, , 2])), c(t(x[, , 3]))) * 255))
      pal <<- list(rgb = rgb, hex = grDevices::rgb(rgb, maxColorValue = 255))
    }
    pal
  }
})

#' @rdname plot_helpers
#' @export
plot_quantize <- local({
  keys <- integer()
  idx  <- integer()
  function(nr) {
    pal <- plot_palette()$rgb
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

#' @rdname plot_helpers
#' @export
plot_save_png <- function(p, file, width = 8, height = 16 / 3, dpi = 187.5) {
  plot_needs("ragg", "fastpng", "png")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  cap <- ragg::agg_capture(width = width, height = height, units = "in", res = dpi,
                           background = "white")
  tryCatch({
    print(p)
    img <- cap(native = TRUE)
  }, finally = grDevices::dev.off())
  fastpng::write_png(plot_quantize(img), file, palette = plot_palette()$hex,
                     use_filter = FALSE, compression_level = 9L)
  invisible(file)
}

#' @rdname plot_helpers
#' @export
plot_conditions <- function(value, median, ptile) {
  # No snow where the median is no snow: the percentile is the middle of the
  # tie at zero (~50th), which reads as more than it says.
  if (isTRUE(value == 0) && isTRUE(median == 0)) return("No snow (typical for this date)")
  pct  <- if (isTRUE(median > 0)) paste0(round(100 * value / median), "% of Median")
  rank <- if (!is.na(ptile)) {
    plot_needs("scales")
    paste0("~", scales::ordinal(round(ptile)), " Percentile")
  }
  if (is.null(pct) && is.null(rank)) return("No normal for this date")
  if (is.null(pct)) rank else if (is.null(rank)) pct else paste0(pct, " (", rank, ")")
}

#' @rdname plot_helpers
#' @export
register_plot_fonts <- function() {
  plot_needs("systemfonts")
  fonts <- systemfonts::system_fonts()
  if (any(fonts$family == "Helvetica")) return(invisible(NULL))

  nimbus <- fonts[fonts$family == "Nimbus Sans", , drop = FALSE]
  if (nrow(nimbus) == 0L) {
    warning(
      "Helvetica not found and Nimbus Sans not available; ",
      "font rendering may differ from macOS. ",
      "Install fonts-urw-base35 to fix this.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  get_path <- function(style) {
    p <- nimbus$path[nimbus$style == style]
    if (length(p) == 0L) nimbus$path[[1L]] else p[[1L]]
  }

  systemfonts::register_font(
    name       = "Helvetica",
    plain      = get_path("Regular"),
    bold       = get_path("Bold"),
    italic     = get_path("Oblique"),
    bolditalic = get_path("Bold Oblique")
  )
  invisible(NULL)
}

#' Plot sizes
#'
#' `PLOT_TITLE_SIZE`: the title size (pt) of every plot, set so the longest
#' HUC name clears the University of Montana logo. `PLOT_UMT_LOGO_HEIGHT`: that
#' logo's height (in).
#'
#' @format Numeric scalars.
#' @export
PLOT_TITLE_SIZE <- 11.5

#' @rdname PLOT_TITLE_SIZE
#' @export
PLOT_UMT_LOGO_HEIGHT <- 0.3
