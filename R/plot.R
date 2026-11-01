#' Shared plot helpers
#'
#' The frame, bands, legends, logos and palette shared by the SNODAS and
#' SNOTEL plots, so every plot reads as one set. Per-dataset plots stay in the
#' pipelines and are built on these. The plotting packages (ggplot2, gtable,
#' patchwork, png, ragg, fastpng, scales, systemfonts) are suggested, not
#' imported; each helper stops with the package to install.
#'
#' * `plot_fix_frame()`: pins the left axis column and x-axis title row of a
#'   ggplot or patchwork grob, so panels sit in the same place on every plot.
#' * `plot_bands()`: the normal-condition bands, USDM classes (`usdm = TRUE`)
#'   or five percentile classes, with breaks, colors and legend keys.
#' * `plot_month_breaks()`, `plot_month_labels()`: a monthly date axis.
#' * `plot_legend_theme()`: legend sizing shared by every plot.
#' * `plot_logo()`: the MCO (`"mco"`) or University of Montana (`"umt"`) logo,
#'   read once per session from the package's `inst/plot/`.
#' * `plot_logo_grob()`: a logo grob that resamples itself to the device
#'   resolution when drawn (its `makeContent()` method is registered).
#' * `plot_current_legend()`: the "Current Conditions" key.
#' * `plot_legend_row()`: rebuilds the bottom legend row with the logos.
#' * `plot_palette()`, `plot_quantize()`, `plot_save_png()`: 8-bit palette
#'   PNGs on one fixed 256-color palette.
#' * `plot_conditions()`: the subtitle's "% of Median (~Nth Percentile)" text.
#' * `register_plot_fonts()`: makes "Helvetica" available (Nimbus Sans on
#'   Linux).
#'
#' @param g A `ggplotGrob()` or `patchworkGrob()`.
#' @param axis_width,xlab_height Grid units.
#' @param usdm `TRUE` for USDM classes. No default.
#' @param limits Date axis limits.
#' @param x Dates to label; for the `makeContent()` method, a `normals_logo`
#'   grob.
#' @param which `"mco"` or `"umt"`.
#' @param img An image array (rows x cols x RGB\[A\], 0-1).
#' @param height Logo height in inches.
#' @param ... Passed to [grid::rasterGrob()].
#' @param obs_lwd Line width of the observed series.
#' @param nr A `nativeRaster`.
#' @param p A plot.
#' @param file Output path.
#' @param width,dpi Size and resolution.
#' @param value,median,ptile The basin's value, its normal's median and its
#'   percentile on the plotted date.
#' @name plot_helpers
NULL

#' @rdname plot_helpers
#' @export
plot_fix_frame <- function(g, axis_width = grid::unit(0.45, "in"),
                           xlab_height = grid::unit(0.2, "in")) not_implemented()

#' @rdname plot_helpers
#' @export
plot_bands <- function(usdm) not_implemented()

#' @rdname plot_helpers
#' @export
plot_month_breaks <- function(limits) not_implemented()

#' @rdname plot_helpers
#' @export
plot_month_labels <- function(x) not_implemented()

#' @rdname plot_helpers
#' @export
plot_legend_theme <- function() not_implemented()

#' @rdname plot_helpers
#' @export
plot_logo <- function(which) not_implemented()

#' @rdname plot_helpers
#' @export
plot_logo_grob <- function(img, height, ...) not_implemented()

#' @rdname plot_helpers
#' @export
plot_current_legend <- function(obs_lwd) not_implemented()

#' @rdname plot_helpers
#' @export
plot_legend_row <- function(g) not_implemented()

#' @rdname plot_helpers
#' @export
plot_palette <- function() not_implemented()

#' @rdname plot_helpers
#' @export
plot_quantize <- function(nr) not_implemented()

#' @rdname plot_helpers
#' @export
plot_save_png <- function(p, file, width = 8, height = 16 / 3, dpi = 187.5) not_implemented()

#' @rdname plot_helpers
#' @export
plot_conditions <- function(value, median, ptile) not_implemented()

#' @rdname plot_helpers
#' @export
register_plot_fonts <- function() not_implemented()

#' @rdname plot_helpers
#' @exportS3Method grid::makeContent
makeContent.normals_logo <- function(x) not_implemented()

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
