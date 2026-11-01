# The plot helpers against the code they were ported from
# (fixtures/reference-plot.R, verbatim from mco-snodas and mco-snotel). Each
# representative plot is drawn twice, once with the reference helpers and once
# with the package's, written as the pipelines write it (ragg, the fixed
# palette, fastpng), and the decoded RGBA arrays must be identical.
#
# Fonts: both renders use "Helvetica" as the pipelines do, after
# register_plot_fonts(). Where neither Helvetica nor Nimbus Sans is installed,
# both fall back to the same system font, so the comparison still holds.

plot_pkgs <- c("ggplot2", "gtable", "patchwork", "png", "ragg", "fastpng",
               "scales", "systemfonts", "lubridate")

skip_if_no_plot_pkgs <- function() {
  for (p in plot_pkgs) skip_if_not_installed(p)
}

# The reference helpers, in their own environment. Their logo and palette
# readers take paths relative to the working directory (the pipeline's repo
# root), so the run happens in a directory holding copies of the package's
# files; makeContent.snodas_logo() is registered as the pipeline's global
# definition would be found.
reference_file <- normalizePath(test_path("fixtures", "reference-plot.R"))

reference_env <- local({
  env <- NULL
  function() {
    if (is.null(env)) {
      env <<- new.env(parent = globalenv())
      sys.source(reference_file, envir = env, keep.source = FALSE)
      registerS3method("makeContent", "snodas_logo", env$makeContent.snodas_logo,
                       envir = asNamespace("grid"))
    }
    env
  }
})

in_reference_dir <- function(code) {
  dir <- tempfile("normals-plot-")
  dir.create(dir)
  for (f in c("mco_logo.png", "umt_logo.png", "plot_palette.png"))
    file.copy(system.file("plot", f, package = "normals"), file.path(dir, f))
  old <- setwd(dir)
  on.exit({ setwd(old); unlink(dir, recursive = TRUE) })
  force(code)
}

# The helpers a plot is built from, by role: the reference's or the package's.
helper_set <- function(which) {
  if (which == "reference") {
    r <- reference_env()
    list(bands = r$snodas_bands, month_breaks = r$snodas_month_breaks,
         month_labels = r$snodas_month_labels, legend_theme = r$snodas_legend_theme,
         current_legend = function(obs_lwd, likely_range = FALSE)
           if (likely_range) r$snotel_outlook_legend(obs_lwd) else r$snodas_current_legend(obs_lwd),
         fix_frame = r$fix_plot_frame, legend_row = r$snodas_legend_row,
         conditions = r$snodas_conditions, save_png = r$snodas_save_png,
         title_size = r$SNODAS_TITLE_SIZE)
  } else {
    list(bands = plot_bands, month_breaks = plot_month_breaks,
         month_labels = plot_month_labels, legend_theme = plot_legend_theme,
         current_legend = plot_current_legend,
         fix_frame = plot_fix_frame, legend_row = plot_legend_row,
         conditions = plot_conditions, save_png = plot_save_png,
         title_size = PLOT_TITLE_SIZE)
  }
}

# A synthetic water year: normal quantiles from a seasonal curve, and this
# year's SWE through mid-February.
synthetic_year <- function() {
  date <- seq(as.Date("2025-10-01"), as.Date("2026-09-30"), by = "day")
  t <- seq_along(date)
  base <- pmax(0, sin(pi * (t - 30) / 260)) * 20
  qs <- c(q0010 = 0.2, q0020 = 0.3, q0050 = 0.45, q0100 = 0.6, q0200 = 0.75, q0300 = 0.88,
          q0500 = 1, q0700 = 1.12, q0800 = 1.25, q0900 = 1.4, q0950 = 1.55, q0980 = 1.7,
          q0990 = 1.8)
  x <- data.frame(date = date)
  for (q in names(qs)) x[[q]] <- base * qs[[q]]
  obs <- date <= as.Date("2026-02-15")
  x$swe <- ifelse(obs, base * (0.8 + 0.1 * sin(t / 9)), NA)
  x$swe_ptile <- ifelse(obs, 30 + 10 * sin(t / 9), NA)
  x
}

# mco-snodas plot_snodas() / mco-snotel plot_snotel() in base R, built on the
# helper set h. With `outlook`, the likely range is drawn and keyed as
# mco-snotel does.
seasonal_plot <- function(h, x, usdm, outlook = NULL) {
  gg <- asNamespace("ggplot2")
  b <- h$bands(usdm)
  rec <- x[!is.na(x$swe), ]
  rec <- rec[nrow(rec), ]
  ribbons <- do.call(rbind, lapply(names(b$bands), function(nm)
    data.frame(date = x$date, ribbon = nm,
               ymin = x[[b$bands[[nm]][1]]], ymax = x[[b$bands[[nm]][2]]])))
  ribbons$ribbon <- factor(ribbons$ribbon, levels = rev(names(b$bands)), ordered = TRUE)
  lines <- rbind(data.frame(date = x$date, name = "1st percentile", value = x$q0010),
                 data.frame(date = x$date, name = "99th percentile", value = x$q0990))
  obs <- x[!is.na(x$swe), ]
  ol <- if (!is.null(outlook))
    rbind(data.frame(date = rec$date, q0100 = rec$swe, q0900 = rec$swe), outlook)

  p <- gg$ggplot() +
    gg$labs(title = "Seasonal SWE for Somewhere Basin (HUC: 10020001)",
            subtitle = paste0(format(rec$date, "%B %d, %Y"), " — ",
                              h$conditions(rec$swe, rec$q0500, rec$swe_ptile)),
            caption = "Smoothed History") +
    gg$geom_ribbon(data = ribbons, mapping = gg$aes(x = date, ymin = ymin, ymax = ymax, fill = ribbon),
                   alpha = 0.7, colour = "white", linewidth = 0.15, key_glyph = "rect") +
    gg$geom_line(data = lines, mapping = gg$aes(x = date, y = value, color = name)) +
    gg$geom_line(data = x, mapping = gg$aes(x = date, y = q0500, linetype = "Median"),
                 colour = "gray50") +
    gg$scale_fill_manual(name = "Historical Conditions WY2004–2025",
                         values = rev(b$colors), limits = rev(b$legend),
                         breaks = if (usdm) c("D4", "D3", "D2", "D1", "D0", "W4", "W3", "W2", "W1", "W0")
                                  else b$legend,
                         na.value = "#FFFFFF80",
                         guide = gg$guide_legend(order = 1, nrow = if (usdm) 2 else 1, byrow = TRUE,
                                                 title.position = "top",
                                                 override.aes = list(alpha = 0.7))) +
    gg$scale_color_manual(name = NULL,
                          values = c(`99th percentile` = unname(b$colors[length(b$colors)]),
                                     `1st percentile` = unname(b$colors[1])),
                          breaks = c("99th percentile", "1st percentile"), guide = "none") +
    gg$scale_y_continuous(limits = c(0, NA), expand = gg$expansion(c(0, 0.05), 0),
                          labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                          name = "Snow Water Equivalent (acre-feet)") +
    gg$scale_x_date(breaks = h$month_breaks, labels = h$month_labels,
                    expand = gg$expansion(0, 0), name = NULL)
  if (!is.null(ol))
    p <- p +
      gg$geom_ribbon(data = ol, mapping = gg$aes(x = date, ymin = q0100, ymax = q0900,
                                                 linetype = "Likely range"),
                     fill = "gray15", alpha = 0.22, colour = NA, key_glyph = "rect") +
      gg$geom_line(data = ol, mapping = gg$aes(x = date, y = q0100), colour = "gray25",
                   linewidth = 0.3, linetype = "dotted") +
      gg$geom_line(data = ol, mapping = gg$aes(x = date, y = q0900), colour = "gray25",
                   linewidth = 0.3, linetype = "dotted")
  p <- p +
    gg$geom_line(data = obs, mapping = gg$aes(x = date, y = swe),
                 color = "white", linewidth = 1.8, lineend = "round") +
    gg$geom_line(data = obs, mapping = gg$aes(x = date, y = swe, linetype = "Observed"),
                 color = "black", linewidth = 1, lineend = "round") +
    h$current_legend(obs_lwd = 1, likely_range = !is.null(ol)) +
    gg$theme_bw(base_family = "Helvetica") +
    gg$theme(legend.margin = gg$margin(),
             panel.grid.minor.x = gg$element_blank(),
             legend.position = "bottom",
             plot.title.position = "plot",
             plot.title = gg$element_text(face = "bold", size = h$title_size, hjust = 0),
             plot.subtitle = gg$element_text(face = "bold", hjust = 0),
             plot.caption = gg$element_text(hjust = 0, color = "gray40"),
             plot.caption.position = "plot",
             axis.title = gg$element_text(face = "bold")) +
    h$legend_theme()
  h$legend_row(h$fix_frame(gg$ggplotGrob(p), xlab_height = grid::unit(5, "pt")))
}

# A two-panel patchwork after mco-snodas plot_hypsomeswe(): USDM classes
# stacked by elevation with this year's line, and an area panel; legends
# collected. Exercises plot_fix_frame() on a patchworkGrob ("-1" names).
patchwork_plot <- function(h) {
  gg <- asNamespace("ggplot2")
  b <- h$bands(TRUE)
  bins <- seq(4000, 9000, by = 250)
  q <- outer(pmax(0, (bins - 4000) / 5000), c(0.2, 0.3, 0.45, 0.6, 0.75, 0.88, 1.12, 1.25,
                                              1.4, 1.55, 1.7, 1.8)) * 1000
  cls <- names(b$bands)
  stack <- do.call(rbind, lapply(seq_along(cls), function(i)
    data.frame(bin = bins, class = cls[i], xmin = q[, i], xmax = q[, i + 1])))
  stack <- stack[stack$xmax > stack$xmin, ]
  med <- data.frame(bin = bins, x = rowMeans(q[, 6:7]))
  obs <- data.frame(bin = bins, x = med$x * (0.7 + 0.2 * sin(bins / 700)))

  left <- gg$ggplot() +
    gg$geom_rect(data = stack, mapping = gg$aes(xmin = xmin, xmax = xmax, ymin = bin - 100,
                                                ymax = bin + 100, fill = class),
                 alpha = 0.7, colour = "white", linewidth = 0.15, key_glyph = "rect",
                 show.legend = c(fill = TRUE, linetype = FALSE)) +
    gg$geom_path(data = med, mapping = gg$aes(x = x, y = bin, linetype = "Median"), colour = "gray50") +
    gg$geom_path(data = obs, mapping = gg$aes(x = x, y = bin), colour = "white", linewidth = 1.8) +
    gg$geom_path(data = obs, mapping = gg$aes(x = x, y = bin, linetype = "Observed"), linewidth = 1) +
    gg$scale_fill_manual(name = "Historical Conditions WY2004–2025",
                         values = b$colors, limits = names(b$bands),
                         breaks = c("D4", "D3", "D2", "D1", "D0", "W4", "W3", "W2", "W1", "W0"),
                         guide = gg$guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                                 title.position = "top",
                                                 override.aes = list(alpha = 0.7))) +
    h$current_legend(1) +
    gg$scale_y_continuous(name = "Elevation (feet)", expand = gg$expansion(0, 0)) +
    gg$scale_x_continuous(name = "Snow Water Equivalent (acre-feet)",
                          expand = gg$expansion(c(0, 0.05), 0)) +
    gg$theme(plot.margin = gg$margin(0, 22, 0, 2.75))
  right <- gg$ggplot(data.frame(bin = bins, area = 1 + sin(bins / 900)^2)) +
    gg$geom_rect(gg$aes(xmin = 0, xmax = area, ymin = bin - 100, ymax = bin + 100), fill = "gray45") +
    gg$scale_x_continuous("Area (sq. mi.)", expand = gg$expansion(c(0, 0.05), 0)) +
    gg$scale_y_continuous(NULL, expand = gg$expansion(0, 0)) +
    gg$theme(axis.text.y = gg$element_blank(), axis.ticks.y = gg$element_blank(),
             plot.margin = gg$margin(0, 0, 0, 2.75))

  p <- (left + right +
          patchwork::plot_layout(widths = c(3, 1), guides = "collect") +
          patchwork::plot_annotation(
            title = "HypsomeSWE for Somewhere Basin (HUC: 10020001)",
            subtitle = paste0("February 15, 2026 — ", h$conditions(5, 10, 20.4)),
            caption = "Fitted History",
            theme = gg$theme(plot.title = gg$element_text(size = h$title_size),
                             plot.caption = gg$element_text(hjust = 0, color = "gray40"),
                             plot.margin = gg$margin(5.5, 5.5, 5.5, 2.75)))) &
    gg$theme(text = gg$element_text(family = "Helvetica"),
             panel.background = gg$element_rect(fill = "white", color = NA),
             panel.border = gg$element_rect(fill = NA, color = "gray20"),
             panel.grid.major = gg$element_line(color = "gray92", linewidth = 0.3),
             panel.grid.minor = gg$element_blank(),
             plot.background = gg$element_rect(fill = "white", color = NA),
             legend.key = gg$element_rect(fill = "white", color = NA),
             axis.ticks = gg$element_line(color = "gray20"),
             plot.title = gg$element_text(face = "bold"),
             plot.subtitle = gg$element_text(face = "bold"),
             legend.position = "bottom",
             axis.title = gg$element_text(face = "bold")) &
    h$legend_theme()
  h$legend_row(h$fix_frame(patchwork::patchworkGrob(p)))
}

# Draw `build` with both helper sets and return the two decoded PNGs.
render_both <- function(build) in_reference_dir({
  out <- list()
  for (which in c("reference", "package")) {
    h <- helper_set(which)
    f <- paste0(which, ".png")
    h$save_png(build(h), f)
    out[[which]] <- png::readPNG(f)
  }
  out
})

test_that("seasonal plot, percentile classes: pixel-identical to mco-snodas", {
  skip_if_no_plot_pkgs()
  suppressWarnings(register_plot_fonts())
  x <- synthetic_year()
  img <- render_both(function(h) seasonal_plot(h, x, usdm = FALSE))
  expect_identical(dim(img$package)[1:2], c(1000L, 1500L))
  expect_identical(img$package, img$reference)
})

test_that("seasonal plot, USDM classes with the outlook key: pixel-identical to mco-snotel", {
  skip_if_no_plot_pkgs()
  suppressWarnings(register_plot_fonts())
  x <- synthetic_year()
  rest <- x$date > as.Date("2026-02-15") & x$date <= as.Date("2026-07-15")
  outlook <- data.frame(date = x$date[rest], q0100 = x$q0100[rest] * 0.9,
                        q0900 = x$q0900[rest] * 1.1)
  img <- render_both(function(h) seasonal_plot(h, x, usdm = TRUE, outlook = outlook))
  expect_identical(img$package, img$reference)
})

test_that("patchwork plot: pixel-identical to mco-snodas", {
  skip_if_no_plot_pkgs()
  suppressWarnings(register_plot_fonts())
  img <- render_both(patchwork_plot)
  expect_identical(img$package, img$reference)
})

test_that("non-drawing helpers match the reference", {
  skip_if_no_plot_pkgs()
  r <- reference_env()
  expect_identical(plot_bands(TRUE), r$snodas_bands(TRUE))
  expect_identical(plot_bands(FALSE), r$snodas_bands(FALSE))
  expect_identical(plot_month_labels(as.Date("2026-01-01")), r$snodas_month_labels(as.Date("2026-01-01")))
  starts <- as.Date("2025-10-01") + c(0, 1, 14, 30, 31, 45, 92)
  for (s in as.list(starts)) {
    lim <- c(s, as.Date("2026-09-30"))
    expect_identical(plot_month_breaks(lim), r$snodas_month_breaks(lim))
  }
  in_reference_dir({
    expect_identical(plot_palette(), r$snodas_palette())
    expect_identical(plot_logo("mco"), r$snodas_logo())
    expect_identical(plot_logo("umt"), r$snodas_umt_logo())
  })
  cases <- list(c(0, 0, 50), c(5, 10, 20.4), c(5, 0, 70), c(5, NA, NA), c(0, 10, NA),
                c(5, 10, NA), c(12, 10, 61.5), c(NA, 10, 3))
  for (cs in cases)
    expect_identical(plot_conditions(cs[1], cs[2], cs[3]), r$snodas_conditions(cs[1], cs[2], cs[3]))
  expect_identical(PLOT_TITLE_SIZE, r$SNODAS_TITLE_SIZE)
  expect_identical(PLOT_UMT_LOGO_HEIGHT, r$UMT_LOGO_HEIGHT)
})
