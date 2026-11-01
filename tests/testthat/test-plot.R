test_that("plot_bands() has no default and returns matching classes", {
  expect_error(plot_bands())
  for (usdm in c(TRUE, FALSE)) {
    b <- plot_bands(usdm)
    expect_named(b, c("bands", "breaks", "colors", "legend", "ncol"))
    expect_identical(names(b$colors), names(b$bands))
    expect_length(b$breaks, length(b$bands) + 1L)
    expect_true(all(unlist(b$bands) %in% names(NORMALS_QUANTILES)))
    expect_true(all(b$legend %in% names(b$bands)))
  }
  expect_identical(plot_bands(FALSE)$legend[1], "1st–5th")
})

test_that("plot_month_breaks() puts a break on every first of the month in the limits", {
  lim <- as.Date(c("2025-10-15", "2026-01-31"))
  expect_identical(plot_month_breaks(lim), as.Date(c("2025-11-01", "2025-12-01", "2026-01-01")))
  lim <- as.Date(c("2025-10-01", "2026-09-30"))
  expect_identical(plot_month_breaks(lim), seq(as.Date("2025-10-01"), by = "month", length.out = 12))
  expect_identical(plot_month_labels(as.Date("2026-02-01")), "Feb\n2026")
})

test_that("plot_conditions() reads as the subtitle does", {
  skip_if_not_installed("scales")
  expect_identical(plot_conditions(5, 10, 20.4), "50% of Median (~20th Percentile)")
  expect_identical(plot_conditions(0, 0, 50), "No snow (typical for this date)")
  expect_identical(plot_conditions(5, 0, 71), "~71st Percentile")
  expect_identical(plot_conditions(5, 10, NA), "50% of Median")
  expect_identical(plot_conditions(5, NA, NA), "No normal for this date")
})

test_that("logos and palette load from the package", {
  skip_if_not_installed("png")
  mco <- plot_logo("mco")
  expect_identical(dim(mco)[c(1, 3)], c(400L, 4L))
  expect_identical(dim(plot_logo("umt"))[3] %in% 3:4, TRUE)
  expect_error(plot_logo("nope"))
  pal <- plot_palette()
  expect_identical(dim(pal$rgb), c(256L, 3L))
  expect_length(pal$hex, 256L)
  expect_identical(plot_logo("mco"), mco)   # cached
})

test_that("plot_logo_grob() resamples itself to the device", {
  skip_if_not_installed("png")
  skip_if_not_installed("ragg")
  g <- plot_logo_grob(plot_logo("umt"), height = 0.3)
  expect_s3_class(g, "normals_logo")
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  ragg::agg_png(f, width = 2, height = 1, units = "in", res = 100)
  drawn <- grid::makeContent(g)
  grDevices::dev.off()
  r <- drawn$children[[1]]
  expect_s3_class(r, "rastergrob")
  expect_identical(nrow(r$raster), 30L)
})

test_that("plot_quantize() maps colors to their nearest palette index", {
  skip_if_not_installed("png")
  pal <- plot_palette()$rgb
  # nativeRaster is 0xAABBGGRR; two exact palette colors and white-ish.
  col <- function(r, g, b) bitwOr(bitwOr(r, bitwShiftL(g, 8L)), bitwShiftL(b, 16L)) - 16777216L
  i <- c(1L, 100L)
  nr <- structure(c(col(pal[i[1], 1], pal[i[1], 2], pal[i[1], 3]),
                    col(pal[i[2], 1], pal[i[2], 2], pal[i[2], 3])),
                  dim = c(1L, 2L), class = "nativeRaster")
  expect_identical(plot_quantize(nr), matrix(i - 1L, nrow = 1))
})

test_that("plot_save_png() writes an 8-bit palette PNG at the shared size", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ragg")
  skip_if_not_installed("fastpng")
  skip_if_not_installed("png")
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_identical(plot_save_png(p, f), f)
  img <- png::readPNG(f)
  expect_identical(dim(img)[1:2], c(1000L, 1500L))
  hex <- unique(grDevices::rgb(c(img[, , 1]), c(img[, , 2]), c(img[, , 3])))
  expect_true(all(hex %in% plot_palette()$hex))
})

test_that("plot_current_legend() keys Observed and Median, optionally the likely range", {
  skip_if_not_installed("ggplot2")
  s <- plot_current_legend(1)
  expect_identical(s$get_limits(), c("Observed", "Median"))
  s <- plot_current_legend(1, likely_range = TRUE)
  expect_identical(s$get_limits(), c("Observed", "Median", "Likely range"))
})

test_that("the makeContent() method is registered", {
  expect_true(is.function(utils::getS3method("makeContent", "normals_logo", envir = asNamespace("grid"))))
})
