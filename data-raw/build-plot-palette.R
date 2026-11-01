# Build the fixed 256-colour palette every plot is written with
# (inst/plot/plot_palette.png; plot_save_png() in R/plot.R). From mco-snodas
# scripts/build-plot-palette.R (62fe19a).
#
# One palette for all plots means a colour is the same in every plot, so the
# shared parts (legends, logos, axes) are pixel-identical across plots and
# frames. Rebuild it when plot colours or styling change, from the package
# root:
#
#   Rscript data-raw/build-plot-palette.R SAMPLE_DIR
#
# SAMPLE_DIR holds full-colour sample plots (not on the palette they build:
# ggsave() with ragg::agg_png, not plot_save_png()), all the same size, e.g.
# a pipeline's plots of a dozen basins from mountain to plains on five dates
# across the season, every flavor, every plot type. mco-snodas renders such a
# sample with its stage_5_viz(); its SAMPLE_ZONES and SAMPLE_DATES were
#   17010205 17010305 10030203 10 1701 1002 14010001 05130106 18040006
#   13020101 11 04010201, on 2025-11-15 2026-01-15 2026-03-15 2026-05-01
#   2026-06-15.
# The script tiles 36 of them into one image, has pngquant (needed only here)
# choose 256 colours for it, and writes them as a 16 x 16 PNG.

N_TILES <- 36

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: Rscript data-raw/build-plot-palette.R SAMPLE_DIR")
if (!nzchar(Sys.which("pngquant"))) stop("pngquant is needed to build the palette")
sample_dir <- args[1]
out <- file.path("inst", "plot", "plot_palette.png")
if (!dir.exists(dirname(out))) stop("Run from the package root")

files <- list.files(sample_dir, "\\.png$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No PNGs in ", sample_dir)
set.seed(1)
files <- if (length(files) > N_TILES) sample(files, N_TILES) else files

imgs <- lapply(files, \(f) png::readPNG(f)[, , 1:3])
h <- dim(imgs[[1]])[1]
w <- dim(imgs[[1]])[2]
if (!all(vapply(imgs, \(x) identical(dim(x)[1:2], c(h, w)), logical(1))))
  stop("Sample plots must all be the same size")
n_col <- ceiling(sqrt(length(imgs)))
n_row <- ceiling(length(imgs) / n_col)
tiles <- array(1, c(n_row * h, n_col * w, 3))   # white where the grid is short
for (i in seq_along(imgs)) {
  r <- (i - 1) %/% n_col
  c <- (i - 1) %% n_col
  tiles[r * h + seq_len(h), c * w + seq_len(w), ] <- imgs[[i]]
}
rm(imgs)

montage <- tempfile(fileext = ".png")
quant   <- tempfile(fileext = ".png")
png::writePNG(tiles, montage)
rm(tiles)
status <- system2("pngquant", c("256", "--nofs", "--speed", "1", "--force", "--output", shQuote(quant), shQuote(montage)))
if (status != 0) stop("pngquant failed (exit ", status, ")")

q <- png::readPNG(quant)
cols <- unique(round(cbind(c(q[, , 1]), c(q[, , 2]), c(q[, , 3])) * 255))
cols <- cols[order(cols[, 1], cols[, 2], cols[, 3]), , drop = FALSE]
pal <- array(1, c(16, 16, 3))
for (k in seq_len(nrow(cols))) pal[(k - 1) %/% 16 + 1, (k - 1) %% 16 + 1, ] <- cols[k, ] / 255
png::writePNG(pal, out)
unlink(c(montage, quant))
message(nrow(cols), " colours from ", length(files), " plots written to ", out)
