# The 50 x 50 pixel fixture (API.md section 7, Q5): a checksummed GitHub
# release asset, made by `Rscript data-raw/make-fixtures.R --large`. Fetched
# when online; skipped otherwise. NORMALS_GRID50 may point to a local copy.

grid50_asset <- function() {
  url <- "https://github.com/mt-climate-office/normals/releases/download/fixtures-v1/grid50.rds"
  sha <- readRDS(test_path("fixtures", "grid50-sha256.rds"))
  path <- Sys.getenv("NORMALS_GRID50")
  if (!nzchar(path)) {
    skip_on_cran()
    skip_if_not_installed("curl")
    skip_if_offline("github.com")
    path <- tempfile(fileext = ".rds")
    ok <- tryCatch(utils::download.file(url, path, mode = "wb", quiet = TRUE) == 0,
                   error = \(e) FALSE, warning = \(w) FALSE)
    if (!ok) skip("the grid50.rds release asset is not published")
  }
  if (golden_sha256(path) != sha) fail("grid50.rds does not match its sha256 in fixtures/grid50-sha256.rds")
  readRDS(path)
}

test_that("gridded normals match mco-snodas (50 x 50 release asset)", {
  fx <- grid50_asset()
  grid_golden_check(fx)
  grid_golden_check_smoothed(fx)
})
