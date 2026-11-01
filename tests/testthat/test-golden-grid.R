# mco-snodas gridded normals (calc_snodas_normals(), annual_fit_block(),
# annual_materialize()): API.md 6.2, on the synthetic pixel block.

test_that("gridded ECDF, ZIG and ZIG-smooth normals and percentiles match mco-snodas (20 x 20)", {
  grid_golden_check(golden("grid"))
})
