test_that("NORMALS_QUANTILES is the 19 probabilities, named q%04d", {
  p <- c(0, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5,
         0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.98, 0.99, 1)
  expect_identical(unname(NORMALS_QUANTILES), p)
  expect_identical(names(NORMALS_QUANTILES), sprintf("q%04d", round(p * 1000)))
  expect_identical(names(NORMALS_QUANTILES)[c(1, 4, 10, 19)],
                   c("q0000", "q0025", "q0500", "q1000"))
})

test_that("flavor vocabularies agree", {
  expect_identical(names(NORMALS_FLAVORS), c("empirical", "fitted", "smoothed"))
  expect_identical(unname(NORMALS_FLAVORS), c("ecdf", "zig", "zig_smooth"))
  expect_identical(names(NORMALS_FLAVOR_CAPTIONS), names(NORMALS_FLAVORS))
})

test_that("stubs stop until the API is approved", {
  expect_error(normals_table(), "not implemented")
  expect_error(family_zig(), "not implemented")
})
