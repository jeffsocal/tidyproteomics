test_that("hdf() returns the top portion of a tibble as a data.frame", {
  x <- tibble::tibble(a = 1:10, b = 11:20)
  y <- data.frame(a = 1:5, b = 11:15)
  expect_equal(hdf(x), y)
})
