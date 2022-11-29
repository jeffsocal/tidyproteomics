test_that("invlog2() returns the inverse log2 of a value", {
  y <- 12345
  x <- log2(y)
  expect_equal(invlog2(x), y)
})
