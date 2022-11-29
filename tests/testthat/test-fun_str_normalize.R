test_that("str_normalize() returns SQL like column names", {
  x <- tibble::tibble(`Col A` = 1:10, `Col B` = 11:20)
  y <- c('col_a', 'col_b')
  z <- str_normalize(colnames(x))
  expect_equal(z, y)
})
