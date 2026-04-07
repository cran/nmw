test_that("ResTest returns correct structure", {
  set.seed(42)
  ResTab <- data.frame(
    WRES = rnorm(100),
    CWRE = rnorm(100),
    IWRE = rnorm(100)
  )
  result <- ResTest(ResTab)
  expect_true(is.matrix(result))
  expect_true("Z-value" %in% colnames(result))
  expect_true("Expected" %in% colnames(result))
  expect_true("WRES Cnt" %in% colnames(result))
  expect_true("CWRE Cnt" %in% colnames(result))
  expect_true("IWRE Cnt" %in% colnames(result))
})

test_that("ResTest handles custom columns", {
  set.seed(42)
  ResTab <- data.frame(
    WRES = rnorm(50),
    CWRE = rnorm(50),
    IWRE = rnorm(50)
  )
  result <- ResTest(ResTab, TestCols = c("WRES"))
  expect_true("WRES Cnt" %in% colnames(result))
  expect_false("CWRE Cnt" %in% colnames(result))
})

test_that("mlr2 returns correct structure", {
  set.seed(42)
  y <- rnorm(20)
  x <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  result <- mlr2(y, x)
  expect_true(is.list(result))
  expect_equal(length(result), 7)
  expect_equal(result[[3]], 20)  # n
  expect_equal(result[[4]], 3)   # p (intercept + 2 predictors)
})

test_that("mlr2 with standardization works", {
  set.seed(42)
  y <- rnorm(20)
  x <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  result <- mlr2(y, x, standardize = 3)
  expect_true(is.list(result))
  expect_true(grepl("Standardization", names(result)[1]))
})

test_that("mlr2 handles mismatched lengths", {
  y <- rnorm(10)
  x <- data.frame(x1 = rnorm(5), x2 = rnorm(5))
  expect_message(result <- mlr2(y, x), "different")
  expect_null(result)
})
