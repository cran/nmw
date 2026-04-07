test_that("RemoveNA removes rows with NA", {
  df <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  result <- RemoveNA(df)
  expect_equal(nrow(result), 1)
  expect_equal(result$a, 1)
  expect_equal(result$b, 4)
})

test_that("RemoveNA keeps all rows when no NA", {
  df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  result <- RemoveNA(df)
  expect_equal(nrow(result), 3)
})

test_that("RmvFixed removes constant columns", {
  m <- data.frame(a = c(1, 1, 1), b = c(1, 2, 3), c = c(5, 5, 5))
  result <- nmw:::RmvFixed(m)
  # When only one column remains, R drops to a vector
  expect_equal(as.numeric(result), c(1, 2, 3))
})

test_that("RmvFixed keeps multiple varying columns", {
  m <- data.frame(a = c(1, 1, 1), b = c(1, 2, 3), c = c(5, 6, 7))
  result <- nmw:::RmvFixed(m)
  expect_equal(ncol(result), 2)
  expect_true("b" %in% colnames(result))
  expect_true("c" %in% colnames(result))
})

test_that("RmvZero removes zero rows/cols from symmetric matrix", {
  m <- matrix(c(1, 0, 0.5, 0, 0, 0, 0.5, 0, 2), nrow = 3)
  result <- nmw:::RmvZero(m)
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("RmvZero handles all-nonzero matrix", {
  m <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  result <- nmw:::RmvZero(m)
  expect_equal(dim(result), c(2, 2))
})

test_that("RmvCol removes specified columns", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  result <- nmw:::RmvCol(df, c("b"))
  expect_equal(ncol(result), 2)
  expect_false("b" %in% colnames(result))
})

test_that("RenCol renames specified column", {
  cols <- c("a", "b", "c")
  result <- nmw:::RenCol(cols, "b", "B_new")
  expect_equal(result, c("a", "B_new", "c"))
})

test_that("ExpandDoseHist expands ADDL records", {
  dose_tab <- data.frame(
    TIME = c(0),
    AMT = c(100),
    II = c(24),
    ADDL = c(2)
  )
  result <- ExpandDoseHist(dose_tab)
  expect_equal(nrow(result), 3)  # original + 2 ADDL
  expect_equal(sort(result$TIME), c(0, 24, 48))
})

test_that("MatchEnd matches string endings", {
  names <- c("file.csv", "data.txt", "report.csv", "notes.doc")
  result <- MatchEnd(names, ".csv")
  expect_equal(result, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("MatchEnd is case insensitive by default", {
  names <- c("file.CSV", "data.txt")
  result <- MatchEnd(names, ".csv")
  expect_equal(result, c(TRUE, FALSE))
})

test_that("MatchEnd respects case when IgnoreCase=FALSE", {
  names <- c("file.CSV", "data.csv")
  result <- MatchEnd(names, ".csv", IgnoreCase = FALSE)
  expect_equal(result, c(FALSE, TRUE))
})
