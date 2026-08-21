test_that("AllNA returns TRUE for all-NA vector", {
  expect_true(nmw:::AllNA(c(NA, NA, NA)))
})

test_that("AllNA returns FALSE for mixed vector", {
  expect_false(nmw:::AllNA(c(1, NA, 3)))
})

test_that("AllSame returns TRUE for identical values", {
  expect_true(nmw:::AllSame(c(5, 5, 5)))
})

test_that("AllSame returns FALSE for different values", {
  expect_false(nmw:::AllSame(c(1, 2, 3)))
})

test_that("AllSame handles single element", {
  expect_true(nmw:::AllSame(c(42)))
})

test_that("FuncDep detects functional dependency", {
  df <- data.frame(
    ID = c(1, 1, 2, 2),
    SEX = c(0, 0, 1, 1),
    TIME = c(0, 1, 0, 1)
  )
  expect_true(nmw:::FuncDep(df, "ID", "SEX"))
  expect_false(nmw:::FuncDep(df, "ID", "TIME"))
})

test_that("NMVarStat returns correct structure", {
  df <- data.frame(
    ID = c(1L, 1L, 2L, 2L),
    TIME = c(0, 1, 0, 1),
    DV = c(10.5, 8.3, 12.1, 9.7),
    MDV = c(0L, 0L, 0L, 0L)
  )
  result <- NMVarStat(df)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)  # 4 variables
  expect_true("nNA" %in% colnames(result))
  expect_true("nUniq" %in% colnames(result))
  expect_true("ALLInt" %in% colnames(result))
})

test_that("NMVarStat counts NA correctly", {
  df <- data.frame(
    ID = c(1L, 1L, 2L, 2L),
    TIME = c(0, NA, 0, 1),
    DV = c(NA, NA, 12.1, 9.7),
    MDV = c(0L, 0L, 0L, 0L)
  )
  result <- NMVarStat(df)
  expect_equal(result["TIME", "nNA"], 1)
  expect_equal(result["DV", "nNA"], 2)
})

test_that("NMIDStat returns correct structure", {
  df <- data.frame(
    ID = c(1, 1, 2, 2, 2),
    TIME = c(0, 1, 0, 1, 2),
    DV = c(0, 10, 0, 12, 8),
    MDV = c(1, 0, 1, 0, 0),
    AMT = c(100, 0, 200, 0, 0),
    EVID = c(1, 0, 1, 0, 0)
  )
  result <- NMIDStat(df)
  expect_true(is.list(result))
  expect_equal(length(result), 3)
  expect_true(result[[2]])  # sorted by ID
  expect_true(result[[3]])  # sorted by time
  expect_equal(result[[1]]["1", "nRec"], 2)
  expect_equal(result[[1]]["2", "nRec"], 3)
})

test_that("run.test.nm returns p-value in [0, 1]", {
  res <- c(1, -1, 1, -1, 1, -1)
  p <- run.test.nm(res)
  expect_true(p >= 0 && p <= 1)
})

test_that("run.test.nm detects non-random pattern", {
  res <- c(1, 2, 3, 4, 5, -1, -2, -3, -4, -5)
  p <- run.test.nm(res)
  expect_true(p < 0.05)
})

test_that("run.p returns 1 for m=0, r=1", {
  expect_equal(nmw:::run.p(0, 5, 1), 1)
})

test_that("run.p returns 0 for invalid inputs", {
  expect_equal(nmw:::run.p(5, 3, 2), 0)  # m > n
})
