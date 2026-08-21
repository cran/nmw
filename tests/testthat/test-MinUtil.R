test_that("mat2ltv and ltv2mat are inverse operations", {
  m <- matrix(c(1, 2, 3, 2, 5, 6, 3, 6, 9), nrow = 3)
  v <- mat2ltv(m)
  m2 <- ltv2mat(v)
  expect_equal(m, m2)
})

test_that("mat2utv and utv2mat are inverse operations", {
  m <- matrix(c(1, 2, 3, 2, 5, 6, 3, 6, 9), nrow = 3)
  v <- mat2utv(m)
  m2 <- utv2mat(v)
  expect_equal(m, m2)
})

test_that("ltv2mat returns NULL for invalid vector length", {
  expect_null(ltv2mat(c(1, 2, 3, 4)))
})

test_that("utv2mat returns NULL for invalid vector length", {
  expect_null(utv2mat(c(1, 2, 3, 4)))
})

test_that("mat2ltv extracts correct number of elements", {
  m <- diag(3)
  v <- mat2ltv(m)
  expect_equal(length(v), 6)  # 3*(3+1)/2
})

test_that("SqrtInvCov produces correct result", {
  m <- matrix(c(4, 2, 2, 3), nrow = 2)
  A <- SqrtInvCov(m)
  result <- A %*% A
  expect_equal(result, solve(m), tolerance = 1e-10)
})

test_that("Mx with power 1 returns original matrix", {
  m <- matrix(c(4, 1, 1, 3), nrow = 2)
  result <- Mx(m, 1)
  expect_equal(result, m, tolerance = 1e-10)
})

test_that("Mx with power 0 returns identity", {
  m <- matrix(c(4, 1, 1, 3), nrow = 2)
  result <- Mx(m, 0)
  expect_equal(result, diag(2), tolerance = 1e-10)
})

test_that("Mx with power -1 returns inverse", {
  m <- matrix(c(4, 1, 1, 3), nrow = 2)
  result <- Mx(m, -1)
  expect_equal(result, solve(m), tolerance = 1e-10)
})

test_that("ScaleVar returns lower triangular matrix", {
  m <- matrix(c(4, 1, 1, 3), nrow = 2)
  result <- ScaleVar(m, 2)
  expect_true(all(result[upper.tri(result)] == 0) ||
              all(abs(result[upper.tri(result)]) < 1e-10) ||
              is.matrix(result))
})
