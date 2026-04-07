test_that("CountEXTParams counts parameters correctly", {
  EXT <- data.frame(
    ITERATION = 1,
    THETA1 = 0.1, THETA2 = 0.2, THETA3 = 0.3,
    OMEGA.1.1. = 0.04, OMEGA.2.1. = 0.01, OMEGA.2.2. = 0.09,
    SIGMA.1.1. = 0.01,
    OBJ = 100
  )
  result <- nmw:::CountEXTParams(EXT)
  expect_equal(result$nTheta, 3)
  expect_equal(result$nEta, 2)   # 3 omega elements -> 2x2 matrix
  expect_equal(result$nEps, 1)   # 1 sigma element -> 1x1 matrix
})

test_that("CountEXTParams handles single parameter", {
  EXT <- data.frame(
    ITERATION = 1,
    THETA1 = 0.5,
    OMEGA.1.1. = 0.1,
    SIGMA.1.1. = 0.01,
    OBJ = 50
  )
  result <- nmw:::CountEXTParams(EXT)
  expect_equal(result$nTheta, 1)
  expect_equal(result$nEta, 1)
  expect_equal(result$nEps, 1)
})

test_that("CalcTaLDForReport calculates time after latest dose", {
  FDATA <- data.frame(
    ID = c(1, 1, 1, 2, 2),
    TIME = c(0, 1, 4, 0, 2),
    AMT = c(100, 0, 0, 200, 0),
    MDV = c(1, 0, 0, 1, 0)
  )
  result <- nmw:::CalcTaLDForReport(FDATA, nrow(FDATA))
  expect_equal(length(result$TaLD), 5)
  expect_equal(result$TaLD[1], 0)  # dosing record
  expect_equal(result$TaLD[2], 1)  # 1h after dose at TIME=0
  expect_equal(result$TaLD[3], 4)  # 4h after dose at TIME=0
  expect_equal(result$TaLD[4], 0)  # dosing record for ID=2
  expect_equal(result$TaLD[5], 2)  # 2h after dose at TIME=0
})

test_that("CalcTaLDForReport handles multiple dosing occasions", {
  FDATA <- data.frame(
    ID = c(1, 1, 1, 1),
    TIME = c(0, 1, 24, 25),
    AMT = c(100, 0, 100, 0),
    MDV = c(1, 0, 1, 0)
  )
  result <- nmw:::CalcTaLDForReport(FDATA, nrow(FDATA))
  expect_equal(result$TaLD[2], 1)   # 1h after first dose
  expect_equal(result$TaLD[4], 1)   # 1h after second dose
})

test_that("OFV_SCREEN_LAYOUT has correct dimensions", {
  layout <- nmw:::OFV_SCREEN_LAYOUT
  expect_equal(dim(layout), c(4, 4))
  expect_equal(layout[1, ], c(0, 1, 0, 1))  # full screen for first panel
})
