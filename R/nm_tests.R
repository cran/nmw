#' Run Test Probability
#'
#' Calculates P(Run count <= r | m, n) for the runs test.
#'
#' @param m count of fewer species (minimum value = 0)
#' @param n count of more frequent species (minimum value = 1)
#' @param r count of runs (minimum value = 1)
#' @return numeric, probability
#' @keywords internal
run.p <- function(m, n, r) {
  if (m == 0 & r == 1) {
    return(1)
  } else if (m > n | m < 1 | n < 1 | r < 2 | (r > min(m + n, 2 * m + 1))) {
    return(0)
  }

  sumfu <- 0
  for (u in 2:r) {
    if (u %% 2 == 0) {
      k <- u / 2
      fu <- 2 * choose(m - 1, k - 1) * choose(n - 1, k - 1)
    } else {
      k <- (u + 1) / 2
      fu <- choose(m - 1, k - 1) * choose(n - 1, k - 2) + choose(m - 1, k - 2) * choose(n - 1, k - 1)
    }
    sumfu <- sumfu + fu
  }
  return(sumfu / choose(m + n, m))
}


#' Runs Test for NONMEM Residuals
#'
#' Performs a runs test on residuals. Zeros are omitted.
#'
#' @param RES numeric vector of residuals
#' @return numeric, two-sided p-value
#' @export
run.test.nm <- function(RES) {
  RES <- RES[RES != 0]
  t <- length(RES)
  r <- RES > 0
  m <- sum(r)
  if (t > 1) {
    j <- 2:t
    run <- sum(abs(r[j] - r[j - 1])) + 1
  } else {
    run <- 1
  }
  m <- min(m, t - m)
  n <- max(m, t - m)
  if (run > 1) {
    p <- run.p(m, n, run)
    if (p > 0.5) {
      p <- 1 - run.p(m, n, run - 1)
    }
  } else {
    p <- 0.5^(t - 1)
  }
  return(p)
}


#' Residual Count Test Using Binomial Distribution
#'
#' Tests whether residual counts at various Z-value thresholds follow
#' expected binomial distributions.
#'
#' @param ResTab data.frame with residual columns
#' @param TestCols character vector of column names to test
#' @param ZVals numeric vector of Z-value thresholds
#' @param PctVals numeric vector of percentile values
#' @return matrix of test results
#' @export
ResTest <- function(ResTab, TestCols = c("WRES", "CWRE", "IWRE"),
                    ZVals = c(0, 1, 2, 3),
                    PctVals = c(0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) {
  nPct <- length(PctVals)
  for (i in 1:nPct) {
    ZVals <- union(ZVals, abs(qnorm(PctVals[i])))
  }
  ZVals <- sort(ZVals)
  nZVals <- length(ZVals)

  nRec <- dim(ResTab)[1]
  nTestCols <- length(TestCols)
  CountCols <- c("Z-value", "Percent", "Expected", "LB", "UB")
  for (i in 1:nTestCols) {
    CountCols <- c(CountCols, paste(TestCols[i], "Cnt"), paste(TestCols[i], "p-val"))
  }

  Count <- matrix(nrow = 2 * nZVals, ncol = length(CountCols))
  colnames(Count) <- CountCols

  for (i in 1:nZVals) {
    Count[i, "Z-value"] <- round(-1 * ZVals[nZVals + 1 - i], digits = 3)
    Count[i + nZVals, "Z-value"] <- round(ZVals[i], digits = 3)
    NegProb <- pnorm(-1 * ZVals[nZVals + 1 - i])
    PosProb <- 1 - pnorm(ZVals[i])
    Count[i, "Percent"] <- round(NegProb * 100, digits = 2)
    Count[i + nZVals, "Percent"] <- round(PosProb * 100, digits = 2)
    Count[i, "Expected"] <- round(NegProb * nRec, digits = 1)
    Count[i + nZVals, "Expected"] <- round(PosProb * nRec, digits = 1)
    Count[i, "LB"] <- qbinom(0.025, nRec, NegProb)
    Count[i, "UB"] <- qbinom(0.975, nRec, NegProb) - 1
    Count[i + nZVals, "LB"] <- qbinom(0.025, nRec, PosProb)
    Count[i + nZVals, "UB"] <- qbinom(0.975, nRec, PosProb) - 1

    for (j in 1:length(TestCols)) {
      NegCount <- length(ResTab[ResTab[, paste(TestCols[j])] < (-1 * ZVals[nZVals + 1 - i]),
                                paste(TestCols[j])])
      PosCount <- length(ResTab[ResTab[, paste(TestCols[j])] > ZVals[i],
                                paste(TestCols[j])])

      NegPVal <- pbinom(NegCount, nRec, NegProb)
      PosPVal <- pbinom(PosCount, nRec, PosProb)

      Count[i, paste(TestCols[j], "Cnt")] <- NegCount
      Count[i + nZVals, paste(TestCols[j], "Cnt")] <- PosCount

      Count[i, paste(TestCols[j], "p-val")] <- round(min(NegPVal, 1 - NegPVal), digits = 3)
      Count[i + nZVals, paste(TestCols[j], "p-val")] <- round(min(PosPVal, 1 - PosPVal), digits = 3)
    }
  }
  return(Count)
}


#' Multiple Linear Regression with Influence Diagnostics
#'
#' Performs multiple linear regression and computes influence diagnostics
#' including DFFITS, DFBETAS, Cook's Distance, COVRATIO, and R-Student.
#'
#' @param y numeric vector of response
#' @param x.raw data.frame of predictors
#' @param standardize integer, standardization method (0=none, 1=center, 2=scale by mean, 3=z-score)
#' @return list with model estimates and influence diagnostics
#' @export
mlr2 <- function(y, x.raw, standardize = 0) {
  x.raw <- RemoveNA(x.raw)

  if (kappa(x.raw) > 999 & standardize == 0) {
    message(paste("Condition Number is ", kappa(x.raw), ". Consider standardization !", sep = ""))
  }
  if (length(y) != length(x.raw[, 1])) {
    message("Numbers of rows of x matrix and y vector are different.")
    return(NULL)
  }

  n <- length(y)
  namelist <- c("Intercept", names(x.raw))

  x.avg <- matrix(rep(colMeans(x.raw, na.rm = TRUE), each = n), nrow = n)
  if (standardize == 3) {
    x.sd <- matrix(rep(apply(x.raw, 2, sd, na.rm = TRUE), each = n), nrow = n)
    x <- (x.raw - x.avg) / x.sd
  } else if (standardize == 2) {
    x <- x.raw / x.avg
  } else if (standardize == 1) {
    x <- x.raw - x.avg
  } else {
    x <- x.raw
  }

  Intercept <- 1
  x <- as.matrix(cbind(Intercept, x))

  b <- MASS::ginv(t(x) %*% x) %*% t(x) %*% y
  p <- length(b)

  y.hat <- x %*% b
  e <- y - y.hat
  SSE <- sum(e^2)
  MSE <- SSE / (n - p)
  b.se <- sqrt(diag(as.numeric(MSE) * MASS::ginv(t(x) %*% x)))
  b.t <- b / b.se
  b.p <- pt(b.t, n - p)
  for (i in 1:p) {
    if (!is.na(b.p[i]) & b.p[i] > 0.5) b.p[i] <- 1 - b.p[i]
  }
  b.p <- 2 * b.p

  res1 <- data.frame(namelist, cbind(b, b.se, b.t, b.p))
  names(res1) <- c("Variable", "Estimate", "SE", "T", "p-value")

  h <- as.matrix(diag(x %*% MASS::ginv(t(x) %*% x) %*% t(x)))
  sr <- e / sqrt(MSE * (1 - h))
  MSEi <- ((n - p) * MSE - e^2 / (1 - h)) / (n - p - 1)
  sdr <- e / sqrt(MSEi * (1 - h))

  DFFITS <- sqrt(h / (1 - h)) * e / sqrt(MSEi * (1 - h))

  bi <- matrix(nrow = n, ncol = p)
  for (i in 1:n) {
    z <- x[-i, ]
    bi[i, ] <- MASS::ginv(t(z) %*% z) %*% t(z) %*% y[-i]
  }
  bm <- matrix(rep(t(b), n), byrow = TRUE, ncol = p)
  DFBETAS <- (bm - bi) / sqrt(MSEi %*% diag(MASS::ginv(t(x) %*% x)))

  COVRATIO <- matrix(nrow = n)
  for (i in 1:n) {
    COVRATIO[i] <- det(MSEi[i] * MASS::ginv(t(x[-i, ]) %*% x[-i, ])) / det(MSE * MASS::ginv(t(x) %*% x))
  }

  D <- e^2 / (1 - h)^2 * h / (p * MSE)

  res2 <- data.frame(cbind(y.hat, e, sdr, h, D, COVRATIO, DFFITS, DFBETAS))
  names(res2) <- c("Yhat", "Residual", "R-Student", "hat", "Cook's D", "COV-Ratio", "DFFITS", namelist)

  result <- list(res1, res2, n, p, n - p, SSE, MSE)
  if (standardize == 1 | standardize == 2 | standardize == 3) {
    names(result) <- c("Model Estimates with Standardization", "Influence Diagnostics with DFBETAs",
                       "n", "Parameter Count", "Degree of Freedom", "SSE", "MSE")
  } else {
    names(result) <- c("Model Estimates", "Influence Diagnostics with DFBETAs",
                       "n", "Parameter Count", "Degree of Freedom", "SSE", "MSE")
  }

  return(result)
}
