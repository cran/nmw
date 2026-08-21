ctl_path <- function(f) system.file("doc", f, package = "nmw")

# helper: a self-contained control stream as a character vector
mk_ctl <- function(pred, theta = "(0,1)", omega = "0.1", sigma = "0.1",
                    est = "METHOD=COND", input = "$INPUT ID TIME DV",
                    subr = NULL, pk = FALSE) {
  c("$PROB test", input, if (!is.null(subr)) subr else "$DATA d",
    if (pk) c("$PK", pred, "$ERROR", "  Y = F + EPS(1)") else c("$PRED", pred),
    paste("$THETA", theta), paste("$OMEGA", omega), paste("$SIGMA", sigma),
    paste("$EST", est))
}

test_that("ParseCtl reads THETA/OMEGA/SIGMA/EST from the golden FOCEI stream", {
  cs <- nmw:::ParseCtl(ctl_path("THEO-FOCEI.CTL"))
  expect_equal(as.numeric(cs$theta_init), c(2, 50, 0.1))
  expect_equal(cs$lb, c(0, 0, 0))
  expect_equal(cs$ub, rep(1e6, 3))
  expect_equal(cs$om_init,
               matrix(c(0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2), 3))
  expect_equal(cs$sg_init, diag(c(0.1, 0.1)))
  expect_equal(cs$method, "COND")
})

test_that("$EST METHOD maps to ZERO / COND / LAPL", {
  expect_equal(nmw:::ParseCtl(ctl_path("THEO-FO.CTL"))$method, "ZERO")
  expect_equal(nmw:::ParseCtl(ctl_path("THEO-FOCEI.CTL"))$method, "COND")
  expect_equal(nmw:::ParseCtl(ctl_path("THEO-LAPLI.CTL"))$method, "LAPL")
})

test_that("$INPUT drops DROP/SKIP items but keeps ID/TIME/DV", {
  cs <- nmw:::ParseCtl(ctl_path("THEO-FOCEI.CTL"))
  expect_false("AMT" %in% cs$input$used)
  expect_true(all(c("ID", "TIME", "DV") %in% cs$input$used))
})

test_that("ctl2nmw returns an nmwCode object with the expected fields", {
  gen <- ctl2nmw(ctl_path("THEO-FOCEI.CTL"))
  expect_s3_class(gen, "nmwCode")
  expect_type(gen$code, "character")
  expect_null(gen$PRED)
  expect_equal(gen$nTheta, 3L)
  expect_equal(gen$nEta, 3L)
  expect_equal(gen$nEps, 2L)
})

test_that("generated PRED matches the hand-written example PRED (COND and LAPL)", {
  FGDref <- deriv(~DOSE/(TH2*exp(ETA2))*TH1*exp(ETA1)/(TH1*exp(ETA1) - TH3*exp(ETA3))*
                   (exp(-TH3*exp(ETA3)*TIME) - exp(-TH1*exp(ETA1)*TIME)),
                  c("ETA1", "ETA2", "ETA3"),
                  function.arg = c("TH1","TH2","TH3","ETA1","ETA2","ETA3","DOSE","TIME"),
                  func = TRUE, hessian = TRUE)
  Href <- deriv(~F + F*EPS1 + EPS2, c("EPS1", "EPS2"),
                function.arg = c("F", "EPS1", "EPS2"), func = TRUE)
  PREDref <- function(THETA, ETA, DATAi) {
    r <- FGDref(THETA[1], THETA[2], THETA[3], ETA[1], ETA[2], ETA[3], DOSE = 320, DATAi[, "TIME"])
    G <- attr(r, "gradient"); H <- attr(Href(r, 0, 0), "gradient")
    if (e$METHOD == "LAPL") {
      D <- attr(r, "hessian")
      Res <- cbind(r, G, H, D[,1,1], D[,2,1], D[,2,2], D[,3,1], D[,3,2], D[,3,3])
      colnames(Res) <- c("F","G1","G2","G3","H1","H2","D11","D21","D22","D31","D32","D33")
    } else {
      Res <- cbind(r, G, H); colnames(Res) <- c("F","G1","G2","G3","H1","H2")
    }
    Res
  }
  DATAi <- as.matrix(data.frame(TIME = c(0.25, 1, 2, 4, 8, 12), DV = 1))
  for (M in c("COND", "LAPL")) {
    e <- nmw:::e; old <- e$METHOD; on.exit(e$METHOD <- old, add = TRUE); e$METHOD <- M
    gen <- ctl2nmw(ctl_path(if (M == "LAPL") "THEO-LAPLI.CTL" else "THEO-FOCEI.CTL"), eval = TRUE)
    for (th in list(c(1.5, 32, 0.09), c(2, 50, 0.1))) {
      for (et in list(c(0, 0, 0), c(0.1, -0.2, 0.05))) {
        a <- gen$PRED(th, et, DATAi); b <- PREDref(th, et, DATAi)
        expect_equal(colnames(a), colnames(b))
        expect_equal(unname(a), unname(b), tolerance = 1e-8)
      }
    }
  }
})

test_that("combined additive+proportional error wires THETAs that live only in Y", {
  cs <- c("$PROB comb", "$INPUT ID TIME DV", "$DATA d", "$PRED",
          "  V = THETA(2)*EXP(ETA(1))", "  F = THETA(1)/V*EXP(-TIME)",
          "  Y = F + F*EPS(1) + EPS(2)",
          "$THETA (0,10)(0,5)", "$OMEGA 0.1", "$SIGMA 0.1 0.1", "$EST METHOD=COND INTER")
  gen <- ctl2nmw(cs, eval = TRUE)
  expect_equal(gen$nEps, 2L)
  DATAi <- as.matrix(data.frame(TIME = c(1, 2), DV = 1))
  e <- nmw:::e; old <- e$METHOD; on.exit(e$METHOD <- old, add = TRUE); e$METHOD <- "COND"
  M <- gen$PRED(c(10, 5), 0.1, DATAi)
  expect_equal(colnames(M), c("F", "G1", "H1", "H2"))
  expect_equal(as.numeric(M[, "H1"]), as.numeric(M[, "F"]))   # dY/dEPS1 = F
  expect_equal(as.numeric(M[, "H2"]), c(1, 1))                # dY/dEPS2 = 1
})

test_that("out-of-scope constructs raise informative errors", {
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*EXP(ETA(1))\n  Y=F+EPS(1)", theta = "(0,1) FIX")),
               "FIX")
  expect_error(ctl2nmw(mk_ctl("  CL=THETA(1)", subr = "$SUBROUTINE ADVAN2", pk = TRUE)),
               "scope|closed-form")
  expect_error(ctl2nmw(mk_ctl("  IF(TIME.GT.0) F=THETA(1)\n  Y=F+EPS(1)")),
               "conditional")
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*EXP(ETA(1))\n  Y=F+EPS(1)", est = "METHOD=SAEM")),
               "not supported")
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*EXP(ETA(3))\n  Y=F+EPS(1)")),
               "ETA")
})

test_that("sequential reassignment and forward/self references raise informative errors", {
  expect_error(ctl2nmw(mk_ctl("  CL=THETA(1)\n  CL=CL*2\n  F=CL*TIME\n  Y=F+EPS(1)")),
               "assigned more than once")
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*TIME\n  W=W+1\n  Y=F+W*EPS(1)")),
               "before it is defined")
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*TIME\n  A=B+1\n  B=A+1\n  Y=F+A*EPS(1)")),
               "before it is defined")
})

test_that("unbounded non-positive THETA init widens LB to -1e6 with a note", {
  gen <- ctl2nmw(mk_ctl("  F=THETA(1)+THETA(2)*TIME\n  Y=F+EPS(1)",
                        theta = "-0.5 (0,2)"))
  expect_equal(gen$THETAinit, c(-0.5, 2))
  expect_equal(gen$LB, c(-1e6, 0))
  expect_equal(gen$UB, c(1e6, 1e6))
  expect_true(any(grepl("non-positive", gen$notes)))
})

test_that("an initial estimate on or outside explicit $THETA bounds errors informatively", {
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*TIME\n  Y=F+EPS(1)", theta = "(1, 0.5, 2)")),
               "LB < init < UB")
  expect_error(ctl2nmw(mk_ctl("  F=THETA(1)*TIME\n  Y=F+EPS(1)", theta = "(0, 0, 2)")),
               "LB < init < UB")
})

test_that("diagonal $OMEGA is densified to a full block with a note", {
  gen <- ctl2nmw(mk_ctl("  F=THETA(1)*EXP(ETA(1))+THETA(2)*EXP(ETA(2))\n  Y=F+F*EPS(1)",
                        theta = "(0,1)(0,2)", omega = "0.1 0.2"))
  expect_equal(gen$OMinit, matrix(c(0.1, 0, 0, 0.2), 2))
  expect_true(any(grepl("FULL block", gen$notes)))
})

test_that("FO round-trip reproduces the NONMEM golden OFV", {
  skip_on_cran()
  DataAll <- Theoph
  colnames(DataAll) <- c("ID", "BWT", "DOSE", "TIME", "DV")
  DataAll[, "ID"] <- as.numeric(as.character(DataAll[, "ID"]))
  gen <- ctl2nmw(ctl_path("THEO-FO.CTL"), eval = TRUE)
  InitStep(DataAll, THETAinit = gen$THETAinit, OMinit = gen$OMinit, SGinit = gen$SGinit,
           LB = gen$LB, UB = gen$UB, Pred = gen$PRED, METHOD = gen$METHOD)
  est <- EstStep()
  out <- readLines(ctl_path("THEO-FO.OUT.txt"))
  gl <- tail(grep("WITHOUT CONSTANT", out, value = TRUE), 1)
  gold <- as.numeric(regmatches(gl, regexpr("[0-9.]+[ ]*$", gl)))
  expect_equal(est[["Optim"]]$value, gold, tolerance = 1e-3)
})
