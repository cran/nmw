copy_example_run <- function() {
  source_dir <- system.file("extdata", package = "nmw")
  stopifnot(nzchar(source_dir))

  run_dir <- tempfile("nmw-report-test-")
  dir.create(run_dir)
  copied <- file.copy(list.files(source_dir, full.names = TRUE), run_dir)
  stopifnot(all(copied))
  run_dir
}


drop_sdtab_columns <- function(run_dir, columns) {
  path <- file.path(run_dir, "sdtab")
  table_header <- readLines(path, n = 1L)
  values <- read.table(path, header = TRUE, skip = 1L,
                       na.strings = "***********")
  values[intersect(columns, colnames(values))] <- NULL

  writeLines(table_header, path)
  suppressWarnings(write.table(values, path, append = TRUE, quote = FALSE,
                               row.names = FALSE, col.names = TRUE))
}


test_that("output report honors a custom filename and returns metadata", {
  run_dir <- copy_example_run()
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  result <- NULL
  expect_message(
    result <- nmw_report_output(run_dir, file = "custom-output.pdf",
                                model = "ignored"),
    "generated")

  expect_true(file.exists(file.path(run_dir, "custom-output.pdf")))
  expect_identical(result$file, "custom-output.pdf")
  expect_gt(result$pages, 0L)
  expect_gt(result$lines, 0L)
})


test_that("residual report omits missing residual diagnostics", {
  run_dir <- copy_example_run()
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  drop_sdtab_columns(run_dir, "IWRE")
  result <- NULL
  expect_message(
    result <- nmw_report_resid(run_dir, file = "missing-iwre.pdf",
                               model = "2006"),
    "generated")

  expect_true(file.exists(file.path(run_dir, "missing-iwre.pdf")))
  expect_setequal(result$residuals, c("WRES", "CWRE"))
  expect_identical(result$omittedResiduals, "IWRE")
  expect_gt(result$pages, 0L)

  drop_sdtab_columns(run_dir, c("WRES", "CWRES"))
  no_residuals <- NULL
  expect_message(
    no_residuals <- nmw_report_resid(run_dir, file = "no-residuals.pdf",
                                     model = "2006"),
    "no usable residuals")

  expect_true(file.exists(file.path(run_dir, "no-residuals.pdf")))
  expect_length(no_residuals$residuals, 0L)
  expect_setequal(no_residuals$omittedResiduals, c("WRES", "CWRE", "IWRE"))
})
