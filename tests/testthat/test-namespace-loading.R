test_that("the simulation generic is explicitly available in the namespace", {
  expect_identical(get("simulate", envir=asNamespace("sommer")), stats::simulate)
  expect_identical(
    getS3method("simulate", "mmes.glmm.laplace", envir=asNamespace("sommer")),
    get("simulate.mmes.glmm.laplace", envir=asNamespace("sommer")))
})

test_that("the installed namespace loads with only base packages attached", {
  package <- find.package("sommer")
  skip_if_not(file.exists(file.path(package, "Meta", "package.rds")),
    "Requires an installed package rather than a source-tree namespace")
  script <- tempfile(fileext=".R")
  output <- tempfile(fileext=".log")
  on.exit(unlink(c(script, output)), add=TRUE)
  writeLines(c(
    "stopifnot(identical(search(), c('.GlobalEnv', 'Autoloads', 'package:base')))",
    paste0("namespace <- loadNamespace('sommer', lib.loc=",
      encodeString(dirname(package), quote='"'), ")"),
    "stopifnot(identical(get('simulate', namespace), stats::simulate))",
    "stopifnot(is.function(utils::getS3method('simulate', 'mmes.glmm.laplace',",
    "  envir=namespace)))"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "--default-packages=NULL", shQuote(script)),
    stdout=output, stderr=output)
  expect_equal(status, 0L,
    info=paste(readLines(output, warn=FALSE), collapse="\n"))
})
