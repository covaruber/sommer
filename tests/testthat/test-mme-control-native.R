test_that("native control defaults and parsers preserve their contracts", {
  .native_run("test-mme-control", "mme_control.h", strict=TRUE)
})

test_that("native control transport respects descriptor and option precedence", {
  source <- test_path("..", "native", "test-mme-control-transport.cpp")
  header <- test_path("..", "..", "src", "mme_control_r.h")
  .native_require(all(file.exists(c(source, header))),
    "Native checks require the package source tree")
  .native_require(nzchar(Sys.which("c++")), "C++ compiler unavailable")
  root <- normalizePath(test_path("..", ".."))
  previous <- options("sommer.mme.testMode", "sommer.mme.memoryPolicy")
  directory <- getwd()
  (function(){
    previous <- setwd(root)
    on.exit(setwd(previous), add=TRUE)
    sys.source("tests/native/test-mme-control-transport.R",
      envir=new.env(parent=globalenv()))
  })()
  expect_identical(lapply(names(previous), getOption), unname(previous))
  expect_identical(getwd(), directory)
  succeed()
})
