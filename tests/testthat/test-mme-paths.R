test_that("MME path policies preserve dispatch and resource boundaries", {
  root <- test_path("..", "..")
  source <- file.path(root, "tests", "native", "test-mme-paths.cpp")
  header <- file.path(root, "src", "mme_paths.h")
  skip_if(!file.exists(header), "Native path checks require the package source tree")
  compiler <- Sys.which("c++")
  skip_if(!nzchar(compiler), "C++ compiler unavailable")
  binary <- tempfile("sommer-mme-paths-")
  on.exit(unlink(binary), add=TRUE)
  output <- system2(compiler,
    c("-std=c++11", "-Wall", "-Wextra", "-pedantic", shQuote(source), "-o", shQuote(binary)),
    stdout=TRUE, stderr=TRUE)
  expect_null(attr(output, "status"), info=paste(output, collapse="\n"))
  if(!file.exists(binary)) return(invisible(NULL))
  output <- system2(binary, stdout=TRUE, stderr=TRUE)
  expect_null(attr(output, "status"), info=paste(output, collapse="\n"))
})