test_that("the README quick start runs and returns usable uncertainty", {
  source <- test_path("..", "..", "README.md")
  skip_if(!file.exists(source), "README requires the package source tree")
  lines <- readLines(source, warn=FALSE)
  start <- which(lines == "```r")
  expect_length(start, 1L)
  end <- which(seq_along(lines) > start & lines == "```")[1L]
  example <- new.env(parent=globalenv())
  invisible(capture.output(eval(parse(text=lines[seq.int(start + 1L, end - 1L)]),
    envir=example)))
  expect_s3_class(example$fit, "mmes")
  expect_true(isTRUE(example$fit$convergence))
  expect_true(all(is.finite(example$genotypes$pvals$std.error)))
  expect_gt(nrow(example$genotypes$pvals), 0L)
})
