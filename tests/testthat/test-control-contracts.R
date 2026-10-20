test_that("every private PCG setting has exactly one public owner", {
  mapping <- sommer:::.mme_pcg_mapping()
  control <- mmes_control()
  expect_setequal(names(mapping), names(sommer:::.mme_pcg_control()))
  expect_identical(anyDuplicated(mapping), 0L)
  expect_true(all(mapping %in% sommer:::.mme_control_paths(control)))
  for(name in names(mapping)){
    expect_identical(sommer:::.mme_pcg_public(name), mapping[[name]])
  }
  engine <- sommer:::.mmes_compile_control(control)
  for(name in setdiff(names(mapping), "compression")){
    expect_equal(engine$pcg[[name]],
      sommer:::.mme_path_value(control, mapping[[name]]), info=name)
  }
  expect_error(sommer:::.mme_pcg_public("unknown"), "Unknown internal")
})

test_that("mutable controls and recursive transport do not use stale compilation", {
  control <- mmes_control(solve=list(method="pcg"),
    execution=list(seed=9L), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  control$execution$seed <- 27L
  control$score$traceProbes <- 3L
  control$convergence$maxTraceProbes <- 3L
  control$output$setup <- TRUE
  data <- data.frame(y=c(1, 2, 1.5, 3, 2.5, 4), id=factor(rep(1:3, 2)))
  setup <- mmes(y~1, random=~id, data=data, control=control)
  expect_identical(setup$engine$pcg$seed, 27L)
  expect_identical(setup$engine$pcg$traceProbes, 3L)
})
