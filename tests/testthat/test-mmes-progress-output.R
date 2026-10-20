test_that("progress reports iterations without changing numerical results", {
  data(DT_example)
  for(formulation in c("henderson", "direct")){
    fit <- function(verbose){
      mmes(Yield~Env, random=~Name, data=DT_example,
        control=mmes_control(formulation=formulation,
          convergence=list(maxIters=3L),
          diagnostics=list(verbose=verbose, dateWarning=FALSE)))
    }
    output <- capture.output(loud <- fit(TRUE))
    quiet <- fit(FALSE)
    expect_length(grep("^ *iter +LogLik", output), 1L)
    rows <- strsplit(trimws(grep("^ +[0-9]+ +", output, value=TRUE)), " +")
    expect_length(rows, length(loud$llik))
    expect_equal(as.numeric(vapply(rows, `[[`, character(1), 2L)),
                 as.numeric(loud$llik), tolerance=1e-6)
    expect_equal(loud[c("llik", "monitor", "theta", "bu")],
                 quiet[c("llik", "monitor", "theta", "bu")])
  }
})
