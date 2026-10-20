.information_data <- function(){
  dataset <- expand.grid(id=factor(seq_len(10L)), replicate=seq_len(4L))
  dataset$x <- rep(c(-1, 0, 1, 2), each=10L)
  dataset$y <- cos(seq_len(nrow(dataset)) * 0.8) +
    as.integer(dataset$id) / 5 + dataset$x / 3
  dataset
}

.information_fit <- function(mode="full", solver="ldlt", blend=0,
                             acceleration="none", inverse=0L, REML=TRUE){
  arguments <- list(fixed=y~x, random=~id, data=.information_data(),
    control=mmes_control(
      solve=list(method=solver),
      task=list(objective=if(REML) "reml" else "ml"),
      convergence=list(maxIters=4L),
      output=list(inverse=inverse, covariance=mode),
      update=list(blendWeight=blend),
      proposal=list(extrapolation=if(acceleration == "aitken") "aitken" else "none"),
      curvature=list(method=if(acceleration == "hybrid-ai") "hybrid-ai" else "exact-ai"),
      score=if(solver == "pcg") list(traceProbes=12L) else list(),
      logdet=if(solver == "pcg") list(lanczosSteps=40L) else list(),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  suppressMessages(do.call(mmes, arguments))
}

test_that("covariance information controls validate and preserve defaults", {
  expect_identical(mmes_control()$output$covariance, "full")
  for(value in list(NULL, NA_character_, "", c("full", "none"), TRUE, 1)){
    expect_error(mmes_control(output=list(covariance=value)))
  }
  expect_error(mmes_control(output=list(unknownInformation="none")),
    "Unknown output settings")
})

test_that("information omission preserves fits and coefficient uncertainty across solvers", {
  for(solver in c("auto", "ldlt", "cholmod", "pcg")){
    for(REML in if(solver == "pcg") TRUE else c(TRUE, FALSE)){
      reference <- .information_fit(solver=solver, inverse=2L, REML=REML)
      candidate <- .information_fit("none", solver=solver, inverse=2L, REML=REML)
      for(field in c("bu", "C", "Ci", "llik", "theta", "dLu", "monitor",
                    "convergence", "uPevList")){
        expect_identical(candidate[[field]], reference[[field]], info=paste(solver, REML, field))
      }
      expect_null(candidate$InfMat)
      expect_null(candidate$theta_se)
      expect_identical(candidate$control$output$covariance, "none")
      expect_identical(candidate$control$output$covariance, "none")
      expect_equal(covparams_mmes(candidate), covparams_mmes(reference))
      expect_true(all(is.na(candidate$covParNative$StdError)))
      expect_true(all(is.na(candidate$covParNative$Zratio)))
      expect_equal(fitted(candidate), fitted(reference))
      expect_equal(residuals(candidate), residuals(reference))
      expect_output(print(summary(candidate)), "Covariance information omitted")
      expect_equal(covmatrix_mmes(candidate, 1L, se=FALSE),
        covmatrix_mmes(reference, 1L, se=FALSE))
    }
  }
})

test_that("methods requiring omitted covariance information fail explicitly", {
  candidate <- .information_fit("none")
  error <- "refit with output=list\\(covariance='full'\\)"
  expect_error(covparams_mmes_se(candidate), error)
  expect_error(covmatrix_mmes(candidate, 1L), error)
  expect_error(vpredict(candidate, ~V1/(V1+V2)), error)
  for(method in c("satterthwaite", "kr")){
    expect_error(wald_mmes(candidate, denDF=method), error)
  }
  expect_s3_class(wald_mmes(candidate, denDF="none"), "wald.mmes")
  expect_equal(predict(candidate, D="id")$pvals,
    predict(.information_fit(), D="id")$pvals)
  expect_error(predict(candidate, D="id", pairwise=TRUE, df="kr"), error)
  expect_equal(postPEV(candidate, mode=2L)$Ci,
    postPEV(.information_fit(), mode=2L)$Ci)
})

test_that("information omission rejects unsupported engines before fitting", {
  baseControl <- mmes_control(output=list(covariance="none"),
    diagnostics=list(verbose=FALSE))
  for(settings in list(list(formulation="direct"),
                       list(task=list(mode="solve")), NULL)){
    arguments <- list(fixed=y~x, random=~id, data=.information_data())
    expect_error(suppressMessages({
      if(!is.null(settings)){
        arguments$control <- do.call(mmes_control,
          utils::modifyList(unclass(baseControl), settings, keep.null=TRUE))
      }else{
        arguments$control <- baseControl
        arguments$family <- poisson()
      }
      do.call(mmes, arguments)
    }))
  }
})
