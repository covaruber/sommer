.cfree_fit <- function(retain=TRUE, information="full", REML=TRUE){
  dataset <- expand.grid(id=factor(seq_len(8L)), group=factor(seq_len(3L)),
    replicate=seq_len(3L))
  dataset$y <- sin(seq_len(nrow(dataset)) * 0.7) +
    as.integer(dataset$id) / 4 + as.integer(dataset$group) / 3
  suppressMessages(mmes(y~group, random=~id, data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(representation=list(equations="saddle-point"), output=list(coefficientMatrix=retain, covariance=information))), list(convergence=list(maxIters=3L), update=list(blendWeight=0), task=list(objective=if(REML) "reml" else "ml"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))))
}

test_that("coefficient-matrix retention defaults and controls validate", {
  expect_true(mmes_control()$output$coefficientMatrix)
  expect_error(mmes_control(output=list(coefficientMatrix=1)),
    "output\\$coefficientMatrix")
  dataset <- data.frame(y=seq_len(6L), id=factor(rep(seq_len(3L), 2L)))
  expect_error(mmes(y~1, random=~id, data=dataset, control=mmes_control(
    output=list(coefficientMatrix=FALSE, setup=TRUE),
    diagnostics=list(verbose=FALSE))),
    "saddle-point")
  expect_error(mmes(y~1, random=~id, data=dataset, control=mmes_control(
    representation=list(equations="saddle-point"),
    output=list(coefficientMatrix=FALSE, setup=TRUE, inverse=2L),
    diagnostics=list(verbose=FALSE))), "output\\$inverse=0")
})

test_that("C-free saddle fits skip compatibility assembly and preserve statistical outputs", {
  for(REML in c(TRUE, FALSE)){
    for(information in c("full", "none")){
      reference <- .cfree_fit(information=information, REML=REML)
      candidate <- .cfree_fit(FALSE, information=information, REML=REML)
      expect_equal(reference$engineDiagnostics$saddlePoint$compatibilityAssemblies, 1L)
      expect_equal(candidate$engineDiagnostics$saddlePoint$compatibilityAssemblies, 0L)
      expect_false(candidate$engineDiagnostics$saddlePoint$compatibilityCMaterialized)
      expect_equal(candidate$engineDiagnostics$saddlePoint$normalEquationAssemblies, 0L)
      expect_null(candidate$C)
      expect_false(candidate$control$output$coefficientMatrix)
      for(field in c("bu", "b", "u", "llik", "covPar", "theta", "InfMat",
                    "theta_se", "monitor", "convergence")){
        expect_identical(candidate[[field]], reference[[field]], info=field)
      }
      expect_equal(fitted(candidate), fitted(reference))
      expect_equal(residuals(candidate), residuals(reference))
      expect_equal(coef(candidate), coef(reference))
      expect_equal(covparams_mmes(candidate), covparams_mmes(reference))
      expect_output(print(summary(candidate)), "Coefficient matrix C was not retained")
      if(information == "full"){
        expect_equal(covparams_mmes_se(candidate), covparams_mmes_se(reference))
        expect_equal(vpredict(candidate, ~V1/(V1+V2)),
          vpredict(reference, ~V1/(V1+V2)))
      }
    }
  }
})

test_that("C-free objects fail explicitly without reconstructing C", {
  candidate <- unserialize(serialize(.cfree_fit(FALSE), NULL))
  error <- "refit with output=list\\(coefficientMatrix=TRUE\\)"
  expect_error(postPEV(candidate), error)
  expect_error(postVarU(candidate), error)
  expect_error(predict(candidate, D="id"), error)
  expect_error(wald_mmes(candidate), error)
  expect_error(sommer:::post_mme_Cinverse_cpp(candidate, 2L),
    "refit with retainCoefficientMatrix=TRUE")
  contrast <- Matrix::Matrix(matrix(1, 1L, nrow(candidate$bu)), sparse=TRUE)
  expect_error(sommer:::predict_mmes_vcov_cpp(candidate, contrast),
    "refit with retainCoefficientMatrix=TRUE")
  expect_null(candidate$C)
  expect_equal(candidate$engineDiagnostics$saddlePoint$compatibilityAssemblies, 0L)
})
