.reml_estimability_data <- function(){
  set.seed(418)
  d <- expand.grid(environment=factor(c("a", "b", "c")),
                   id=factor(seq_len(8)), replicate=seq_len(3))
  d$block <- factor(ifelse(d$environment == "a", "1",
                           ifelse(as.integer(d$id) <= 4L, "1", "2")))
  d$y <- rnorm(nrow(d)) + rep(c(-1, 0.5, 1), length.out=nrow(d))
  d
}

.reml_estimability_control <- function(...){
  mmes_control(convergence=list(maxIters=4),
               diagnostics=list(verbose=FALSE, dateWarning=FALSE), ...)
}

test_that("REML omission works without runtime BLAS control", {
  local_mocked_bindings(.mmes_blas_controller=function() NULL, .package="sommer")
  d <- .reml_estimability_data()
  control <- .reml_estimability_control()
  expect_warning(fit <- mmes(y~environment,
    random=~vsm(dsm(environment), ism(block)), data=d, control=control),
    "Omitting 1 structurally unestimable")
  expect_identical(fit$threadingDiagnostics$controller, "unavailable")
  expect_null(fit$threadingDiagnostics$blasThreadsRequested)
  expect_true(is.na(fit$threadingDiagnostics$blasThreadsConfigured))
  expect_identical(covparams_mmes(fit, 1)$estimate[1], 0)
})

test_that("REML omits only independent block variances in the fixed space", {
  d <- .reml_estimability_data()
  for(solver in c("ldlt", "cholmod")){
    control <- .reml_estimability_control(solve=list(method=solver), output=list(inverse=2))
    expect_warning(fit <- mmes(y~environment,
      random=~vsm(dsm(environment), ism(block)), data=d, control=control),
      "Omitting 1 structurally unestimable")
    reference <- mmes(y~environment,
      random=~vsm(atm(environment, c("b", "c")), ism(block)), data=d, control=control)
    report <- covparams_mmes_se(fit, 1)
    expect_equal(report$parameter, paste0("variance[", c("a", "b", "c"), "]"))
    expect_identical(report$estimate[1], 0)
    expect_true(is.na(report$StdError[1]))
    expect_true(is.na(report$Zratio[1]))
    expect_identical(report$status[1], "unestimable (omitted)")
    expect_equal(report$estimate[-1], covparams_mmes(reference, 1)$estimate, tolerance=1e-7)
    expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(reference)), tolerance=1e-8)
    expect_equal(as.numeric(fit$b), as.numeric(reference$b), tolerance=1e-8)
    expect_equal(as.numeric(fit$llik), as.numeric(reference$llik), tolerance=1e-8)
    expect_equal(fit$AIC, reference$AIC, tolerance=1e-8)
    expect_equal(fit$BIC, reference$BIC, tolerance=1e-8)
    expect_equal(as.numeric(fit$uList[[1]][,1]), c(0, 0))
    expect_equal(as.numeric(fit$uPevList[[1]][,1]), c(0, 0))
    expect_equal(nrow(fit$coefficientMap), nrow(fit$bu))
    expect_equal(ncol(fit$coefficientMap), nrow(fit$C))
    expect_equal(fit$engineDiagnostics$remlOmissions, fit$remlOmissions)
    expect_identical(summary(fit)$varcomp$status[1], "unestimable (omitted)")
    expect_output(print(summary(fit)), "zeros are not variance estimates")
    covariance <- covmatrix_mmes(fit, 1)
    expect_equal(covariance$covariance, fit$theta[[1]], ignore_attr=TRUE)
    expect_true(all(is.na(covariance$covariance.se[1,])))
    expect_true(all(is.na(covariance$correlation[1,])))
    D <- cbind(diag(length(fit$b)), matrix(0, length(fit$b), length(fit$u)))
    Dr <- cbind(diag(length(reference$b)),
                matrix(0, length(reference$b), length(reference$u)))
    expect_equal(sommer:::.mmes_vcov_solver(fit)(D),
                 sommer:::.mmes_vcov_solver(reference)(Dr), tolerance=1e-8)
    post <- postPEV(fit)
    expect_equal(as.numeric(post$uPevList[[1]][,1]), c(0, 0))
  }
})

test_that("whole null terms are removed while public effects remain zero", {
  d <- .reml_estimability_data()
  control <- .reml_estimability_control(output=list(inverse=2))
  expect_warning(fit <- mmes(y~environment,
    random=~vsm(ism(environment)) + vsm(ism(id)), data=d, control=control),
    "Omitting 1 structurally unestimable")
  reference <- mmes(y~environment, random=~vsm(ism(id)), data=d, control=control)
  expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(reference)), tolerance=1e-8)
  expect_equal(fit$AIC, reference$AIC, tolerance=1e-8)
  expect_true(all(fit$uList[[1]] == 0))
  expect_true(all(fit$theta[[1]] == 0))
  expect_equal(covparams_mmes(fit, 1)$estimate, 0)
  expect_true(is.na(covparams_mmes_se(fit, 1)$StdError))
  expect_equal(nrow(fit$theta_se), sum(lengths(fit$covPar)))

  expect_warning(allNull <- mmes(y~environment,
    random=~vsm(dsm(environment), ism(environment)), data=d, control=control),
    "Omitting 3 structurally unestimable")
  noRandom <- mmes(y~environment, data=d, control=control)
  expect_equal(as.numeric(fitted(allNull)), as.numeric(fitted(noRandom)), tolerance=1e-8)
  expect_true(all(covparams_mmes(allNull, 1)$estimate == 0))
  expect_true(all(is.na(covparams_mmes_se(allNull, 1)$StdError)))
  expect_equal(allNull$AIC, noRandom$AIC, tolerance=1e-8)
})

test_that("shared, fixed, correlated and ML components are not automatically omitted", {
  d <- .reml_estimability_data()
  setup <- function(random, objective="reml", fixed=y~environment){
    mmes(fixed, random=random, data=d, control=.reml_estimability_control(
      task=list(objective=objective), output=list(setup=TRUE)))
  }
  expect_null(setup(~vsm(ism(block)))$remlOmissionInfo)
  expect_null(setup(~vsm(dsm(environment), ism(block)), "ml")$remlOmissionInfo)
  expect_null(setup(~vsm(usm(environment), ism(block)))$remlOmissionInfo)
  expect_null(setup(~vsm(dsm(environment, fixed=c(TRUE, FALSE)), ism(block)))$remlOmissionInfo)
  expect_null(setup(~vsm(dsm(environment), ism(block), fixedSigma2=TRUE))$remlOmissionInfo)
  expect_null(setup(~vsm(dsm(environment), ism(block)), fixed=y~1)$remlOmissionInfo)
  near <- matrix(as.numeric(d$environment == "a") + 1e-6 * rnorm(nrow(d)), ncol=1)
  expect_null(setup(~vsm(ism(near)))$remlOmissionInfo)
})

test_that("omission uses retained observations and supports a single retained level", {
  d <- .reml_estimability_data()
  d$y[d$environment == "c" & d$block == "2"] <- NA_real_
  control <- .reml_estimability_control(output=list(inverse=1))
  expect_warning(fit <- mmes(y~environment,
    random=~vsm(dsm(environment), ism(block)), data=d, control=control),
    "Omitting 2 structurally unestimable")
  expect_equal(fit$remlOmissions$parameter, c("variance[a]", "variance[c]"))
  expect_true(all(is.na(covparams_mmes_se(fit, 1)$StdError[c(1,3)])))
  expect_equal(dim(fit$theta[[1]]), c(3L,3L))
  expect_equal(ncol(fit$uList[[1]]), 3L)
  expect_equal(nrow(fit$theta_se), sum(lengths(fit$covPar)))
})

test_that("omission composes with latent reduced-rank coefficient maps", {
  d <- .reml_estimability_data()
  control <- .reml_estimability_control(output=list(inverse=2))
  expect_warning(fit <- mmes(y~environment,
    random=~vsm(rrm(environment, 1), ism(id)) +
      vsm(dsm(environment), ism(block)), data=d, control=control),
    "Omitting 1 structurally unestimable")
  reference <- mmes(y~environment,
    random=~vsm(rrm(environment, 1), ism(id)) +
      vsm(atm(environment, c("b", "c")), ism(block)), data=d, control=control)
  expect_true(fit$engineDiagnostics$latentReducedRank)
  expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(reference)), tolerance=1e-7)
  expect_equal(fit$AIC, reference$AIC, tolerance=1e-7)
  expect_equal(ncol(fit$coefficientMap), nrow(fit$C))
  D <- cbind(diag(length(fit$b)), matrix(0, length(fit$b), length(fit$u)))
  Dr <- cbind(diag(length(reference$b)),
              matrix(0, length(reference$b), length(reference$u)))
  expect_equal(sommer:::.mmes_vcov_solver(fit)(D),
               sommer:::.mmes_vcov_solver(reference)(Dr), tolerance=1e-7)
})

test_that("direct and PCG fits use the same structural reduction", {
  d <- .reml_estimability_data()
  for(formulation in c("direct", "henderson")){
    control <- .reml_estimability_control(formulation=formulation,
      solve=if(formulation == "direct") list() else list(method="pcg"),
      output=list(inverse=2))
    expect_warning(fit <- mmes(y~environment,
      random=~vsm(dsm(environment), ism(block)), data=d, control=control),
      "Omitting 1 structurally unestimable")
    reference <- mmes(y~environment,
      random=~vsm(atm(environment, c("b", "c")), ism(block)), data=d, control=control)
    expect_equal(as.numeric(fit$b), as.numeric(reference$b), tolerance=1e-8)
    expect_equal(unname(fit$uList[[1]][,-1,drop=FALSE]),
                 unname(reference$uList[[1]]), tolerance=1e-8)
    expect_equal(as.numeric(fitted(fit)),
                 as.numeric(fit$W %*% fit$coefficientMap %*% reference$bu), tolerance=1e-8)
    expect_equal(fit$AIC, reference$AIC, tolerance=1e-8)
    expect_true(is.na(covparams_mmes_se(fit, 1)$StdError[1]))
  }
})

test_that("small-sample inference rebuilds only retained covariance terms", {
  d <- .reml_estimability_data()
  control <- .reml_estimability_control(output=list(inverse=2))
  expect_warning(fit <- mmes(y~environment,
    random=~vsm(ism(environment)) + vsm(ism(id)), data=d, control=control),
    "Omitting 1 structurally unestimable")
  reference <- mmes(y~environment, random=~vsm(ism(id)), data=d, control=control)
  inputs <- sommer:::.mmes_rebuild_inputs(fit)
  information <- sommer:::.mmes_fixed_information_at(fit, inputs, fit$covPar)
  referenceInformation <- sommer:::.mmes_fixed_information_at(reference,
    sommer:::.mmes_rebuild_inputs(reference), reference$covPar)
  expect_equal(information, referenceInformation, tolerance=1e-8)
  expect_equal(wald_mmes(fit, denDF="satterthwaite")$denDF,
               wald_mmes(reference, denDF="satterthwaite")$denDF, tolerance=1e-6)
})

test_that("omission respects covariance-free and coefficient-free outputs", {
  d <- .reml_estimability_data()
  for(output in list(list(covariance="none"),
                     list(covariance="none", coefficientMatrix=FALSE))){
    control <- .reml_estimability_control(output=output,
      representation=if(isFALSE(output$coefficientMatrix))
        list(equations="saddle-point") else list())
    expect_warning(fit <- mmes(y~environment,
      random=~vsm(dsm(environment), ism(block)), data=d, control=control),
      "Omitting 1 structurally unestimable")
    expect_true(all(is.na(summary(fit)$varcomp$StdError)))
    expect_identical(covparams_mmes(fit, 1)$estimate[1], 0)
    expect_identical(summary(fit)$varcomp$status[1], "unestimable (omitted)")
  }
})
