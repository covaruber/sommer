.hybrid_ai_fixture <- function(){
  set.seed(391)
  dataset <- expand.grid(id=factor(seq_len(24L)), level=factor(seq_len(4L)),
                         replicate=factor(seq_len(3L)))
  covariance <- 0.5^abs(outer(1:4, 1:4, "-"))
  effects <- matrix(rnorm(96L), 24L) %*% chol(covariance)
  dataset$response <- as.integer(dataset$level) / 3 +
    effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
    rnorm(nrow(dataset), sd=0.7)
  dataset
}

.hybrid_ai_at <- function(setup, fit, em=0, computeCi=0L){
  covariance <- setup$covStruct
  for(term in seq_along(covariance)){
    covariance[[term]]$par <- as.numeric(fit$covParWorking[[term]])
  }
  attr(covariance, "acceleration") <- "none"
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    em, 0.9, FALSE, computeCi=computeCi, solver=setup$solver, reml=setup$REML,
    responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
    preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
}

test_that("hybrid AI controls merge, print, and reject unsupported combinations", {
  control <- mmes_control(acceleration="hybrid-ai")
  expect_equal(control$acceleration, "hybrid-ai")
  expect_output(print(control), "hybrid-ai")
  expect_equal(sommer:::.merge_mmes_control(engine=mmes_control(),
    acceleration="hybrid-ai")$acceleration, "hybrid-ai")
  expect_equal(sommer:::.merge_mmes_control(acceleration="hybrid-ai")$acceleration,
               "hybrid-ai")
  dataset <- .hybrid_ai_fixture()
  model <- function(...){
    mmes(response~level, random=~id, data=dataset, blendWeight=0, engine=control,
      verbose=FALSE, dateWarning=FALSE, ...)
  }
  expect_error(model(solver="pcg"), "deterministic")
  expect_error(model(henderson=FALSE), "deterministic")
  expect_error(model(solveOnly=TRUE), "without solveOnly")
  expect_error(model(family=poisson()), "without solveOnly")
  expect_error(model(.internal=list(pqlInner=TRUE)), "without solveOnly")
  expect_error(model(factorScoreAugmentation="fixed-shape"), "factor-score augmentation")
  expect_error(model(factorScoreAugmentation="profile"), "factor-score augmentation")
  expect_error(mmes(response~level, random=~vsm(rrm(level, 1L), ism(id)),
    data=dataset, blendWeight=0, engine=control, verbose=FALSE, dateWarning=FALSE),
    "latent reduced-rank")
})

test_that("hybrid AI preserves converged fits and reduces actual sensitivity solves", {
  dataset <- .hybrid_ai_fixture()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    model <- function(mode, computeCi=0L, returnParam=FALSE){
      mmes(response~level, random=~vsm(usm(level), ism(id)), data=dataset,
        blendWeight=0,
        nIters=100L, tolParConvLL=1e-10, tolParConvNorm=1e-8,
        REML=reml, computeCi=computeCi, returnParam=returnParam,
        engine=mmes_control(solver=solver, acceleration=mode,
          experimental=list(profile="force")),
        verbose=FALSE, dateWarning=FALSE)
    }
    plain <- model("none")
    hybrid <- model("hybrid-ai")
    diagnostics <- hybrid$engineDiagnostics$hybridAI
    expect_true(plain$convergence)
    expect_true(hybrid$convergence)
    expect_true(diagnostics$active)
    expect_gt(diagnostics$bfgsUpdates, 0)
    expect_gt(diagnostics$certifications, 0)
    expect_true(diagnostics$finalInformationExact)
    expect_lte(diagnostics$finalProjectedScoreNorm, 1e-4)
    expect_lt(diagnostics$sensitivityRHSColumns,
              plain$engineDiagnostics$hybridAI$sensitivityRHSColumns)
    expect_gt(diagnostics$sensitivityRHSColumnsSkipped, 0)
    expect_equal(hybrid$engineDiagnostics$profile$counters$sensitivityRHSColumns,
                 diagnostics$sensitivityRHSColumns)
    expect_true(all(diff(c(hybrid$llik)) >= -1e-8))
    expect_equal(tail(c(hybrid$llik), 1L), tail(c(plain$llik), 1L), tolerance=1e-8)
    for(field in c("theta", "covPar", "bu", "theta_se")){
      expect_equal(hybrid[[field]], plain[[field]], tolerance=2e-5)
    }
    reference <- .hybrid_ai_at(model("none", returnParam=TRUE), hybrid)
    for(field in c("InfMat", "theta_se", "dLu")){
      expect_equal(hybrid[[field]], reference[[field]], tolerance=1e-8)
    }
    if(reml){
      for(mode in 1:2){
        uncertainty <- model("hybrid-ai", mode)
        expect_equal(uncertainty$theta_se, hybrid$theta_se, tolerance=1e-8)
        expect_true(uncertainty$engineDiagnostics$hybridAI$finalInformationExact)
        expect_true(all(vapply(uncertainty$uPevList, length, integer(1)) > 0L))
      }
    }
  }
})

test_that("hybrid AI returns exact information at iteration limits and skips fixed directions", {
  dataset <- .hybrid_ai_fixture()
  for(iterations in c(1L, 2L, 3L, 6L)){
    model <- function(mode, returnParam=FALSE){
      mmes(response~level,
        random=~vsm(ar1m(level, rho=0.35, fixed=TRUE), ism(id)), data=dataset,
        nIters=iterations, tolParConvLL=1e-14, tolParConvNorm=1e-14,
        blendWeight=0, stepWeight=0.7, returnParam=returnParam,
        engine=mmes_control(solver="ldlt", acceleration=mode,
          experimental=list(profile="force")),
        verbose=FALSE, dateWarning=FALSE)
    }
    hybrid <- model("hybrid-ai")
    reference <- .hybrid_ai_at(model("none", TRUE), hybrid, em=0)
    expect_true(hybrid$engineDiagnostics$hybridAI$finalInformationExact)
    expect_equal(hybrid$engineDiagnostics$profile$counters$freeCovarianceParameters, 2L)
    expect_false(hybrid$convergence)
    expect_equal(as.numeric(hybrid$covPar[[1]])[2L], 0.35, tolerance=1e-14)
    for(field in c("InfMat", "theta_se", "dLu")){
      expect_equal(hybrid[[field]], reference[[field]], tolerance=1e-9)
    }
    if(iterations == 1L){
      plain <- model("none")
      for(field in c("llik", "monitor", "theta", "bu", "InfMat", "theta_se", "dLu")){
        expect_identical(hybrid[[field]], plain[[field]])
      }
    }else{
      expect_gt(hybrid$engineDiagnostics$hybridAI$sensitivityRHSColumnsSkipped, 0)
      if(iterations == 2L){
        plain <- model("none")
        for(field in c("monitor", "theta", "bu", "InfMat", "theta_se", "dLu")){
          expect_equal(hybrid[[field]], plain[[field]], tolerance=1e-12)
        }
      }
    }
  }
})

test_that("tied covariance parameters explicitly retain the exact optimizer", {
  dataset <- .hybrid_ai_fixture()
  model <- function(mode){
    mmes(response~level, random=~id + replicate, data=dataset,
      blendWeight=0, vcc=data.frame(parameter=c(1L, 2L), group=1L), nIters=8L,
      engine=mmes_control(acceleration=mode), verbose=FALSE, dateWarning=FALSE)
  }
  plain <- model("none")
  hybrid <- model("hybrid-ai")
  expect_false(hybrid$engineDiagnostics$hybridAI$active)
  expect_equal(hybrid$engineDiagnostics$hybridAI$fallbackReason, "constraint-map")
  expect_equal(hybrid$engineDiagnostics$hybridAI$bfgsUpdates, 0L)
  for(field in c("llik", "monitor", "theta", "bu", "InfMat", "theta_se", "dLu")){
    expect_identical(hybrid[[field]], plain[[field]])
  }
})

test_that("hybrid safeguards retain exact output after large proposals and fixed models", {
  dataset <- .hybrid_ai_fixture()
  model <- function(returnParam=FALSE){
    mmes(response~level, random=~vsm(usm(level), ism(id)), data=dataset,
      nIters=12L, blendWeight=0, stepWeight=30,
      tolParConvLL=1e-12, tolParConvNorm=1e-12, returnParam=returnParam,
      engine=mmes_control(acceleration="hybrid-ai", solver="ldlt"),
      verbose=FALSE, dateWarning=FALSE)
  }
  fit <- model()
  expect_gt(fit$engineDiagnostics$lineSearchHalvings, 0L)
  expect_true(all(diff(c(fit$llik)) >= -1e-8))
  expect_true(fit$engineDiagnostics$hybridAI$finalInformationExact)
  reference <- .hybrid_ai_at(model(TRUE), fit, em=0)
  for(field in c("InfMat", "theta_se", "dLu")){
    expect_equal(fit[[field]], reference[[field]], tolerance=1e-8)
  }

  fixed <- mmes(response~level,
    random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=0.5, fixedSigma2=TRUE), data=dataset,
    nIters=4L, blendWeight=0, engine=mmes_control(acceleration="hybrid-ai"),
    verbose=FALSE, dateWarning=FALSE)
  expect_true(fixed$engineDiagnostics$hybridAI$finalInformationExact)
  expect_equal(as.numeric(fixed$theta_se), rep(0, 4L))
  expect_equal(unlist(fixed$covPar, use.names=FALSE), c(1, 0.5), tolerance=1e-12)
})

test_that("hybrid final information preserves covariance SE and small-sample Wald consumers", {
  skip_on_cran()
  data(DT_yatesoats, package="enhancer", envir=environment())
  model <- function(mode){
    mmes(Y~V*N, random=~B+B:MP, data=DT_yatesoats,
      nIters=100L, tolParConvLL=1e-12, tolParConvNorm=1e-10,
      blendWeight=0,
      engine=mmes_control(acceleration=mode), verbose=FALSE, dateWarning=FALSE)
  }
  plain <- model("none")
  hybrid <- model("hybrid-ai")
  expect_true(hybrid$convergence)
  expect_equal(covparams_mmes_se(hybrid), covparams_mmes_se(plain), tolerance=1e-5)
  for(method in c("satterthwaite", "kr")){
    expect_equal(wald_mmes(hybrid, denDF=method), wald_mmes(plain, denDF=method),
                 tolerance=1e-4)
  }
})
