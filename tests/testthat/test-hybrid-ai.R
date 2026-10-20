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
  hybridControl <- mmes_control(curvature=list(method="hybrid-ai"))
  expect_equal(hybridControl$curvature$method, "hybrid-ai")
  expect_output(print(hybridControl), "hybrid-ai")
  dataset <- .hybrid_ai_fixture()
  model <- function(overrides=list()){
    control <- utils::modifyList(unclass(hybridControl), overrides,
      keep.null=TRUE)
    control <- utils::modifyList(control, list(update=list(blendWeight=0),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)
    mmes(response~level, random=~id, data=dataset,
      control=do.call(mmes_control, control))
  }
  expect_error(model(list(solve=list(method="pcg"))), "deterministic")
  expect_error(model(list(formulation="direct")), "formulation='direct'")
  expect_error(model(list(task=list(mode="solve"))), "task\\$mode='solve'")
  expect_error(mmes(response~level, random=~id, data=dataset, family=poisson(),
    control=hybridControl), "PQL")
  expect_error(mmes(response~level, random=~id, data=dataset,
    control=hybridControl, .internal=list(pqlInner=TRUE)), "PQL")
  expect_error(model(list(representation=list(factorScores="fixed-shape"))),
    "representation\\$factorScores")
  expect_error(model(list(representation=list(factorScores="profile"))),
    "representation\\$factorScores")
  expect_error(mmes(response~level, random=~vsm(rrm(level, 1L), ism(id)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(hybridControl), list(update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))),
    "latent reduced-rank")
})

test_that("hybrid AI preserves converged fits and reduces actual sensitivity solves", {
  dataset <- .hybrid_ai_fixture()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    model <- function(mode, computeCi=0L, returnParam=FALSE){
      mmes(response~level, random=~vsm(usm(level), ism(id)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(proposal=list(extrapolation=if(mode == "aitken") "aitken" else "none"), curvature=list(method=if(mode == "hybrid-ai") "hybrid-ai" else "exact-ai"))), list(update=list(blendWeight=0), convergence=list(maxIters=100L, likelihoodTol=1e-10, parameterTol=1e-8), task=list(objective=if(reml) "reml" else "ml"), output=list(inverse=computeCi, setup=returnParam), diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))
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
      mmes(response~level, random=~vsm(ar1m(level, rho=0.35, fixed=TRUE), ism(id)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(proposal=list(extrapolation=if(mode == "aitken") "aitken" else "none"), curvature=list(method=if(mode == "hybrid-ai") "hybrid-ai" else "exact-ai"))), list(convergence=list(maxIters=iterations, likelihoodTol=1e-14, parameterTol=1e-14), update=list(blendWeight=0), proposal=list(stepWeight=0.7), output=list(setup=returnParam), diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))
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
      vcc=data.frame(parameter=c(1L, 2L), group=1L),
      control=mmes_control(
        proposal=list(extrapolation=if(mode == "aitken") "aitken" else "none"),
        curvature=list(method=if(mode == "hybrid-ai") "hybrid-ai" else "exact-ai"),
        update=list(blendWeight=0), convergence=list(maxIters=8L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
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
