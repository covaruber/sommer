test_that("the default uses undamped AI with on-demand recovery", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=30), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(setup$blendWeight, rep(0, 30L))
  expect_equal(setup$stepWeight, rep(1, 30L))
  expect_equal(setup$blend, "em")
  expect_true(isTRUE(attr(setup$covStruct, "directionBlend")))
  expect_equal(attr(setup$covStruct, "blend"), "em")
  expect_true(isTRUE(attr(setup$covStruct, "adaptiveRecovery")))
})

test_that("guarded Aitken acceleration preserves the converged deterministic fit", {
  data(DT_example)
  fit <- function(acceleration){
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=80, likelihoodTol=1e-10, parameterTol=1e-10), update=list(blendWeight=0), proposal=list(extrapolation=if(acceleration == "aitken") "aitken" else "none"), curvature=list(method=if(acceleration == "hybrid-ai") "hybrid-ai" else "exact-ai"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  plain <- fit("none")
  accelerated <- fit("aitken")
  expect_equal(plain$engineDiagnostics$acceleratedProposals, 0L)
  expect_gt(accelerated$engineDiagnostics$acceleratedProposals, 0L)
  expect_true(all(diff(as.numeric(accelerated$llik)) >= -1e-8))
  expect_equal(unlist(accelerated$theta), unlist(plain$theta), tolerance=1e-5)
  expect_equal(as.numeric(accelerated$bu), as.numeric(plain$bu), tolerance=1e-5)
  expect_gte(accelerated$engineDiagnostics$coefficientEvaluations,
    ncol(accelerated$llik) + accelerated$engineDiagnostics$lineSearchHalvings)
})

test_that("acceleration rejects unsupported solver modes", {
  data(DT_example)
  expect_error(mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(update=list(blendWeight=0), proposal=list(extrapolation=if("aitken" %in% c("aitken", "hybrid-ai")) "aitken" else "none"), curvature=list(method=if("aitken" == "hybrid-ai") "hybrid-ai" else "exact-ai"), solve=list(method="pcg"))), "deterministic")
  expect_error(mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", update=list(blendWeight=0), proposal=list(extrapolation=if("aitken" %in% c("aitken", "hybrid-ai")) "aitken" else "none"), curvature=list(method=if("aitken" == "hybrid-ai") "hybrid-ai" else "exact-ai"))), "formulation='direct'")
})
