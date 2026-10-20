test_that("blend weights validate and recycle to the iteration count", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=5L), update=list(blendWeight=c(0, 0.4)), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))

  expect_equal(setup$blendWeight, c(0, 0.4, 0, 0.4, 0))
  expect_true(isTRUE(attr(setup$covStruct, "directionBlend")))

  expect_error(mmes(Yield~Env, data=DT_example, control=mmes_control(update=list(blendWeight=1.1))),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, control=mmes_control(update=list(blendWeight=Inf))),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, control=mmes_control(update=list(blendWeight="0.5"))),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, control=mmes_control(convergence=list(maxIters=1L), update=list(blendWeight=c(1, NA)))),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, control=mmes_control(convergence=list(maxIters=1e20))), "positive integer")
})

test_that("scalar blend counts generate capped exponential decay schedules", {
  data(DT_example)
  setup <- function(weight, iterations=6L){
    mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(convergence=list(maxIters=iterations), update=list(blendWeight=weight), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  expect_equal(setup(3)$blendWeight,
    c(exp(seq(log(1), log(0.03), length.out=3L)), rep(0.03, 3L)))
  expect_equal(setup(2)$blendWeight, c(1, rep(0.03, 5L)))
  expect_equal(setup(10, 4L)$blendWeight,
    exp(seq(log(1), log(0.03), length.out=4L)))
  expect_equal(setup(10, 1L)$blendWeight, 1)
  expect_equal(setup(1)$blendWeight, rep(1, 6L))
  expect_equal(setup(0.4)$blendWeight, rep(0.4, 6L))
  expect_equal(setup(0)$blendWeight, rep(0, 6L))
  expect_error(setup(2.5), "whole-number")
  for(weight in list(c(2, 3), c(0.5, 2), c(2, 0.5))){
    expect_error(setup(weight, 1L), "vectors must contain values between 0 and 1")
  }
  for(weight in list(NA_real_, NaN, -1, numeric(), 2+1i)){
    expect_error(setup(weight), "blendWeight")
  }
  ml <- mmes(Yield~Env, data=DT_example, control=mmes_control(update=list(blendWeight=3), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=3L), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(ml$blendWeight, exp(seq(log(1), log(0.03), length.out=3L)))
})

test_that("pxem selects expansion without changing pure AI behavior", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(update=list(recovery="pxem", blendWeight=3), convergence=list(maxIters=4L), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(attr(setup$covStruct, "blend"), "pxem")
  expect_equal(setup$blend, "pxem")
  expect_equal(setup$blendWeight, c(1, sqrt(0.03), 0.03, 0.03))
  fit <- function(blend) mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(update=list(recovery=blend, blendWeight=0), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  ai <- fit("em")
  pxAI <- fit("pxem")
  for(component in c("llik", "bu", "covPar", "InfMat")){
    expect_equal(pxAI[[component]], ai[[component]])
  }
  expect_equal(pxAI$inputArgs$blend, "pxem")
  ml <- mmes(Yield~Env, data=DT_example, control=mmes_control(update=list(recovery="pxem", blendWeight=3), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=2L), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(isTRUE(attr(ml$covStruct, "directionBlend")))
})

test_that("positive EM-direction blending rejects unsupported engines", {
  data(DT_example)
  fit <- function(control=mmes_control(convergence=list(maxIters=2L),
                                       update=list(blendWeight=1)),
                  family=stats::gaussian()){
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, family=family,
      control=control)
  }

  expect_error(fit(mmes_control(solve=list(method="pcg"),
    update=list(blendWeight=1), convergence=list(maxIters=2L))), "update\\$blendWeight")
  expect_error(fit(family=binomial()), "PQL")
  expect_error(fit(mmes_control(proposal=list(extrapolation="aitken"),
    update=list(blendWeight=1), convergence=list(maxIters=2L))),
    "Aitken")
  expect_error(fit(mmes_control(representation=list(factorScores="fixed-shape"),
    update=list(blendWeight=1), convergence=list(maxIters=2L))),
    "representation\\$factorScores")
  expect_error(fit(mmes_control(representation=list(equations="saddle-point"),
    update=list(blendWeight=1), convergence=list(maxIters=2L))),
    "saddle-point")

  ai <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=2L), solve=list(method="pcg"), update=list(blendWeight=0), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(ai$blendWeight, c(0, 0))
  expect_false(isTRUE(attr(ai$covStruct, "directionBlend")))
})

test_that("unsupported EM modes default to pure AI only when weights are omitted", {
  data(DT_example)
  modes <- list(
    list(solve=list(method="pcg")),
    list(proposal=list(extrapolation="aitken")),
    list(representation=list(equations="saddle-point")))
  for(mode in modes){
    control <- utils::modifyList(unclass(mmes_control(
      convergence=list(maxIters=3L), output=list(setup=TRUE),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))), mode,
      keep.null=TRUE)
    setup <- mmes(Yield~Env, fixed=Yield~Env, random=~Name, data=DT_example,
      control=do.call(mmes_control, control))
    expect_equal(setup$blendWeight, rep(0, 3L))
    expect_false(isTRUE(attr(setup$covStruct, "directionBlend")))
  }
})
