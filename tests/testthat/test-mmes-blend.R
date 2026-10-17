test_that("blend weights validate and recycle to the iteration count", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
    nIters=5L, blendWeight=c(0, 0.4), returnParam=TRUE,
    verbose=FALSE, dateWarning=FALSE)

  expect_equal(setup$blendWeight, c(0, 0.4, 0, 0.4, 0))
  expect_true(isTRUE(attr(setup$covStruct, "directionBlend")))

  expect_error(mmes(Yield~Env, data=DT_example, blendWeight=1.1),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, blendWeight=Inf),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, blendWeight="0.5"),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, nIters=1L, blendWeight=c(1, NA)),
    "blendWeight")
  expect_error(mmes(Yield~Env, data=DT_example, nIters=1e20), "positive integer")
})

test_that("scalar blend counts generate capped exponential decay schedules", {
  data(DT_example)
  setup <- function(weight, iterations=6L){
    mmes(Yield~Env, random=~Name, data=DT_example, nIters=iterations,
      blendWeight=weight, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
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
  ml <- mmes(Yield~Env, data=DT_example, blendWeight=3, REML=FALSE,
    nIters=3L, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
  expect_equal(ml$blendWeight, exp(seq(log(1), log(0.03), length.out=3L)))
})

test_that("fitted objects retain normalized blend settings for rebuilds", {
  data(DT_example)
  fit <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
    nIters=3L, blendWeight=0, verbose=FALSE, dateWarning=FALSE)
  expect_equal(fit$inputArgs$blend, "em")
  expect_equal(fit$inputArgs$blendWeight, rep(0, 3L))
})

test_that("deprecated emWeight aliases blendWeight and rejects collisions", {
  data(DT_example)
  expect_warning(
    setup <- mmes(Yield~Env, random=~Name, data=DT_example,
      nIters=3L, emWeight=c(0.2, 0.6), returnParam=TRUE,
      verbose=FALSE, dateWarning=FALSE),
    "deprecated"
  )
  expect_equal(setup$blendWeight, c(0.2, 0.6, 0.2))
  expect_error(mmes(Yield~Env, data=DT_example,
    blendWeight=0, emWeight=0), "Supply only blendWeight")
})

test_that("pxem selects expansion without changing pure AI behavior", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, data=DT_example, blend="pxem",
    blendWeight=3, nIters=4L, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
  expect_equal(attr(setup$covStruct, "blend"), "pxem")
  expect_equal(setup$blend, "pxem")
  expect_equal(setup$blendWeight, c(1, sqrt(0.03), 0.03, 0.03))
  fit <- function(blend) mmes(Yield~Env, random=~Name, data=DT_example,
    blend=blend, blendWeight=0, nIters=3L, verbose=FALSE, dateWarning=FALSE)
  ai <- fit("em")
  pxAI <- fit("pxem")
  for(component in c("llik", "bu", "covPar", "InfMat")){
    expect_equal(pxAI[[component]], ai[[component]])
  }
  expect_equal(pxAI$inputArgs$blend, "pxem")
  ml <- mmes(Yield~Env, data=DT_example, blend="pxem", blendWeight=3,
    REML=FALSE, nIters=2L, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
  expect_true(isTRUE(attr(ml$covStruct, "directionBlend")))
})

test_that("solveOnly ignores unused positive blend weights", {
  data(DT_example)
  solved <- mmes(Yield~Env,
    random=~vsm(ism(Name), sigma2=2, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=1, fixedSigma2=TRUE),
    data=DT_example, solveOnly=TRUE, solver="pcg",
    verbose=FALSE, dateWarning=FALSE)
  expect_s3_class(solved, "mmesSolve")
})

test_that("positive EM-direction blending rejects unsupported engines", {
  data(DT_example)
  fit <- function(...) {
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
      nIters=2L, blendWeight=1, ...)
  }

  expect_error(fit(solver="pcg"), "solver='pcg'")
  expect_error(fit(family=binomial()), "PQL")
  expect_error(fit(acceleration="aitken"), "acceleration")
  expect_error(fit(factorScoreAugmentation="fixed-shape"),
    "factorScoreAugmentation")
  expect_error(fit(saddlePoint=TRUE), "saddlePoint='force'")

  ai <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
    nIters=2L, solver="pcg", blendWeight=0, returnParam=TRUE,
    verbose=FALSE, dateWarning=FALSE)
  expect_equal(ai$blendWeight, c(0, 0))
  expect_false(isTRUE(attr(ai$covStruct, "directionBlend")))
})

test_that("unsupported EM modes default to pure AI only when weights are omitted", {
  data(DT_example)
  modes <- list(list(solver="pcg"),
    list(acceleration="aitken"), list(saddlePoint=TRUE))
  for(mode in modes){
    setup <- do.call(mmes, c(list(fixed=Yield~Env, random=~Name, data=DT_example,
      nIters=3L, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE), mode))
    expect_equal(setup$blendWeight, rep(0, 3L))
    expect_false(isTRUE(attr(setup$covStruct, "directionBlend")))
  }
})
