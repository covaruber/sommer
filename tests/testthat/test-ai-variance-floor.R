test_that("AI reaches variance floors without damping unrelated or tied coordinates", {
  helper <- .native_ai_recovery()
  lower <- c(-10, -Inf, -10, -Inf)
  point <- c(-4, 0, -3, 2)
  information <- diag(c(1e-10, 2, 1, 1))
  score <- c(-1e-3, 1, -20, 0)
  tangent <- diag(4)[, 1:3]
  for(reduced in c(FALSE, TRUE)){
    step <- helper$feasible_ai_probe(information, score, tangent, point, lower,
      rep(1, 4), reduced)
    expect_true(step$reliable)
    expect_equal(as.numeric(point + step$direction), c(-10, 0.5, -10, 2))
    expect_equal(step$ridge, 0)
    expect_equal(as.numeric(step$cappedDirection), c(0, 0.5, 0, 0))
    expect_equal(ncol(step$nextTangent), 1L)
  }
  tied <- cbind(c(1, 0, 1, 0), c(0, 1, 0, 0))
  step <- helper$feasible_ai_probe(information, score, tied, point, lower, rep(1, 4))
  expect_true(step$reliable)
  expect_equal(step$direction[1], step$direction[3])
  expect_gte(min(point[c(1, 3)] + step$direction[c(1, 3)] - lower[c(1, 3)]), -1e-10)
  expect_equal(step$direction[4], 0)
  expect_equal(ncol(step$nextTangent), 1L)

  floorPoint <- c(-10, 0, -10, 2)
  outward <- helper$feasible_ai_probe(information, score, tangent, floorPoint, lower, rep(1, 4))
  expect_equal(as.numeric(outward$direction), c(0, 0.5, 0, 0))
  inward <- helper$feasible_ai_probe(diag(4), c(1, 0, 1, 0), tangent,
    floorPoint, lower, rep(1, 4))
  expect_equal(as.numeric(inward$direction), c(1, 0, 1, 0))
})

test_that("small undamped AI trust invokes a regularized direction", {
  helper <- .native_ai_recovery()
  step <- helper$feasible_ai_probe(diag(2), c(100, 1), diag(2),
    c(0, 0), c(-Inf, -Inf), c(1, 1))
  expect_true(step$reliable)
  expect_equal(step$attempts, 1L)
  expect_gt(step$ridge, 0)
  expect_lte(max(abs(step$direction)), 1)
  healthy <- helper$feasible_ai_probe(diag(2), c(0.5, 1), diag(2),
    c(0, 0), c(-Inf, -Inf), c(1, 1))
  expect_equal(healthy$ridge, 0)
  expect_equal(as.numeric(healthy$direction), c(0.5, 1))
})

test_that("correlated Newton overshoots do not pin marginally supported variances", {
  helper <- .native_ai_recovery()
  information <- matrix(c(1, 0.99, 0.99, 1), 2)
  step <- helper$feasible_ai_probe(information, c(-1, 1), diag(2),
    c(0, 0), c(-10, -Inf), c(1, 1))
  expect_true(step$reliable)
  expect_gt(step$ridge, 0)
  expect_gt(step$direction[1], -10)
})

test_that("diagonal variance floors converge consistently on direct and MME paths", {
  dataset <- expand.grid(id=factor(seq_len(12)), environment=factor(letters[1:2]),
    replicate=seq_len(4))
  dataset$y <- (dataset$environment == "a") * 2 * sin(as.numeric(dataset$id)) +
    c(-1, -0.4, 0.4, 1)[dataset$replicate]
  for(objective in c("reml", "ml")){
    fits <- lapply(c("cholmod", "ldlt", "direct"), function(solver){
      mmes(y~environment, random=~vsm(dsm(environment), ism(id), sigma2=1),
        data=dataset, control=mmes_control(
          formulation=if(solver == "direct") "direct" else "henderson",
          solve=list(method=if(solver == "direct") "auto" else solver),
          task=list(objective=objective),
          convergence=list(maxIters=20, likelihoodTol=1e-8, parameterTol=1e-8),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    })
    for(fit in fits){
      recovery <- fit$engineDiagnostics$recovery
      expect_true(fit$convergence)
      expect_lte(length(fit$llik), 20)
      expect_equal(recovery$stopReason, "converged")
      expect_lte(recovery$finalProjectedScoreNorm, recovery$scoreTolerance)
      expect_equal(as.numeric(fit$covParWorking[[1]][2]), log(1e-6), tolerance=1e-8)
      expect_true(all(diff(as.numeric(fit$llik)) >= -1e-9))
      expect_equal(unlist(fit$covParWorking), unlist(fits[[1]]$covParWorking), tolerance=1e-7)
    }
  }
})

test_that("latent RR certifies outward diagonal boundary scores with the correct sign", {
  dataset <- expand.grid(id=factor(seq_len(8)), block=factor(seq_len(4)),
    environment=factor(letters[1:3]), replicate=seq_len(2))
  dataset$y <- c(0, 1, 2)[dataset$environment] +
    sin(as.numeric(dataset$id)) * c(1, 0.7, 1.4)[dataset$environment] +
    (dataset$environment == "a") * cos(as.numeric(dataset$block)) +
    c(-0.4, 0.4)[dataset$replicate]
  for(objective in c("reml", "ml")){
    fit <- mmes(y~environment,
      random=~vsm(rrm(environment, 1), ism(id), sigma2=1) +
        vsm(dsm(environment), ism(block), sigma2=1),
      data=dataset, control=mmes_control(task=list(objective=objective),
        convergence=list(maxIters=30), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    recovery <- fit$engineDiagnostics$recovery
    expect_true(fit$convergence)
    expect_lte(length(fit$llik), 20)
    expect_equal(recovery$stopReason, "converged")
    expect_lte(recovery$finalProjectedScoreNorm, recovery$scoreTolerance)
    expect_equal(as.numeric(fit$covParWorking[[1]][2]), 0)
    expect_equal(as.numeric(fit$covParWorking[[2]][2:3]), rep(log(1e-6), 2), tolerance=1e-8)
    expect_true(all(fit$engineDiagnostics$scoreWorking[6:7] < 0))
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-9))
  }
})

test_that("fixed variance ratios below the numerical floor remain exact", {
  dataset <- expand.grid(id=factor(seq_len(8)), environment=factor(letters[1:2]),
    replicate=seq_len(2))
  dataset$y <- sin(as.numeric(dataset$id)) + c(-0.4, 0.4)[dataset$replicate]
  for(formulation in c("henderson", "direct")){
    fit <- mmes(y~environment,
      random=~vsm(dsm(environment, values=c(1, 1e-10), fixed=TRUE), ism(id), sigma2=1),
      data=dataset, control=mmes_control(formulation=formulation,
        convergence=list(maxIters=3), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(as.numeric(fit$covParWorking[[1]][2]), log(1e-10), tolerance=1e-12)
  }
})

test_that("a supported variance ratio can leave its initial floor", {
  dataset <- expand.grid(id=factor(seq_len(12)), environment=factor(letters[1:2]),
    replicate=seq_len(4))
  dataset$y <- 2 * sin(as.numeric(dataset$id)) + c(-1, -0.4, 0.4, 1)[dataset$replicate]
  for(formulation in c("henderson", "direct")){
    fit <- mmes(y~environment,
      random=~vsm(dsm(environment, values=c(1, 1e-6)), ism(id), sigma2=1),
      data=dataset, control=mmes_control(formulation=formulation,
        convergence=list(maxIters=40), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_true(fit$convergence)
    expect_equal(as.numeric(fit$covParWorking[[1]][2]), 0, tolerance=1e-3)
    expect_lte(fit$engineDiagnostics$recovery$finalProjectedScoreNorm,
      fit$engineDiagnostics$recovery$scoreTolerance)
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-9))
  }
})
