test_that("fixed-shape heterogeneous/homogeneous FA augmentation matches the marginal MME", {
  set.seed(82)
  data <- expand.grid(env=factor(letters[1:4]), id=factor(seq_len(7L)))
  data$y <- rnorm(nrow(data))
  for(model in c("fa", "rr")){
    shape <- if(model == "fa") fam(data$env, 2L, fixed=rep(TRUE, 10L)) else
      fam(data$env, 2L, fixed=rep(TRUE, 7L), variance="homogeneous")
    fit <- function(mode){
      mmes(y~env, random=~vsm(shape, ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=60, likelihoodTol=1e-12, parameterTol=1e-12), output=list(inverse=0), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE), representation=list(factorScores=mode)))
    }
    marginal <- fit("none")
    augmented <- fit("fixed-shape")
    expect_true(marginal$convergence)
    expect_true(augmented$convergence)
    expect_equal(tail(as.numeric(augmented$llik), 1L), tail(as.numeric(marginal$llik), 1L), tolerance=1e-8)
    expect_equal(unlist(augmented$theta[[1L]]), unlist(marginal$theta[[1L]]), tolerance=1e-5)
    expect_equal(as.numeric(augmented$bu), as.numeric(marginal$bu), tolerance=1e-5)
    expect_equal(as.numeric(fitted(augmented)), as.numeric(fitted(marginal)), tolerance=1e-5)
    expect_equal(dim(augmented$C)[1L], length(augmented$bu) +
      (if(model == "fa") 2L else 2L) * nlevels(data$id))
    expect_identical(augmented$CRepresentation, "factor-score-augmented")
    expect_error(mmes(y~env, random=~vsm(if(model == "fa") fam(env, 2L) else
      fam(env, 2L, variance="homogeneous"), ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=2), representation=list(factorScores="fixed-shape"))), "requires all covariance-shape")
  }
})

test_that("compound-symmetry augmentation matches marginal fits", {
  set.seed(86)
  q <- 4L
  nId <- 35L
  data <- expand.grid(env=factor(letters[1:q]), id=factor(seq_len(nId)))
  shapes <- list(
    csm(data$env, rho=0.25, fixed=TRUE),
    csm(data$env, rho=0.2, variance="heterogeneous",
        values=c(1, 1.4, 0.8, 1.2), fixed=rep(TRUE, 4L))
  )

  covariances <- list(
    {K <- matrix(0.25, q, q); diag(K) <- 1; K},
    {v <- c(1, 1.4, 0.8, 1.2); K <- outer(sqrt(v), sqrt(v), `*`) * 0.2; diag(K) <- v; K}
  )
  for(index in seq_along(shapes)){
    shape <- shapes[[index]]
    latent <- t(chol(covariances[[index]])) %*%
      matrix(rnorm(q * nId), q, nId) * sqrt(0.5)
    data$y <- 2 + as.numeric(data$env) +
      latent[cbind(as.integer(data$env), as.integer(data$id))] +
      rnorm(nrow(data), sd=sqrt(0.5))
    fit <- function(mode){
      mmes(y~env, random=~vsm(shape, ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=100, likelihoodTol=1e-8, parameterTol=1e-8), output=list(inverse=0), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE), representation=list(factorScores=mode)))
    }
    marginal <- fit("none")
    augmented <- fit("fixed-shape")
    expect_true(marginal$convergence)
    expect_true(augmented$convergence)
    expect_equal(tail(as.numeric(augmented$llik), 1L),
                 tail(as.numeric(marginal$llik), 1L), tolerance=1e-8)
    expect_equal(unlist(augmented$theta[[1L]]),
                 unlist(marginal$theta[[1L]]), tolerance=1e-4)
    expect_equal(as.numeric(augmented$bu), as.numeric(marginal$bu),
                 tolerance=1e-4)
    expect_identical(augmented$CRepresentation, "factor-score-augmented")
  }
})

test_that("augmented prediction contrasts map to the original-effect coefficients", {
  set.seed(83)
  data <- expand.grid(env=factor(letters[1:3]), id=factor(seq_len(5L)))
  data$y <- rnorm(nrow(data))
  fit <- mmes(y~env, random=~vsm(fam(env, 1L, fixed=rep(TRUE, 5L)), ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=4), output=list(inverse=0), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE), representation=list(factorScores="fixed-shape")))
  contrast <- Matrix::Diagonal(nrow(fit$bu))
  dimnames(contrast) <- list(rownames(fit$bu), rownames(fit$bu))
  augmentedContrast <- sommer:::.mmes_engine_contrast(fit, contrast)
  expect_equal(as.numeric(contrast %*% fit$bu),
               as.numeric(augmentedContrast %*% fit$bu_augmented), tolerance=1e-10)
  prediction <- predict(fit, D=as.matrix(contrast))
  expect_equal(as.numeric(prediction$pvals$predicted.value),
               as.numeric(fit$bu), tolerance=1e-10)
})

test_that("profile mode optimizes free FA shapes using augmented REML", {
  set.seed(9)
  data <- expand.grid(env=factor(letters[1:2]), id=factor(seq_len(8L)))
  data$y <- rnorm(nrow(data))
  shape <- fam(data$env, 1L, fixed=c(TRUE, TRUE, FALSE))
  initial <- mmes(y~env, random=~vsm(shape, ism(id)), rcov=~units, data=data, control=mmes_control(output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  profile <- mmes(y~env, random=~vsm(shape, ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=1), output=list(inverse=0), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE), representation=list(factorScores="profile")))
  expect_identical(profile$factorScoreAugmentation, "profile")
  expect_identical(profile$CRepresentation, "factor-score-augmented")
  expect_gt(profile$factorScoreProfile$profileEvaluations, 1L)
  expect_gt(profile$factorScoreProfile$warmStartedEvaluations, 0L)
  expect_equal(profile$factorScoreProfile$totalSymbolicAnalyses, 1L)
  expect_true(profile$factorScoreProfile$objective <=
              profile$factorScoreProfile$startObjective)
  expect_gt(max(abs(profile$factorScoreProfile$parameters -
                    initial$covStruct[[1L]]$par[-1L][initial$covStruct[[1L]]$free[-1L]])), 1e-5)
  expect_true(all(is.finite(covparams_mmes(profile)$estimate)))
  expect_true(anyNA(profile$theta_se))
  if(profile$factorScoreProfile$convergence != 0L){
    expect_false(profile$convergence)
  }
  expect_error(mmes(y~env, random=~vsm(fam(env, 1L), ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=2), representation=list(factorScores="profile"), solve=list(method="pcg"), diagnostics=list(verbose=FALSE))), "representation\\$factorScores")

  rrProfile <- mmes(y~env, random=~vsm(fam(env, 1L, fixed=c(TRUE, FALSE),
    variance="homogeneous"), ism(id)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=1), output=list(inverse=0), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE), representation=list(factorScores="profile")))
  expect_identical(rrProfile$factorScoreAugmentation, "profile")
  expect_gt(rrProfile$factorScoreProfile$profileEvaluations, 1L)
  expect_lt(rrProfile$factorScoreProfile$objective,
            rrProfile$factorScoreProfile$startObjective)
})
