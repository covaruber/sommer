test_that("fixed-shape FA/RR augmentation matches the marginal MME", {
  set.seed(82)
  data <- expand.grid(env=factor(letters[1:4]), id=factor(seq_len(7L)))
  data$y <- rnorm(nrow(data))
  for(model in c("fa", "rr")){
    shape <- if(model == "fa") fam(data$env, 2L, fixed=rep(TRUE, 10L)) else
      rrm(data$env, 2L, fixed=rep(TRUE, 7L))
    fit <- function(mode){
      mmes(y~env, random=~vsm(shape, ism(id)), rcov=~units,
        data=data, nIters=60, tolParConvLL=1e-9, tolParConvNorm=1e-9,
        computeCi=0, solver="ldlt", verbose=FALSE,
        dateWarning=FALSE, factorScoreAugmentation=mode)
    }
    marginal <- fit("none")
    augmented <- fit("fixed-shape")
    expect_equal(tail(as.numeric(augmented$llik), 1L), tail(as.numeric(marginal$llik), 1L), tolerance=1e-8)
    expect_equal(unlist(augmented$theta[[1L]]), unlist(marginal$theta[[1L]]), tolerance=1e-5)
    expect_equal(as.numeric(augmented$bu), as.numeric(marginal$bu), tolerance=1e-5)
    expect_equal(as.numeric(fitted(augmented)), as.numeric(fitted(marginal)), tolerance=1e-5)
    expect_equal(dim(augmented$C)[1L], length(augmented$bu) +
      (if(model == "fa") 2L else 2L) * nlevels(data$id))
    expect_identical(augmented$CRepresentation, "factor-score-augmented")
    expect_error(mmes(y~env, random=~vsm(if(model == "fa") fam(env, 2L) else
      rrm(env, 2L), ism(id)), rcov=~units, data=data, nIters=2,
      factorScoreAugmentation="fixed-shape"), "requires every")
  }
})

test_that("augmented prediction contrasts map to the original-effect coefficients", {
  set.seed(83)
  data <- expand.grid(env=factor(letters[1:3]), id=factor(seq_len(5L)))
  data$y <- rnorm(nrow(data))
  fit <- mmes(y~env, random=~vsm(fam(env, 1L, fixed=rep(TRUE, 5L)), ism(id)),
    rcov=~units, data=data, nIters=4, computeCi=0, solver="ldlt",
    verbose=FALSE, dateWarning=FALSE, factorScoreAugmentation="fixed-shape")
  contrast <- Matrix::Diagonal(nrow(fit$bu))
  dimnames(contrast) <- list(rownames(fit$bu), rownames(fit$bu))
  augmentedContrast <- sommer:::.mmes_engine_contrast(fit, contrast)
  expect_equal(as.numeric(contrast %*% fit$bu),
               as.numeric(augmentedContrast %*% fit$bu_augmented), tolerance=1e-10)
  prediction <- predict(fit, D=as.matrix(contrast))
  expect_equal(as.numeric(prediction$pvals$predicted.value),
               as.numeric(fit$bu), tolerance=1e-10)
})

test_that("fixed-shape augmentation uses the CHOLMOD dense-Gu Schur path", {
  set.seed(84)
  nId <- 14L
  ids <- paste0("id", seq_len(nId))
  marker <- matrix(sample(c(-1, 0, 1), nId * 70L, replace=TRUE), nrow=nId)
  relationship <- tcrossprod(scale(marker)) / ncol(marker) + diag(0.25, nId)
  precision <- solve(relationship)
  dimnames(precision) <- list(ids, ids)
  attr(precision, "inverse") <- TRUE
  data <- expand.grid(env=factor(paste0("E", seq_len(6L))),
                      id=factor(ids, levels=ids))
  set.seed(85)
  data$y <- rnorm(nrow(data))
  shape <- fam(data$env, 2L, fixed=rep(TRUE, 16L))
  fit <- function(augmentation, solver){
    mmes(y~env, random=~vsm(shape, ism(id), Gu=precision),
      rcov=~units, data=data, nIters=60, tolParConvLL=1e-9,
      tolParConvNorm=1e-9, computeCi=0, solver=solver,
      verbose=FALSE, dateWarning=FALSE,
      factorScoreAugmentation=augmentation)
  }
  marginal <- fit("none", "cholmod")
  augmented <- fit("fixed-shape", "auto")
  expect_true(marginal$engineDiagnostics$blockSchurActive)
  expect_equal(marginal$engineDiagnostics$blockSchurGroups, 1L)
  expect_equal(marginal$engineDiagnostics$blockSchurBorder, 6L)
  expect_identical(augmented$solver, "cholmod")
  expect_true(augmented$engineDiagnostics$blockSchurActive)
  expect_equal(augmented$engineDiagnostics$blockSchurGroups, 6L)
  expect_equal(augmented$engineDiagnostics$blockSchurBorder, 34L)
  expect_identical(as.character(augmented$solver), "cholmod")
  expect_gt(Matrix::nnzero(precision) / (nId * nId), 0.2)
  expect_equal(tail(as.numeric(augmented$llik), 1L),
               tail(as.numeric(marginal$llik), 1L), tolerance=2e-6)
  expect_equal(as.numeric(augmented$bu), as.numeric(marginal$bu), tolerance=2e-5)
  expect_equal(as.numeric(fitted(augmented)), as.numeric(fitted(marginal)),
               tolerance=2e-5)
  expect_gt(nrow(augmented$C), nrow(marginal$C))
  expect_identical(augmented$CRepresentation, "factor-score-augmented")

  profiled <- mmes(y~env, random=~vsm(fam(env, 2L),
      ism(id), Gu=precision), rcov=~units, data=data, nIters=3,
    computeCi=0, solver="auto", verbose=FALSE, dateWarning=FALSE,
    factorScoreAugmentation="profile")
  expect_identical(as.character(profiled$solver), "cholmod")
  expect_true(profiled$engineDiagnostics$blockSchurActive)
  expect_gt(profiled$factorScoreProfile$profileEvaluations, 1L)
})

test_that("profile mode optimizes free FA shapes using augmented REML", {
  set.seed(9)
  data <- expand.grid(env=factor(letters[1:3]), id=factor(seq_len(8L)))
  data$y <- rnorm(nrow(data))
  initial <- mmes(y~env, random=~vsm(fam(env, 1L), ism(id)),
    rcov=~units, data=data, returnParam=TRUE, verbose=FALSE,
    dateWarning=FALSE)
  profile <- mmes(y~env, random=~vsm(fam(env, 1L), ism(id)),
    rcov=~units, data=data, nIters=5, computeCi=0, solver="ldlt",
    verbose=FALSE, dateWarning=FALSE, factorScoreAugmentation="profile")
  expect_identical(profile$factorScoreAugmentation, "profile")
  expect_identical(profile$CRepresentation, "factor-score-augmented")
  expect_gt(profile$factorScoreProfile$profileEvaluations, 1L)
  expect_gt(profile$factorScoreProfile$warmStartedEvaluations, 0L)
  expect_equal(profile$factorScoreProfile$totalSymbolicAnalyses, 1L)
  expect_true(profile$factorScoreProfile$objective <=
              profile$factorScoreProfile$startObjective)
  expect_gt(max(abs(profile$factorScoreProfile$parameters -
                    initial$covStruct[[1L]]$par[-1L])), 1e-5)
  expect_true(all(is.finite(covparams_mmes(profile)$estimate)))
  expect_true(anyNA(profile$theta_se))
  if(profile$factorScoreProfile$convergence != 0L){
    expect_false(profile$convergence)
  }
  expect_error(mmes(y~env, random=~vsm(fam(env, 1L), ism(id)),
    rcov=~units, data=data, nIters=2, factorScoreAugmentation="profile",
    solver="pcg", verbose=FALSE), "requires Gaussian REML")

  rrProfile <- mmes(y~env, random=~vsm(rrm(env, 1L), ism(id)),
    rcov=~units, data=data, nIters=5, computeCi=0, solver="ldlt",
    verbose=FALSE, dateWarning=FALSE, factorScoreAugmentation="profile")
  expect_identical(rrProfile$factorScoreAugmentation, "profile")
  expect_gt(rrProfile$factorScoreProfile$profileEvaluations, 1L)
  expect_lt(rrProfile$factorScoreProfile$objective,
            rrProfile$factorScoreProfile$startObjective)
})

test_that("profile mode updates free FA parameters through augmented REML fits", {
  set.seed(86)
  data <- expand.grid(env=factor(letters[1:3]), id=factor(seq_len(8L)))
  data$y <- rnorm(nrow(data))
  profile <- mmes(y~env, random=~vsm(fam(env, 1L), ism(id)),
    rcov=~units, data=data, nIters=5, computeCi=0, solver="ldlt",
    verbose=FALSE, dateWarning=FALSE, factorScoreAugmentation="profile")
  expect_identical(profile$factorScoreAugmentation, "profile")
  expect_identical(profile$CRepresentation, "factor-score-augmented")
  expect_gt(profile$factorScoreProfile$profileEvaluations, 1L)
  expect_gt(profile$factorScoreProfile$warmStartedEvaluations, 0L)
  expect_lt(profile$factorScoreProfile$objective,
            profile$factorScoreProfile$startObjective)
  expect_true(all(is.finite(covparams_mmes(profile)$estimate)))
  expect_true(anyNA(profile$theta_se))
  expect_error(mmes(y~env, random=~vsm(fam(env, 1L), ism(id)),
    rcov=~units, data=data, nIters=2, factorScoreAugmentation="profile",
    solver="pcg", verbose=FALSE), "Gaussian REML Henderson")
})