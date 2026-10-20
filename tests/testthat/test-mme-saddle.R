.saddle_dataset <- function(){
  set.seed(892)
  dataset <- expand.grid(id=factor(seq_len(8L)), level=factor(seq_len(3L)),
                         replicate=seq_len(3L))
  dataset$response <- as.integer(dataset$level) / 4 +
    rep(rnorm(8L), 9L) + rnorm(nrow(dataset), sd=0.6)
  dataset
}

.saddle_at <- function(setup, saddle=FALSE, mode=0L){
  attr(setup$covStruct, "saddlePoint") <- if(saddle) "force" else "off"
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv,
    setup$covStruct, 0, 1, FALSE, computeCi=mode, solver=setup$solver,
    reml=setup$REML)
}

test_that("augmented numerical identities pass standalone native checks", {
  .native_run("test-mme-saddle", c("mme_saddle.h", "mme_exact.h"), eigen=TRUE)
})

test_that("saddle controls are opt-in and reject unsupported routes", {
  expect_equal(mmes_control()$representation$equations, "normal")
  expect_equal(mmes_control(representation=list(equations="saddle-point"))$
    representation$equations, "saddle-point")
  dataset <- .saddle_dataset()
  saddleControl <- mmes_control(representation=list(equations="saddle-point"),
    update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  model <- function(control=list(), family=stats::gaussian(), internal=NULL){
    merged <- utils::modifyList(unclass(saddleControl), control, keep.null=TRUE)
    mmes(response~level, random=~id, data=dataset, family=family,
      control=do.call(mmes_control, merged), .internal=internal)
  }
  expect_error(model(list(solve=list(method="pcg")), family=stats::gaussian()), "representation\\$equations")
  expect_error(model(list(formulation="direct")), "representation\\$equations")
  expect_error(model(list(task=list(mode="solve"))), "representation\\$equations")
  expect_error(model(list(curvature=list(method="hybrid-ai"))), "representation\\$equations")
  expect_error(model(family=poisson()), "representation\\$equations")
  expect_error(model(list(representation=list(factorScores="fixed-shape"))),
    "representation\\$equations")
  expect_error(model(internal=list(pqlInner=TRUE)), "representation\\$equations")
  setup <- mmes(response~level, random=~id, data=dataset, control=mmes_control(output=list(setup=TRUE), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  reference <- .saddle_at(setup)
  attr(setup$covStruct, "saddlePoint") <- "auto"
  auto <- sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv,
    setup$covStruct, 0, 1, FALSE, solver=setup$solver)
  expect_identical(auto$bu, reference$bu)
  expect_null(auto$engineDiagnostics$saddlePoint)
  expect_error(mmes(response~level, random=~vsm(rrm(level, 1L), ism(id)), data=dataset,
    control=saddleControl),
    "reduced rank")
})

test_that("saddle likelihood score and AI match independent dense marginal formulas", {
  dataset <- .saddle_dataset()
  for(reml in c(TRUE, FALSE)){
    setup <- mmes(response~level, random=~id, data=dataset, control=mmes_control(output=list(setup=TRUE), task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    fit <- .saddle_at(setup, TRUE, 2L)
    X <- as.matrix(setup$X)
    Z <- as.matrix(setup$Z[[1L]])
    y <- as.numeric(fit$yWorking)
    scale <- fit$Cscale
    su <- as.numeric(fit$covPar[[1L]]) / scale
    se <- as.numeric(fit$covPar[[2L]]) / scale
    derivatives <- list(su * tcrossprod(Z), se * diag(nrow(X)))
    V <- derivatives[[1L]] + derivatives[[2L]]
    Vi <- solve(V)
    fixed <- crossprod(X, Vi %*% X)
    P <- Vi - Vi %*% X %*% solve(fixed, t(X) %*% Vi)
    Py <- P %*% y
    quadratic <- as.numeric(crossprod(y, Py))
    likelihood <- -0.5 * (as.numeric(determinant(V, logarithm=TRUE)$modulus) +
      if(reml) as.numeric(determinant(fixed, logarithm=TRUE)$modulus) else 0)
    likelihood <- likelihood - 0.5 * quadratic
    expect_equal(as.numeric(fit$llik), likelihood, tolerance=1e-8)
    score <- vapply(derivatives, function(derivative){
      as.numeric(crossprod(Py, derivative %*% Py)) -
        sum(diag((if(reml) P else Vi) %*% derivative))
    }, numeric(1))
    working <- do.call(cbind, lapply(derivatives, function(derivative) derivative %*% Py))
    twiceAI <- crossprod(working, P %*% working)
    expect_equal(as.numeric(fit$dLu), -score, tolerance=1e-7)
    expect_equal(unname(fit$InfMat), unname(twiceAI), tolerance=1e-7)
  }
})

test_that("weighted and repeated residual blocks retain exact saddle derivatives", {
  dataset <- .saddle_dataset()
  for(weight in c("none", "diagonal", "nondiagonal")){
    for(reml in c(TRUE, FALSE)){
      arguments <- list(fixed=response~level, random=~id,
        rcov=~vsm(ar1m(level, rho=0.3), ism(units)), data=dataset,
        control=mmes_control(task=list(objective=if(reml) "reml" else "ml"),
          update=list(blendWeight=0), output=list(setup=TRUE),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      if(weight != "none"){
        arguments$W <- if(weight == "diagonal") diag(seq(0.7, 1.3, length.out=nrow(dataset))) else
          Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
            diagonals=list(rep(0.05, nrow(dataset)-1L), rep(1.2, nrow(dataset)),
                           rep(0.05, nrow(dataset)-1L)))
      }
      setup <- do.call(mmes, arguments)
      reference <- .saddle_at(setup)
      candidate <- .saddle_at(setup, TRUE)
      for(field in c("llik", "bu", "C", "dLu", "InfMat", "theta_se")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-7,
          info=paste(weight, reml, field))
      }
    }
  }
})

test_that("saddle convergence and coefficient inverse modes preserve post-fit consumers", {
  dataset <- .saddle_dataset()
  model <- function(saddle, mode, reml=TRUE){
    mmes(response~level, random=~id, data=dataset, control=mmes_control(
      representation=list(equations=if(saddle) "saddle-point" else "normal"),
      output=list(inverse=mode), task=list(objective=if(reml) "reml" else "ml"),
      update=list(blendWeight=0),
      convergence=list(maxIters=60L, likelihoodTol=1e-10, parameterTol=1e-8),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  for(reml in c(TRUE, FALSE)) for(mode in 0:2){
    reference <- model(FALSE, mode, reml)
    candidate <- model(TRUE, mode, reml)
    expect_equal(candidate$convergence, reference$convergence)
    expect_true(candidate$convergence)
    for(field in c("llik", "covPar", "bu", "theta_se", "Ci", "uPevList")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-6,
        info=paste(reml, mode, field))
    }
    expect_equal(candidate$engineDiagnostics$saddlePoint$compatibilityAssemblies, 1L)
    expect_lte(candidate$engineDiagnostics$saddlePoint$symbolicAnalyses, if(reml) 1 else 2)
    D <- Matrix::Diagonal(length(candidate$bu))[1:3, , drop=FALSE]
    expect_equal(sommer:::predict_mmes_vcov_cpp(candidate, D),
                 sommer:::predict_mmes_vcov_cpp(reference, D), tolerance=1e-6)
    saved <- tempfile(fileext=".rds")
    saveRDS(candidate, saved)
    restored <- readRDS(saved)
    unlink(saved)
    expect_equal(sommer:::predict_mmes_vcov_cpp(restored, D),
                 sommer:::predict_mmes_vcov_cpp(candidate, D), tolerance=1e-12)
  }
})

test_that("saddle supports fixed and tied parameters and models without random effects", {
  dataset <- .saddle_dataset()
  for(reml in c(TRUE, FALSE)){
    setup <- mmes(response~level, data=dataset, control=mmes_control(task=list(objective=if(reml) "reml" else "ml"), output=list(setup=TRUE), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    reference <- .saddle_at(setup)
    candidate <- .saddle_at(setup, TRUE)
    for(field in c("llik", "bu", "dLu", "InfMat", "theta_se")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
    }
  }
  for(tied in c(FALSE, TRUE)){
    model <- function(saddle){
      mmes(response~level, random=~vsm(dsm(level), ism(id)), data=dataset, vcc=if(tied) data.frame(parameter=c(2L, 3L), group=1L) else NULL, control=mmes_control(representation=list(equations=if(saddle) "saddle-point" else "normal"), update=list(blendWeight=0), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    if(!tied){
      # Fixed descriptor coordinates still contribute exact final information.
      setup <- mmes(response~level, random=~vsm(ar1m(level, rho=0.3, fixed=TRUE), ism(id)), data=dataset, control=mmes_control(update=list(blendWeight=0), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      reference <- .saddle_at(setup)
      candidate <- .saddle_at(setup, TRUE)
    }else{
      reference <- model(FALSE)
      candidate <- model(TRUE)
    }
    for(field in c("llik", "bu", "dLu", "InfMat", "theta_se")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-7)
    }
  }
})

test_that("saddle rejects compact grids rather than materializing hidden covariance fill", {
  dataset <- expand.grid(trial=factor(1:2), row=factor(1:3), column=factor(1:3))
  dataset <- dataset[-1L, ]
  dataset$response <- sin(seq_len(nrow(dataset)))
  expect_error(mmes(response~trial, rcov=~vsm(ar1m(row), ar1m(column), ism(units)), data=dataset, control=mmes_control(representation=list(equations="saddle-point"), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "grid/ragged")
})

test_that("saddle prediction and emmeans preserve their public coefficient contract", {
  dataset <- .saddle_dataset()
  model <- function(saddle){
    mmes(response~level, random=~id, data=dataset, control=mmes_control(representation=list(equations=if(saddle) "saddle-point" else "normal"), convergence=list(maxIters=40L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  reference <- model(FALSE)
  candidate <- model(TRUE)
  for(classify in c("level", "id")){
    expect_equal(predict(candidate, D=classify)$pvals, predict(reference, D=classify)$pvals,
      tolerance=1e-7)
  }
  skip_if_not_installed("emmeans")
  expect_equal(as.data.frame(emmeans::emmeans(candidate, ~level)),
    as.data.frame(emmeans::emmeans(reference, ~level)), tolerance=1e-7)
})
