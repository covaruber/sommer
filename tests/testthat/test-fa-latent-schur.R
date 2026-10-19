test_that("latent heterogeneous/homogeneous FA Schur preserves AI trajectories and uncertainty", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 4
  precision <- diag(6L) + matrix(0.15, 6L, 6L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  for(model in c("fa", "rr")) for(rank in 1:2) for(inverseMode in 0:2){
    shape <- fam(dataset$environment, rank,
                 variance=if(model == "fa") "heterogeneous" else "homogeneous")
    fitModel <- function(solver){
      mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision) +
          vsm(dsm(environment), ism(replicate)) + vsm(ism(id)), rcov=~vsm(dsm(environment), ism(units)), data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=4L), output=list(inverse=inverseMode), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("cholmod")
    expect_true(candidate$engineDiagnostics$factorSchurActive)
    expect_false(candidate$engineDiagnostics$blockChainActive)
    expect_equal(candidate$engineDiagnostics$blockSchurGroups, 5L)
    expect_equal(dim(candidate$C), dim(reference$C))
    for(component in c("llik", "theta", "covPar", "bu", "InfMat", "theta_se")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-6)
    }
    expect_equal(covparams_mmes_se(candidate)$StdError,
                 covparams_mmes_se(reference)$StdError, tolerance=1e-6)
    if(inverseMode == 2L){ expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-6) }
    if(inverseMode > 0L){ expect_equal(candidate$uPevList, reference$uPevList, tolerance=1e-6) }
  }
})

test_that("latent Schur handles multiple random terms, diagonal weights and missing records", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(7L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  set.seed(651)
  dataset$response <- rnorm(nrow(dataset)) + as.integer(dataset$id) / 5
  dataset <- dataset[-c(2L, 11L, 40L), ]
  precision <- diag(7L) + matrix(0.1, 7L, 7L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  secondShape <- fam(dataset$environment, 1L, variance="homogeneous")
  secondShape$covFactor$free[] <- FALSE
  for(weighted in c(FALSE, TRUE)){
    weightMatrix <- if(weighted) diag(seq(0.7, 1.3, length.out=nrow(dataset))) else diag(nrow(dataset))
    fitModel <- function(solver){
      mmes(response~environment, random=~vsm(fam(environment, 2L), ism(id), Gu=precision) +
          vsm(secondShape, ism(replicate)), rcov=~vsm(dsm(environment), ism(units)), data=dataset, W=weightMatrix, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=4L), output=list(inverse=2L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("cholmod")
    expect_true(candidate$engineDiagnostics$factorSchurActive)
    for(component in c("llik", "theta", "bu", "InfMat", "theta_se")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-6)
    }
    expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-6)
  }
})

test_that("latent Schur rejects incompatible residuals, ML and near-boundary specifics", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  precision <- diag(6L) + matrix(0.1, 6L, 6L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  fitModel <- function(shape, residual=~units, reml=TRUE, solver="cholmod"){
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision), rcov=residual, data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=2L), output=list(inverse=0L), update=list(blendWeight=0), task=list(objective=if(reml) "reml" else "ml"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  shape <- fam(dataset$environment, 1L)
  expect_true(fitModel(shape)$engineDiagnostics$factorSchurActive)
  correlated <- fitModel(shape, ~vsm(ar1m(environment), ism(record)))
  expect_false(correlated$engineDiagnostics$factorSchurActive)
  expect_equal(correlated$llik,
    fitModel(shape, ~vsm(ar1m(environment), ism(record)), solver="ldlt")$llik, tolerance=1e-6)
  overlapping <- shape
  overlapping$Z[1L, 2L] <- 1
  rejected <- fitModel(overlapping)
  expect_false(rejected$engineDiagnostics$factorSchurActive)
  expect_equal(rejected$llik, fitModel(overlapping, solver="ldlt")$llik, tolerance=1e-6)
  previousMultiFactor <- options(sommer.mme.factorSchurMultiFactor=NULL)
  on.exit(options(previousMultiFactor), add=TRUE)
  extraFactorDefault <- mmes(response~environment, random=~vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu=precision), data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=2L), output=list(inverse=0L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(extraFactorDefault$engineDiagnostics$factorSchurActive)
  options(sommer.mme.factorSchurMultiFactor=TRUE)
  extraFactor <- mmes(response~environment, random=~vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu=precision), data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=2L), output=list(inverse=0L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  extraFactorReference <- mmes(response~environment, random=~vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu=precision), data=dataset, control=mmes_control(solve=list(method="ldlt"), convergence=list(maxIters=2L), output=list(inverse=0L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(extraFactor$engineDiagnostics$factorSchurActive)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se")){
    expect_equal(extraFactor[[component]], extraFactorReference[[component]], tolerance=1e-6)
  }
  expect_false(fitModel(shape, reml=FALSE)$engineDiagnostics$factorSchurActive)
  boundary <- fam(dataset$environment, 1L, specific=c(rep(1, 4L), 1e-10))
  boundary$covFactor$free[] <- FALSE
  fallback <- fitModel(boundary)
  expect_false(fallback$engineDiagnostics$factorSchurActive)
  reference <- fitModel(boundary, solver="ldlt")
  expect_equal(fallback$llik, reference$llik, tolerance=1e-6)
  expect_equal(fallback$bu, reference$bu, tolerance=1e-6)
  options(sommer.mme.denseMemoryMB=0.001, sommer.mme.memoryPolicy="limit")
  on.exit(options(sommer.mme.memoryPolicy=NULL), add=TRUE)
  expect_false(fitModel(shape)$engineDiagnostics$factorSchurActive)
})

test_that("latent Schur supports a low-rank axis at any position in longer products", {
  previous <- options(sommer.mme.factorSchurMultiFactor=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), site=factor(1:2),
    environment=factor(letters[1:4]), block=factor(1:3))
  dataset$response <- sin(seq_len(nrow(dataset)) / 5) + as.integer(dataset$id) / 6
  relationship <- diag(6L) + matrix(0.08, 6L, 6L)
  dimnames(relationship) <- list(levels(dataset$id), levels(dataset$id))
  attr(relationship, "inverse") <- TRUE
  randomModels <- list(
    ~vsm(fam(environment, 1L), dsm(site), dsm(block), ism(id), Gu=relationship),
    ~vsm(dsm(site), fam(environment, 1L), dsm(block), ism(id), Gu=relationship),
    ~vsm(dsm(site), dsm(block), fam(environment, 1L), ism(id), Gu=relationship))
  default <- mmes(response~site+environment+block, random=randomModels[[1L]], rcov=~units, data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=2L), update=list(blendWeight=0), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(default$engineDiagnostics$factorSchurActive)
  options(sommer.mme.factorSchurMultiFactor=TRUE)
  for(position in seq_along(randomModels)){
    fitModel <- function(solver){
      mmes(response~site+environment+block, random=randomModels[[position]], rcov=~units, data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=3L), update=list(blendWeight=0), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    latent <- fitModel("cholmod")
    reference <- fitModel("ldlt")
    expect_true(latent$engineDiagnostics$factorSchurActive)
    for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se")){
      expect_equal(latent[[component]], reference[[component]], tolerance=1e-6)
    }
  }
})
