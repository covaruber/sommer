test_that("latent FA/RR Schur preserves free-shape AI trajectories and uncertainty", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 4
  precision <- diag(6L) + matrix(0.15, 6L, 6L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  for(model in c("fa", "rr")) for(rank in 1:2) for(inverseMode in 0:2){
    shape <- if(model == "fa") fam(dataset$environment, rank) else rrm(dataset$environment, rank)
    fitModel <- function(solver){
      mmes(response~environment,
        random=~vsm(shape, ism(id), Gu=precision) +
          vsm(dsm(environment), ism(replicate)) + vsm(ism(id)),
        rcov=~vsm(dsm(environment), ism(units)), data=dataset,
        solver=solver, nIters=4L, computeCi=inverseMode,
        getPEV=inverseMode > 0L, verbose=FALSE, dateWarning=FALSE)
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
  secondShape <- rrm(dataset$environment, 1L)
  secondShape$covFactor$free[] <- FALSE
  for(weighted in c(FALSE, TRUE)){
    weightMatrix <- if(weighted) diag(seq(0.7, 1.3, length.out=nrow(dataset))) else diag(nrow(dataset))
    fitModel <- function(solver){
      mmes(response~environment,
        random=~vsm(fam(environment, 2L), ism(id), Gu=precision) +
          vsm(secondShape, ism(replicate)),
        rcov=~vsm(dsm(environment), ism(units)), data=dataset, W=weightMatrix,
        solver=solver, nIters=4L, computeCi=2L, getPEV=FALSE,
        verbose=FALSE, dateWarning=FALSE)
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
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision),
      rcov=residual, data=dataset, solver=solver, nIters=2L, computeCi=0L,
      getPEV=FALSE, REML=reml, verbose=FALSE, dateWarning=FALSE)
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
  extraFactor <- mmes(response~environment,
    random=~vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu=precision),
    data=dataset, solver="cholmod", nIters=2L, computeCi=0L, getPEV=FALSE,
    verbose=FALSE, dateWarning=FALSE)
  expect_false(extraFactor$engineDiagnostics$factorSchurActive)
  expect_false(fitModel(shape, reml=FALSE)$engineDiagnostics$factorSchurActive)
  boundary <- fam(dataset$environment, 1L, specific=c(rep(1, 4L), 1e-10))
  boundary$covFactor$free[] <- FALSE
  fallback <- fitModel(boundary)
  expect_false(fallback$engineDiagnostics$factorSchurActive)
  reference <- fitModel(boundary, solver="ldlt")
  expect_equal(fallback$llik, reference$llik, tolerance=1e-6)
  expect_equal(fallback$bu, reference$bu, tolerance=1e-6)
  options(sommer.mme.denseMemoryMB=0.001)
  expect_false(fitModel(shape)$engineDiagnostics$factorSchurActive)
})