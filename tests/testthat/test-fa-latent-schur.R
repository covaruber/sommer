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

test_that("multi-factor latent admission prices rank, border and marginal connectivity", {
  fitModel <- function(e=6L, n=12L, m=2L, rank=1L, correlated=FALSE,
                       position="first", mode="auto", inverse=0L, specific=NULL){
    dataset <- expand.grid(id=factor(seq_len(n)), environment=factor(seq_len(e)),
      extra=factor(seq_len(m)), replicate=1:2)
    dataset$response <- sin(seq_len(nrow(dataset))/7) + as.integer(dataset$id)/n
    relationship <- diag(n) + matrix(0.1, n, n)
    dimnames(relationship) <- rep(list(levels(dataset$id)), 2L)
    attr(relationship, "inverse") <- TRUE
    shape <- if(is.null(specific)) fam(dataset$environment, rank) else
      fam(dataset$environment, rank, specific=specific)
    extraShape <- if(correlated) ar1m(dataset$extra, rho=0.4) else dsm(dataset$extra)
    random <- if(position == "first") ~vsm(shape, extraShape, ism(id), Gu=relationship) else
      ~vsm(extraShape, shape, ism(id), Gu=relationship)
    mmes(response~environment, random=random, data=dataset,
      control=mmes_control(solve=list(method="cholmod", blockSchur=if(mode == "generic") "off" else "auto",
        factorSchurMultiFactor=if(mode == "dense") "off" else if(mode == "forced") "force" else "auto"),
        output=list(inverse=inverse), convergence=list(maxIters=1L),
        update=list(blendWeight=0),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  for(position in c("first", "last")) for(inverse in c(0L, 2L)){
    automatic <- fitModel(position=position, inverse=inverse)
    expect_true(automatic$engineDiagnostics$factorSchurActive)
    planning <- automatic$engineDiagnostics$factorSchurPlanning
    expect_equal(planning$admission, "work-estimate")
    expect_lte(planning$latentWorkFlops, 0.8 * planning$marginalWorkFlops)
    reference <- fitModel(position=position, mode="generic", inverse=inverse)
    for(component in c("llik", "covPar", "bu", "InfMat", "dLu")){
      expect_equal(automatic[[component]], reference[[component]], tolerance=1e-7)
    }
    if(inverse == 2L){
      expect_equal(as.matrix(automatic$Ci), as.matrix(reference$Ci), tolerance=1e-7)
    }
  }
  diagonal <- fitModel(e=4L, m=3L)
  correlated <- fitModel(e=4L, m=3L, correlated=TRUE)
  expect_false(diagonal$engineDiagnostics$factorSchurActive)
  expect_equal(diagonal$engineDiagnostics$factorSchurPlanning$admission, "work-rejected")
  expect_true(correlated$engineDiagnostics$factorSchurActive)
  expect_equal(correlated$engineDiagnostics$factorSchurPlanning$admission, "work-estimate")
  expect_gt(correlated$engineDiagnostics$factorSchurPlanning$marginalWorkFlops,
    diagonal$engineDiagnostics$factorSchurPlanning$marginalWorkFlops)
  expect_false(fitModel(rank=3L)$engineDiagnostics$factorSchurActive)
  expect_false(fitModel(mode="dense")$engineDiagnostics$factorSchurActive)
  expect_true(fitModel(e=4L, m=3L, mode="forced")$engineDiagnostics$factorSchurActive)
  expect_equal(fitModel(e=8L)$engineDiagnostics$factorSchurPlanning$admission, "small-border")
  nearBoundary <- fitModel(specific=c(0.5, 1e-6, rep(0.7, 4L)))
  boundaryReference <- fitModel(specific=c(0.5, 1e-6, rep(0.7, 4L)), mode="generic")
  expect_false(nearBoundary$engineDiagnostics$factorSchurActive)
  expect_equal(nearBoundary$engineDiagnostics$factorSchurPlanning$admission,
    "work-numerical-rejected")
  expect_gt(nearBoundary$engineDiagnostics$factorSchurPlanning$numericalFallbacks, 0L)
  for(component in c("llik", "bu", "InfMat", "dLu")){
    expect_equal(nearBoundary[[component]], boundaryReference[[component]], tolerance=1e-7)
  }
})
