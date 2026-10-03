test_that("bordered random block chains match LDLT across correlation and inverse modes", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(letters[1:4]),
                         replicate=factor(seq_len(2L)))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 4
  dataset$record <- interaction(dataset$id, dataset$replicate)
  precision <- diag(6L) + matrix(0.15, 6L, 6L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  for(correlation in c(0, 0.3, -0.45)){
    for(inverseMode in 0:2){
      shape <- ar1m(dataset$environment, rho=correlation)
      fitModel <- function(solver){
        mmes(response~environment,
          random=~vsm(shape, ism(id), Gu=precision) +
            vsm(dsm(environment), ism(replicate)) + vsm(ism(id)),
          rcov=~vsm(dsm(environment), ism(units)), data=dataset,
          solver=solver, nIters=4L, computeCi=inverseMode,
          getPEV=inverseMode > 0L, verbose=FALSE, dateWarning=FALSE)
      }
      reference <- fitModel("ldlt")
      chain <- fitModel("cholmod")
      expect_true(chain$engineDiagnostics$blockChainActive)
      expect_equal(chain$engineDiagnostics$blockSchurGroups, 4L)
      for(component in c("llik", "theta", "bu", "InfMat")){
        expect_equal(chain[[component]], reference[[component]], tolerance=1e-7)
      }
      if(inverseMode == 2L){ expect_equal(as.matrix(chain$Ci), as.matrix(reference$Ci), tolerance=1e-7) }
      if(inverseMode > 0L){ expect_equal(chain$uPevList, reference$uPevList, tolerance=1e-7) }
    }
  }
})

test_that("chain solves support missing records and correlated residual trace fallbacks", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(7L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  dataset$response <- cos(seq_len(nrow(dataset))) + as.integer(dataset$id) / 5
  dataset$record <- interaction(dataset$id, dataset$replicate)
  precision <- diag(7L) + matrix(0.1, 7L, 7L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  for(missingRecords in c(FALSE, TRUE)) for(correlatedResidual in c(FALSE, TRUE)){
    activeData <- if(missingRecords) dataset[-c(2L, 11L, 40L), ] else dataset
    residual <- if(correlatedResidual) ~vsm(ar1m(environment), ism(record)) else
      ~vsm(dsm(environment), ism(units))
    fitModel <- function(solver){
      mmes(response~environment, random=~vsm(ar1m(environment), ism(id), Gu=precision),
        rcov=residual, data=activeData, solver=solver, nIters=4L,
        computeCi=2L, getPEV=FALSE, verbose=FALSE, dateWarning=FALSE)
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("cholmod")
    for(component in c("llik", "theta", "bu", "InfMat")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
    }
    expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
    if(!correlatedResidual || !missingRecords){ expect_true(candidate$engineDiagnostics$blockChainActive) }
  }
})

test_that("chain admission rejects nonpaths and respects memory and ML scope", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(letters[1:4]),
                         replicate=seq_len(2L))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  precision <- diag(6L) + matrix(0.1, 6L, 6L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  fitModel <- function(shape, reml=TRUE){
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision), data=dataset,
      solver="cholmod", nIters=2L, computeCi=0L, getPEV=FALSE, REML=reml,
      verbose=FALSE, dateWarning=FALSE)
  }
  arShape <- ar1m(dataset$environment)
  expect_true(fitModel(arShape)$engineDiagnostics$blockChainActive)
  denseShape <- ownm(dataset$environment, K=diag(4L) + matrix(0.2, 4L, 4L))
  expect_false(fitModel(denseShape)$engineDiagnostics$blockChainActive)
  expect_false(fitModel(arShape, FALSE)$engineDiagnostics$blockChainActive)
  options(sommer.mme.denseMemoryMB=0.001)
  expect_false(fitModel(arShape)$engineDiagnostics$blockChainActive)
})

test_that("reordered disconnected paths retain exact border-only cross-component inverses", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(7L)), environment=factor(letters[1:5]),
                         replicate=seq_len(2L))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 4
  precision <- diag(7L) + matrix(0.1, 7L, 7L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  covariance <- diag(5L)
  covariance[1L, 3L] <- covariance[3L, 1L] <- 0.25
  covariance[2L, 4L] <- covariance[4L, 2L] <- -0.25
  shape <- ownm(dataset$environment, K=covariance)
  fitModel <- function(solver){
    mmes(response~environment,
      random=~vsm(shape, ism(id), Gu=precision) + vsm(ism(id)),
      data=dataset, solver=solver, nIters=3L, computeCi=2L,
      getPEV=TRUE, verbose=FALSE, dateWarning=FALSE)
  }
  reference <- fitModel("ldlt")
  candidate <- fitModel("cholmod")
  expect_true(candidate$engineDiagnostics$blockChainActive)
  expect_equal(candidate$engineDiagnostics$blockSchurGroups, 5L)
  for(component in c("llik", "theta", "bu", "InfMat", "uPevList")){
    expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
  }
  expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
})