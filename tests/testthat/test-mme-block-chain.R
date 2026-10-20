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
        mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision) +
            vsm(dsm(environment), ism(replicate)) + vsm(ism(id)), rcov=~vsm(dsm(environment), ism(units)), data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=4L), output=list(inverse=inverseMode), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
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
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision), data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=2L), output=list(inverse=0L), task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  arShape <- ar1m(dataset$environment)
  expect_true(fitModel(arShape)$engineDiagnostics$blockChainActive)
  denseShape <- ownm(dataset$environment, K=diag(4L) + matrix(0.2, 4L, 4L))
  expect_false(fitModel(denseShape)$engineDiagnostics$blockChainActive)
  expect_false(fitModel(arShape, FALSE)$engineDiagnostics$blockChainActive)
  options(sommer.mme.denseMemoryMB=0.001, sommer.mme.memoryPolicy="limit")
  on.exit(options(sommer.mme.memoryPolicy=NULL), add=TRUE)
  expect_false(fitModel(arShape)$engineDiagnostics$blockChainActive)
})

test_that("bounded block-banded AR2 precision matches LDLT across inverse modes", {
  previous <- options(sommer.mme.blockBanded=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(12L)), range=factor(seq_len(10L)))
  dataset$response <- sin(seq_len(nrow(dataset)) / 11) + as.integer(dataset$id) / 6
  default <- mmes(response~1, random=~vsm(ar2m(range), ism(id)), rcov=~units, data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=3L), output=list(inverse=0L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(default$engineDiagnostics$blockBandedActive)
  options(sommer.mme.blockBanded=TRUE)
  for(inverseMode in 0:2){
    fitModel <- function(solver){
      mmes(response~1, random=~vsm(ar2m(range), ism(id)), rcov=~units, data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=3L), output=list(inverse=inverseMode), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("cholmod")
    expect_true(candidate$engineDiagnostics$blockBandedActive)
    expect_identical(candidate$engineDiagnostics$blockBandedWidth, 2L)
    for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
    }
    expect_equal(candidate$uPevList, reference$uPevList, tolerance=1e-7)
    expect_equal(as.matrix(candidate$C), as.matrix(reference$C), tolerance=1e-7)
    if(inverseMode == 2L){ expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7) }
  }
  options(sommer.mme.blockBanded=NA)
  expect_error(mmes(response~1, random=~vsm(ar2m(range), ism(id)), rcov=~units, data=dataset, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=1L), output=list(inverse=0L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "sommer.mme.blockBanded must be TRUE or FALSE")
})
