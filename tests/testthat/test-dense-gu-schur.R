test_that("fixed factor and relationship components split tensor priors without diagonal-only rules", {
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(1:2), trait=factor(1:2))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  precision <- matrix(0, 6L, 6L)
  precision[1:3, 1:3] <- diag(3L) + matrix(0.1, 3L, 3L)
  precision[4:6, 4:6] <- diag(3L) + matrix(0.1, 3L, 3L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  fitModel <- function(solver){
    mmes(response~environment+trait, random=~vsm(ar1m(environment), usm(trait), ism(id), Gu=precision), data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=3L), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  reference <- fitModel("ldlt")
  candidate <- fitModel("cholmod")
  expect_true(candidate$engineDiagnostics$blockSchurActive)
  expect_equal(candidate$engineDiagnostics$blockSchurGroups, 2L)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se", "uPevList")){
    expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
  }
  expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
})

test_that("multiple dense-Gu effects use exact Schur blocks", {
  set.seed(17)
  nInd <- 18
  ids <- paste0("id", seq_len(nInd))
  markers <- matrix(sample(c(-1, 0, 1), nInd * 40, replace=TRUE), nrow=nInd)
  G1 <- tcrossprod(scale(markers)) / ncol(markers) + diag(0.2, nInd)
  G2 <- 0.7 * G1 + 0.3 * diag(nInd)
  Ginv1 <- solve(G1)
  Ginv2 <- solve(G2)
  dimnames(Ginv1) <- dimnames(Ginv2) <- list(ids, ids)
  attr(Ginv1, "inverse") <- TRUE
  attr(Ginv2, "inverse") <- TRUE

  dat <- expand.grid(id=ids, env=paste0("E", 1:3), KEEP.OUT.ATTRS=FALSE)
  dat$id <- factor(dat$id, levels=ids)
  dat$env <- factor(dat$env)
  dat$y <- rnorm(nrow(dat)) + as.numeric(G1 %*% rnorm(nInd))[as.integer(dat$id)]
  random <- ~vsm(dsm(env), ism(id), Gu=Ginv1) +
    vsm(dsm(env), ism(id), Gu=Ginv2)

  cholmodOutput <- capture.output(cholmodFit <- mmes(y~env, random=random, rcov=~units, data=dat, control=mmes_control(solve=list(method="cholmod"), convergence=list(maxIters=8), diagnostics=list(verbose=TRUE, dateWarning=FALSE))))
  ldltFit <- suppressMessages(mmes(y~env, random=random, rcov=~units, data=dat, control=mmes_control(solve=list(method="ldlt"), convergence=list(maxIters=8), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))

  expect_true(any(grepl("Dense block-Schur engine active \\(3 random-effect blocks, 3 border effects\\)", cholmodOutput)))
  expect_false(any(grepl("Materialized cross-group inverse block", cholmodOutput)))
  expect_equal(unlist(cholmodFit$theta), unlist(ldltFit$theta), tolerance=1e-8)
  expect_equal(cholmodFit$u, ldltFit$u, tolerance=1e-8)
  expect_true(cholmodFit$engineDiagnostics$blockSchurActive)
  expect_equal(cholmodFit$engineDiagnostics$blockSchurGroups, 3L)

  mlSchur <- suppressMessages(mmes(y~env, random=random, rcov=~units, data=dat, control=mmes_control(solve=list(method="cholmod"), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=8), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
  mlLdlt <- suppressMessages(mmes(y~env, random=random, rcov=~units, data=dat, control=mmes_control(solve=list(method="ldlt"), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=8), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
  expect_true(mlSchur$engineDiagnostics$dFirstMLActive)
  expect_equal(unlist(mlSchur$theta), unlist(mlLdlt$theta), tolerance=1e-8)
  expect_equal(as.numeric(mlSchur$llik), as.numeric(mlLdlt$llik), tolerance=1e-8)
  expect_equal(mlSchur$bu, mlLdlt$bu, tolerance=1e-8)
})

test_that("partial diagonal tensor products expose exact Schur groups in any factor order", {
  dataset <- expand.grid(id=factor(seq_len(4L)), site=factor(1:2),
    time=factor(1:2), zone=factor(1:2), trait=factor(1:2), batch=factor(1:2))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  precision <- diag(4L) + matrix(0.1, 4L, 4L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  randomModels <- list(
    ~vsm(dsm(site), ar1m(time), usm(trait), ism(id), Gu=precision),
    ~vsm(ar1m(time), dsm(site), usm(trait), dsm(zone), ism(id), Gu=precision),
    ~vsm(dsm(zone), usm(trait), dsm(batch), ar1m(time), dsm(site), ism(id), Gu=precision))
  cases <- list(
    list(model = 1L, reml = TRUE, inverseMode = 0L),
    list(model = 2L, reml = TRUE, inverseMode = 2L),
    list(model = 3L, reml = FALSE, inverseMode = 2L)
  )
  for(case in cases){
    model <- case$model
    reml <- case$reml
    inverseMode <- case$inverseMode
    fitModel <- function(solver){
      mmes(response~site+time+zone+trait+batch, random=randomModels[[model]], data=dataset, control=mmes_control(solve=list(method=solver), task=list(objective=if(reml) "reml" else "ml"), convergence=list(maxIters=2L), output=list(inverse=inverseMode), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("cholmod")
    if(!reml && inverseMode != 1L){
      expect_true(candidate$engineDiagnostics$dFirstMLActive)
    }else{
      expect_true(candidate$engineDiagnostics$blockSchurActive)
      expect_equal(candidate$engineDiagnostics$blockSchurGroups, as.integer(2^model))
    }
    for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
    }
    expect_equal(as.matrix(candidate$C), as.matrix(reference$C), tolerance=1e-7)
    if(inverseMode > 0L){ expect_equal(candidate$uPevList, reference$uPevList, tolerance=1e-7) }
    if(inverseMode == 2L){ expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7) }
  }
})

test_that("correlated ragged residuals merge tensor prior components exactly", {
  dataset <- expand.grid(id=factor(seq_len(4L)), site=factor(1:2),
    time=factor(1:2), zone=factor(1:2), trait=factor(1:2))
  dataset$record <- interaction(dataset$id, dataset$time, dataset$zone, dataset$trait)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  dataset <- dataset[-1L, ]
  precision <- diag(4L) + matrix(0.1, 4L, 4L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  fitModel <- function(solver){
    mmes(response~site+time+zone+trait, random=~vsm(dsm(site), ar1m(time), dsm(zone), usm(trait), ism(id), Gu=precision), rcov=~vsm(ar1m(site, rho=0.35, fixed=TRUE), ism(record)), data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=3L), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  reference <- fitModel("ldlt")
  candidate <- fitModel("cholmod")
  expect_true(candidate$engineDiagnostics$blockSchurActive)
  expect_equal(candidate$engineDiagnostics$blockSchurGroups, 2L)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se", "uPevList")){
    expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
  }
  expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
})
