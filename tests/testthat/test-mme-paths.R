test_that("MME path policies enforce dispatch and resource boundaries", {
  .native_run("test-mme-paths", "mme_paths.h", strict=TRUE)
})

test_that("LDLT selected-inverse dependency caching respects its budget and preserves output", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  data(DT_example)
  fitModel <- function(){
    mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(solve=list(method="ldlt"), convergence=list(maxIters=3L), output=list(inverse=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE), update=list(blendWeight=0)))
  }
  cached <- fitModel()
  expect_gt(cached$engineDiagnostics$CselectedDependencyEntries, 0)
  options(sommer.mme.denseMemoryMB=1e-5, sommer.mme.memoryPolicy="limit")
  on.exit(options(sommer.mme.memoryPolicy=NULL), add=TRUE)
  uncached <- fitModel()
  expect_equal(uncached$engineDiagnostics$CselectedDependencyEntries, 0)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "uPevList")){
    expect_equal(cached[[component]], uncached[[component]], tolerance=1e-10)
  }
  expect_equal(as.matrix(cached$Ci), as.matrix(uncached$Ci), tolerance=1e-10)
})

test_that("generic CHOLMOD trace fallbacks preserve fits under an exhausted cache budget", {
  previous <- options(sommer.mme.denseMemoryMB=1e-5, sommer.mme.memoryPolicy="limit")
  on.exit(options(previous), add=TRUE)
  data(DT_example)
  fitModel <- function(solver){
    mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=3L), output=list(inverse=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE), update=list(blendWeight=0)))
  }
  candidate <- fitModel("cholmod")
  reference <- fitModel("ldlt")
  expect_false(candidate$engineDiagnostics$blockSchurActive)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "uPevList")){
    expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
  }
})
test_that("memoryPolicy speed stops with the required memory and limit falls back", {
  previous <- options(sommer.mme.denseMemoryMB=1e-5, sommer.mme.memoryPolicy=NULL)
  on.exit(options(previous), add=TRUE)
  data(DT_example)
  fitModel <- function(policy){
    mmes(Yield~Env, random=~Name, data=DT_example, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(execution=list(memoryPolicy=policy))), list(solve=list(method="cholmod"), convergence=list(maxIters=2L), output=list(inverse=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE), update=list(blendWeight=0)), keep.null=TRUE)))
  }
  expect_equal(mmes_control()$execution$memoryPolicy, "speed")
  expect_error(mmes_control(execution=list(memoryPolicy="fast")),
    "memoryPolicy must be")
  expect_error(fitModel("speed"), "denseMemoryMB.*memoryPolicy")
  limited <- fitModel("limit")
  expect_equal(limited$control$execution$memoryPolicy, "limit")
  options(sommer.mme.denseMemoryMB=NULL)
  expect_equal(fitModel("speed")$control$execution$memoryPolicy, "speed")
})

test_that("sparse relationship blocks avoid dense Schur and banded engines are opt-in", {
  fitFor <- function(relationship, structure="us"){
    wl <- generate_mme_workload(structure=structure, n_levels=4L, n_ind=100L,
      relationship=relationship, seed=11L)
    mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data, control=mmes_control(convergence=list(maxIters=2L), solve=list(method="cholmod"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  expect_true(fitFor("dense")$engineDiagnostics$blockSchurActive)
  expect_false(fitFor("sparse")$engineDiagnostics$blockSchurActive)
  wl <- generate_mme_workload(structure="ar1", n_levels=50L, n_ind=8L, relationship="none")
  fit <- mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data, control=mmes_control(convergence=list(maxIters=2L), solve=list(method="cholmod"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(fit$engineDiagnostics$blockBandedActive)
})

test_that("blockSchur toggle disables or forces the block engines", {
  wl <- generate_mme_workload(structure="us", n_levels=4L, n_ind=100L, relationship="dense", seed=11L)
  fitFor <- function(mode) mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data,
    control=mmes_control(solve=list(method="cholmod", blockSchur=mode),
      convergence=list(maxIters=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(fitFor("auto")$engineDiagnostics$blockSchurActive)
  off <- fitFor("off")
  expect_false(off$engineDiagnostics$blockSchurActive)
  expect_equal(off$control$solve$blockSchur, "off")
  expect_equal(off$llik, fitFor("auto")$llik, tolerance=1e-7)
  sparse <- generate_mme_workload(structure="us", n_levels=4L, n_ind=100L, relationship="sparse", seed=11L)
  forced <- mmes(sparse$fixed, random=sparse$random, rcov=sparse$rcov, data=sparse$data,
    control=mmes_control(solve=list(method="cholmod", blockSchur="force"),
      convergence=list(maxIters=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(forced$engineDiagnostics$blockSchurActive)
})

test_that("fill guard rejects coupled low-fill blocks and retains tiny independent components", {
  set.seed(303L)
  size <- 128L
  precision <- diag(size)
  precision[row(precision) == col(precision) + 1L] <- -0.2
  upper <- which(upper.tri(precision))
  precision[sample(upper, 160L)] <- -0.1
  precision <- precision + t(precision)
  diag(precision) <- rowSums(abs(precision)) + 1
  dimnames(precision) <- rep(list(paste0("g", seq_len(size))), 2L)
  attr(precision, "inverse") <- TRUE
  dataset <- expand.grid(id=factor(rownames(precision)), environment=factor(1:6), replicate=1:2)
  dataset$response <- rnorm(nrow(dataset))
  fit <- function(mode, relationship){
    mmes(response~environment, random=~vsm(ar1m(environment), ism(id), Gu=relationship), data=dataset,
      control=mmes_control(solve=list(method="cholmod", blockSchur=mode),
        convergence=list(maxIters=2L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  automatic <- fit("auto", precision)
  forced <- fit("force", precision)
  expect_false(automatic$engineDiagnostics$blockSchurActive)
  expect_true(forced$engineDiagnostics$blockSchurActive)
  for(component in c("llik", "covPar", "bu", "InfMat", "dLu")){
    expect_equal(automatic[[component]], forced[[component]], tolerance=1e-7)
  }
  identity <- diag(size)
  dimnames(identity) <- dimnames(precision)
  attr(identity, "inverse") <- TRUE
  expect_true(fit("auto", identity)$engineDiagnostics$blockSchurActive)
})

test_that("independent groups still require enough sparse fill to justify dense Schur", {
  size <- 128L
  ids <- paste0("g", seq_len(size))
  dataset <- expand.grid(id=factor(ids, levels=ids),
    environment=factor(1:4), replicate=1:2)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / size
  sparse <- Matrix::bandSparse(size, k=c(0, 1), symmetric=TRUE,
    diagonals=list(rep(1.5, size), rep(-0.2, size - 1L)))
  dense <- diag(size) + matrix(0.08, size, size)
  dimnames(sparse) <- dimnames(dense) <- list(ids, ids)
  attr(sparse, "inverse") <- attr(dense, "inverse") <- TRUE
  fit <- function(relationship, mode, inverse=0L){
    mmes(response~environment,
      random=~vsm(dsm(environment), ism(id), Gu=relationship), data=dataset,
      control=mmes_control(solve=list(method="cholmod", blockSchur=mode),
        output=list(inverse=inverse), convergence=list(maxIters=2L),
        update=list(blendWeight=0),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  for(inverse in c(0L, 2L)){
    reference <- fit(sparse, "off", inverse)
    automatic <- fit(sparse, "auto", inverse)
    forced <- fit(sparse, "force", inverse)
    expect_false(automatic$engineDiagnostics$blockSchurActive)
    expect_true(forced$engineDiagnostics$blockSchurActive)
    expect_equal(forced$engineDiagnostics$blockSchurGroups, 4L)
    for(candidate in list(automatic, forced)){
      for(component in c("llik", "covPar", "bu", "InfMat", "dLu", "uPevList")){
        expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
      }
      if(inverse == 2L){
        expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
      }
    }
  }
  expect_true(fit(dense, "auto")$engineDiagnostics$blockSchurActive)
})
