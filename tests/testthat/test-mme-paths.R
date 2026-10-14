test_that("MME path policies enforce dispatch and resource boundaries", {
  root <- test_path("..", "..")
  source <- file.path(root, "tests", "native", "test-mme-paths.cpp")
  header <- file.path(root, "src", "mme_paths.h")
  skip_if(!file.exists(header), "Native path checks require the package source tree")
  compiler <- Sys.which("c++")
  skip_if(!nzchar(compiler), "C++ compiler unavailable")
  binary <- tempfile("sommer-mme-paths-")
  on.exit(unlink(binary), add=TRUE)
  output <- system2(compiler,
    c("-std=c++11", "-Wall", "-Wextra", "-pedantic", shQuote(source), "-o", shQuote(binary)),
    stdout=TRUE, stderr=TRUE)
  expect_null(attr(output, "status"), info=paste(output, collapse="\n"))
  if(!file.exists(binary)) return(invisible(NULL))
  output <- system2(binary, stdout=TRUE, stderr=TRUE)
  expect_null(attr(output, "status"), info=paste(output, collapse="\n"))
})

test_that("dense storage budget controls coupled Schur admission", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  dataset <- expand.grid(id=factor(seq_len(8L)), environment=factor(seq_len(3L)), replicate=seq_len(2L))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 8
  relationship <- diag(8L) + matrix(0.1, 8L, 8L)
  dimnames(relationship) <- list(levels(dataset$id), levels(dataset$id))
  attr(relationship, "inverse") <- TRUE
  fitModel <- function(){
    mmes(response~environment, random=~vsm(ar1m(environment), ism(id), Gu=relationship),
      data=dataset, solver="cholmod", nIters=2L, computeCi=0L, getPEV=FALSE, verbose=FALSE)
  }
  dense <- fitModel()
  expect_true(dense$engineDiagnostics$blockSchurActive)
  expect_equal(dense$engineDiagnostics$denseMemoryMB, 1000)
  expect_gt(dense$engineDiagnostics$blockSchurEstimatedDenseMB, 0.001)
  options(sommer.mme.denseMemoryMB=0.001, sommer.mme.memoryPolicy="limit")
  on.exit(options(sommer.mme.memoryPolicy=NULL), add=TRUE)
  sparse <- fitModel()
  expect_false(sparse$engineDiagnostics$blockSchurActive)
  expect_equal(sparse$llik, dense$llik, tolerance=1e-7)
  expect_equal(sparse$theta, dense$theta, tolerance=1e-7)
  expect_equal(sparse$bu, dense$bu, tolerance=1e-7)
  expect_equal(sparse$InfMat, dense$InfMat, tolerance=1e-7)
  for(invalid in list(0, -1, Inf, NA_real_, "large", c(1, 2))){
    options(sommer.mme.denseMemoryMB=invalid)
    expect_error(fitModel(), "sommer.mme.denseMemoryMB must be a positive finite number", fixed=TRUE)
  }
})

test_that("LDLT selected-inverse dependency caching respects its budget and preserves output", {
  previous <- options(sommer.mme.denseMemoryMB=NULL)
  on.exit(options(previous), add=TRUE)
  data(DT_example)
  fitModel <- function(){
    mmes(Yield~Env, random=~Name, data=DT_example, solver="ldlt",
      nIters=3L, computeCi=1L, verbose=FALSE, dateWarning=FALSE)
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
    mmes(Yield~Env, random=~Name, data=DT_example, solver=solver,
      nIters=3L, computeCi=1L, verbose=FALSE, dateWarning=FALSE)
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
    mmes(Yield~Env, random=~Name, data=DT_example, solver="cholmod", nIters=2L,
      computeCi=1L, verbose=FALSE, dateWarning=FALSE,
      engine=mmes_control(solver="cholmod", memoryPolicy=policy))
  }
  expect_equal(mmes_control()$memoryPolicy, "speed")
  expect_error(mmes_control(memoryPolicy="fast"), "memoryPolicy must be")
  expect_error(fitModel("speed"), "denseMemoryMB.*memoryPolicy")
  limited <- fitModel("limit")
  expect_equal(limited$engineDiagnostics$control$memoryPolicy, "limit")
  options(sommer.mme.denseMemoryMB=NULL)
  expect_equal(fitModel("speed")$engineDiagnostics$control$memoryPolicy, "speed")
})

test_that("sparse relationship blocks avoid dense Schur and banded engines are opt-in", {
  fitFor <- function(relationship, structure="us"){
    wl <- generate_mme_workload(structure=structure, n_levels=4L, n_ind=100L,
      relationship=relationship, seed=11L)
    mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data, nIters=2L,
      solver="cholmod", verbose=FALSE, dateWarning=FALSE)
  }
  expect_true(fitFor("dense")$engineDiagnostics$blockSchurActive)
  expect_false(fitFor("sparse")$engineDiagnostics$blockSchurActive)
  wl <- generate_mme_workload(structure="ar1", n_levels=50L, n_ind=8L, relationship="none")
  fit <- mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data, nIters=2L,
    solver="cholmod", verbose=FALSE, dateWarning=FALSE)
  expect_false(fit$engineDiagnostics$blockBandedActive)
})

test_that("blockSchur toggle disables or forces the block engines", {
  wl <- generate_mme_workload(structure="us", n_levels=4L, n_ind=100L, relationship="dense", seed=11L)
  fitFor <- function(mode) mmes(wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data, nIters=2L,
    engine=mmes_control(solver="cholmod", experimental=list(blockSchur=mode)),
    verbose=FALSE, dateWarning=FALSE)
  expect_true(fitFor("auto")$engineDiagnostics$blockSchurActive)
  off <- fitFor("off")
  expect_false(off$engineDiagnostics$blockSchurActive)
  expect_equal(off$engineDiagnostics$control$blockSchur, "off")
  expect_equal(off$llik, fitFor("auto")$llik, tolerance=1e-7)
  sparse <- generate_mme_workload(structure="us", n_levels=4L, n_ind=100L, relationship="sparse", seed=11L)
  forced <- mmes(sparse$fixed, random=sparse$random, rcov=sparse$rcov, data=sparse$data, nIters=2L,
    engine=mmes_control(solver="cholmod", experimental=list(blockSchur="force")),
    verbose=FALSE, dateWarning=FALSE)
  expect_true(forced$engineDiagnostics$blockSchurActive)
})

test_that("fill guard rejects coupled low-fill blocks and retains independent groups", {
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
    mmes(response~environment, random=~vsm(ar1m(environment), ism(id), Gu=relationship),
      data=dataset, nIters=2L, verbose=FALSE, dateWarning=FALSE,
      engine=mmes_control(solver="cholmod", experimental=list(blockSchur=mode)))
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
