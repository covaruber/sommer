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
  options(sommer.mme.denseMemoryMB=0.001)
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
  options(sommer.mme.denseMemoryMB=1e-5)
  uncached <- fitModel()
  expect_equal(uncached$engineDiagnostics$CselectedDependencyEntries, 0)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu", "uPevList")){
    expect_equal(cached[[component]], uncached[[component]], tolerance=1e-10)
  }
  expect_equal(as.matrix(cached$Ci), as.matrix(uncached$Ci), tolerance=1e-10)
})

test_that("generic CHOLMOD trace fallbacks preserve fits under an exhausted cache budget", {
  previous <- options(sommer.mme.denseMemoryMB=1e-5)
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