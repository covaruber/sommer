.path_reporting_fit <- function(solver="cholmod", schur="auto", rr=FALSE,
                                blendWeight=NULL, verbose=TRUE, commonScale="off"){
  dataset <- expand.grid(id=factor(seq_len(5)), environment=factor(letters[1:3]),
                         replicate=seq_len(2))
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 3
  random <- if(rr) ~vsm(rrm(environment, 1), ism(id)) else ~vsm(ism(id))
  mmes(response~environment, random=random, data=dataset,
    control=mmes_control(solve=list(method=solver, blockSchur=schur),
      convergence=list(maxIters=1L),
      update=list(blendWeight=blendWeight, commonScale=commonScale),
      diagnostics=list(verbose=verbose, dateWarning=FALSE)))
}

test_that("the existing path line reports active and inactive engines", {
  output <- capture.output(fit <- .path_reporting_fit(schur="off", rr=TRUE))
  path <- grep("^MME path:", output, value=TRUE)
  expect_length(path, 1L)
  expect_match(path, "solver=cholmod", fixed=TRUE)
  expect_match(path, "latent RR=active", fixed=TRUE)
  expect_match(path, "Schur=inactive (disabled)", fixed=TRUE)
  expect_match(path, "chain=inactive; banded=inactive; factor-Schur=inactive", fixed=TRUE)
  expect_match(path, "residual=diagonal; grid=inactive; repeated-Kronecker=inactive", fixed=TRUE)
  expect_identical(fit$engineDiagnostics$residualEngine, "diagonal")
  quiet <- capture.output(reference <- .path_reporting_fit(schur="off", rr=TRUE,
                                                          verbose=FALSE))
  expect_length(grep("^MME path:", quiet), 0L)
  for(component in c("llik", "bu", "theta", "InfMat")){
    expect_equal(fit[[component]], reference[[component]])
  }
})

test_that("native RR and ordinary backends use the same status vocabulary", {
  native <- capture.output(.path_reporting_fit(schur="off", rr=TRUE, blendWeight=0))
  path <- grep("^MME path:", native, value=TRUE)
  expect_length(path, 1L)
  expect_match(path, "latent RR=active", fixed=TRUE)
  expect_match(path, "Schur=inactive (disabled)", fixed=TRUE)
  ordinary <- capture.output(.path_reporting_fit(solver="ldlt", blendWeight=0))
  path <- grep("^MME path:", ordinary, value=TRUE)
  expect_length(path, 1L)
  expect_match(path, "solver=ldlt; latent RR=inactive", fixed=TRUE)
  expect_match(path, "Schur=inactive (requires CHOLMOD)", fixed=TRUE)
})

test_that("common-scale profiling retains a single detailed path line", {
  output <- capture.output(.path_reporting_fit(schur="off", blendWeight=0,
                                              commonScale="force"))
  path <- grep("^MME path:", output, value=TRUE)
  expect_length(path, 1L)
  expect_match(path, "exact common-scale profile AI", fixed=TRUE)
  expect_match(path, "Schur=inactive (disabled)", fixed=TRUE)
})

test_that("active Schur and automatic decisions agree with fitted diagnostics", {
  for(mode in c("auto", "force")){
    output <- capture.output(fit <- .path_reporting_fit(schur=mode))
    path <- grep("^MME path:", output, value=TRUE)
    expect_length(path, 1L)
    expect_true(fit$engineDiagnostics$blockSchurActive)
    expect_match(path, "Schur=active", fixed=TRUE)
    expect_match(path, "latent RR=inactive", fixed=TRUE)
  }
})

test_that("grid residual status describes the actual residual engine", {
  dataset <- expand.grid(environment=factor(letters[1:2]),
                         row=factor(seq_len(3)), column=factor(seq_len(3)))
  dataset$id <- factor(rep(seq_len(5), length.out=nrow(dataset)))
  dataset$response <- sin(seq_len(nrow(dataset)))
  output <- capture.output(fit <- mmes(response~environment,
    random=~vsm(ism(id)),
    rcov=~dsumm(vsm(ar1m(row, rho=0.2), ar1m(column, rho=0.3), ism(units)),
                by=environment),
    data=dataset, control=mmes_control(convergence=list(maxIters=1L),
      update=list(blendWeight=0),
      diagnostics=list(verbose=TRUE, dateWarning=FALSE))))
  path <- grep("^MME path:", output, value=TRUE)
  expect_length(path, 1L)
  expect_identical(fit$engineDiagnostics$residualEngine, "grid")
  expect_match(path, "residual=grid; grid=active; repeated-Kronecker=inactive",
               fixed=TRUE)
})
