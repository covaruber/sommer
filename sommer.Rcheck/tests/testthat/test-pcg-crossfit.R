.crossfit_data <- function(){
  dataset <- expand.grid(id=factor(seq_len(10L)), replicate=seq_len(4L))
  dataset$x <- rep(c(-1, 0, 1, 2), each=10L)
  dataset$y <- sin(seq_len(nrow(dataset)) * 0.7) +
    as.integer(dataset$id) / 5 + dataset$x / 3
  dataset
}

.crossfit_fit <- function(dataset, state=NULL, export=FALSE, iterations=2L){
  suppressMessages(mmes(y~x, random=~id, data=dataset, control=mmes_control(
    score=list(traceProbes=12L), logdet=list(lanczosSteps=40L),
    solve=list(method="pcg", recycleRank=4L, recycleMethod="krylov",
      recycleState=state, exportRecycle=export),
    convergence=list(maxIters=iterations), update=list(blendWeight=0),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
}

test_that("cross-fit recycle controls reject invalid or incompatible requests", {
  expect_null(mmes_control()$solve$recycleState)
  expect_false(mmes_control()$solve$exportRecycle)
  for(value in list(NA, 1, logical(), c(TRUE, FALSE))){
    expect_error(mmes_control(solve=list(exportRecycle=value)), "exportRecycle")
  }
  expect_error(mmes_control(solve=list(exportRecycle=TRUE)), "positive solve\\$recycleRank")
  expect_error(mmes_control(solve=list(recycleState=list())), "valid exported")
  pcg <- list(recycleRank=4L, recycleMethod="krylov", exportRecycle=TRUE)
  for(solver in c("auto", "ldlt", "cholmod")){
    expect_error(suppressMessages(mmes(y~x, random=~id, data=.crossfit_data(),
      control=mmes_control(solve=c(list(method=solver), pcg),
        output=list(setup=TRUE), diagnostics=list(verbose=FALSE)))),
      "solve\\$method='pcg'")
  }
  expect_error(suppressMessages(mmes(y~x, random=~id, data=.crossfit_data(),
    control=mmes_control(solve=c(list(method="pcg"), pcg),
      task=list(mode="solve"), diagnostics=list(verbose=FALSE)))),
    "task\\$mode='solve'")
})

test_that("explicit Krylov state survives serialization and changed observations", {
  dataset <- .crossfit_data()
  first <- .crossfit_fit(dataset, export=TRUE, iterations=1L)
  state <- unserialize(serialize(first$pcgRecycleState, NULL))
  expect_s3_class(state, "sommer_pcg_recycle")
  expect_identical(state$version, 1L)
  expect_equal(nrow(state$basis), nrow(first$bu))
  expect_gt(ncol(state$basis), 0L)
  expect_lte(ncol(state$basis), 4L)
  expect_true(all(is.finite(state$basis)))
  expect_null(first$pcgRecycleBasis)
  expect_equal(first$engineDiagnostics$pcg$importedRecycleRank, 0L)

  dataset <- dataset[-seq_len(10L), ]
  dataset$y <- dataset$y + cos(seq_len(nrow(dataset))) / 10
  original <- unserialize(serialize(state, NULL))
  reference <- .crossfit_fit(dataset)
  candidate <- .crossfit_fit(dataset, state, export=TRUE)
  expect_identical(state, original)
  expect_null(reference$pcgRecycleState)
  expect_equal(candidate$engineDiagnostics$pcg$importedRecycleRank, ncol(state$basis))
  expect_gt(candidate$engineDiagnostics$pcg$recyclePreparationOperatorApplications, 0)
  expect_gt(candidate$engineDiagnostics$pcg$balancedPreconditionerApplications, 0)
  expect_lte(candidate$pcgSolveMaxRelativeResidual, 1e-11)
  for(field in c("bu", "llik", "theta", "dLu", "InfMat", "theta_se")){
    expect_equal(candidate[[field]], reference[[field]], tolerance=1e-6, info=field)
  }
})

test_that("cross-fit states reject changed coefficient coordinates and malformed bases", {
  dataset <- .crossfit_data()
  state <- .crossfit_fit(dataset, export=TRUE)$pcgRecycleState
  changed <- dataset
  levels(changed$id) <- paste0("new", levels(changed$id))
  expect_error(.crossfit_fit(changed, state), "coefficient coordinates")
  changed <- dataset
  changed$id <- factor(changed$id, levels=rev(levels(changed$id)))
  expect_error(.crossfit_fit(changed, state), "coefficient coordinates")
  malformed <- state
  malformed$basis[1, 1] <- NA_real_
  expect_error(.crossfit_fit(dataset, malformed), "valid exported")
  malformed <- state
  malformed$basis <- malformed$basis[-1, , drop=FALSE]
  expect_error(.crossfit_fit(dataset, malformed), "coefficient coordinates")
  malformed <- state
  malformed$basis[,] <- 0
  expect_error(.crossfit_fit(dataset, malformed), "no independent directions")
  expect_error(mmes_control(solve=list(recycleState=state,
    recycleMethod="krylov", recycleRank=1L)), "rank exceeds")
})
