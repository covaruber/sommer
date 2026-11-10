test_that("profiling reports phase totals without changing the fit", {
  workload <- generate_mme_workload(structure="ar2", n_levels=6L, n_ind=4L,
                                   seed=303L)
  fit <- function(profile="off", traceBatch="off"){
    sommer:::.mmes_dispatch(workload$fixed, random=workload$random, data=workload$data,
      control=sommer:::.mmes_control_internal(score=list(traceBatch=traceBatch),
        convergence=list(maxIters=2L),
        diagnostics=list(profile=profile, verbose=FALSE, dateWarning=FALSE)))
  }
  reference <- fit()
  profiled <- fit("force")
  batched <- fit("force", "force")
  expect_null(reference$engineDiagnostics$profile)
  for(candidate in list(profiled, batched)){
    expect_equal(candidate[c("llik", "theta", "bu", "InfMat", "theta_se")],
                 reference[c("llik", "theta", "bu", "InfMat", "theta_se")],
                 tolerance=1e-8)
  }
  profile <- profiled$engineDiagnostics$profile
  expect_named(profile$seconds, sommer:::.mme_phase_names)
  expect_true(all(is.finite(profile$seconds) & profile$seconds >= 0))
  phases <- setdiff(names(profile$seconds), c("setup", "final"))
  expect_equal(colSums(profile$secondsByIteration)[phases],
               profile$seconds[phases], tolerance=1e-9)
  expect_equal(nrow(profile$secondsByIteration), length(profiled$llik))
})

test_that("diagonal residuals do not build unused block design caches", {
  dataset <- data.frame(id=factor(rep(seq_len(20L), each=4L)))
  dataset$x <- 1 + sin(seq_len(nrow(dataset)))
  dataset$y <- cos(seq_len(nrow(dataset))) + as.integer(dataset$id) / 10
  fit <- sommer:::.mmes_dispatch(y~x, random=~id, data=dataset,
    control=sommer:::.mmes_control_internal(convergence=list(maxIters=2L),
      diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)))
  counters <- fit$engineDiagnostics$profile$counters
  expect_identical(fit$engineDiagnostics$residualEngine, "diagonal")
  expect_identical(counters$residualBlockDesignCacheBlocks, 0L)
  expect_identical(counters$residualBlockDesignCacheVisits, 0)
})

test_that("residual block design caches visit each design entry once and preserve row order", {
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(seq_len(3L)),
    replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$x <- 1 + sin(seq_len(nrow(dataset)))
  dataset$y <- cos(seq_len(nrow(dataset))) + as.integer(dataset$id) / 5
  for(incomplete in c(FALSE, TRUE)) for(weighted in if(incomplete) FALSE else c(FALSE, TRUE)){
    data <- if(incomplete) dataset[-c(1L, 7L), ] else dataset
    weights <- if(weighted) Matrix::Diagonal(nrow(data),
      x=1 + seq_len(nrow(data)) / (2*nrow(data))) else NULL
    fitModel <- function(order){
      args <- list(fixed=y~x, random=~id,
        rcov=~vsm(ar1m(environment, rho=0.3, fixed=TRUE), ism(record)),
        data=data[order, ],
        control=sommer:::.mmes_control_internal(convergence=list(maxIters=3L),
          diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)))
      if(weighted) args$W <- weights[order, order]
      do.call(sommer:::.mmes_dispatch, args)
    }
    reference <- fitModel(seq_len(nrow(data)))
    candidate <- fitModel(rev(seq_len(nrow(data))))
    counters <- reference$engineDiagnostics$profile$counters
    expect_equal(counters$residualBlockDesignCacheBlocks,
      nlevels(data$record))
    expect_equal(counters$residualBlockDesignCacheVisits,
      3*nrow(data))
    for(field in c("llik", "theta", "bu", "InfMat", "theta_se", "Ci")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8,
        info=paste(field, incomplete, weighted))
    }
  }
})
