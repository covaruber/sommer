test_that("profiling reports phase totals without changing the fit", {
  workload <- generate_mme_workload(structure="ar2", n_levels=6L, n_ind=4L,
                                   seed=303L)
  fit <- function(profile="off", traceBatch="off"){
    mmes(workload$fixed, random=workload$random, data=workload$data,
      control=mmes_control(score=list(traceBatch=traceBatch),
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
