test_that("benchmarking framework generates workloads and checks solver parity", {
  workload <- generate_mme_workload(
    name="test_ar1_workload",
    n_ind=20L,
    n_levels=4L,
    structure="ar1",
    seed=123L
  )
  expect_s3_class(workload, "mme_workload")
  expect_equal(workload$dimensions$n_ind, 20L)

  # Run benchmark with two engine configurations
  engines <- list(
    ldlt_base=mmes_control(solve=list(method="ldlt")),
    cholmod=mmes_control(solve=list(method="cholmod"))
  )
  benchmark <- run_mme_benchmark(
    workloads=workload,
    engines=engines,
    replications=1L,
    verifyParity=TRUE,
    tolerance=1e-4,
    nIters=3L,
    verbose=FALSE
  )
  expect_s3_class(benchmark, "sommer_mme_benchmark")
  expect_equal(nrow(benchmark$results), 2L)
  expect_true(all(benchmark$results$parity_ok))
})
