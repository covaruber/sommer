test_that("benchmark records frozen-point AI parity rather than just final estimates", {
  workload <- generate_mme_workload(name="ai-parity", n_ind=5L, n_levels=3L,
    structure="ar1", seed=771L)
  for(solver in c("ldlt", "cholmod")){
    result <- sommer:::.mme_benchmark_ai_assembly(workload,
      sommer:::.mmes_control_internal(solve=list(method=solver)), replications=1L)
    expect_equal(nrow(result), 2L)
    expect_identical(result$assembly, c("sensitivity", "augmented-working"))
    expect_true(all(result$solver == solver))
    expect_true(all(result$ai_max_rel_diff <= 1e-8))
    expect_true(all(result$score_max_rel_diff <= 1e-8))
    expect_true(all(result$coefficient_max_rel_diff <= 1e-8))
    expect_true(all(is.finite(result$elapsed_sec) & result$elapsed_sec >= 0))
  }
})
