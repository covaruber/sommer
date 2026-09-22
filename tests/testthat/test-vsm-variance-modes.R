test_that("correlation variance modes create consistent covariance descriptors", {
  levels3 <- factor(rep(letters[1:3], each=2))
  levels4 <- factor(rep(1:4, each=2))
  levels5 <- factor(rep(1:5, each=2))

  structures <- list(
    csm(levels3, variance="heterogeneous"),
    ar1m(levels4, variance="heterogeneous"),
    ar2m(levels4, variance="heterogeneous"),
    ar3m(levels5, variance="heterogeneous")
  )

  for(structure in structures){
    descriptor <- structure$covFactor
    expect_true(length(descriptor$par) == length(descriptor$free))
    expect_true(length(descriptor$par) == length(descriptor$par_names))
  }

  for(structure in structures[2:4]){
    descriptor <- structure$covFactor
    covariance <- descriptor$evaluator$fun(descriptor$par)
    expect_true(isSymmetric(covariance))
    expect_true(all(is.finite(covariance)))
    expect_gt(min(eigen(covariance, symmetric=TRUE, only.values=TRUE)$values), 0)
  }
})

test_that("heterogeneous variance mode validates its fixed parameters", {
  levels4 <- factor(rep(1:4, each=2))

  expect_error(
    ar1m(levels4, variance="heterogeneous", fixed=FALSE),
    "length order + q - 1"
  )
  expect_error(
    csm(levels4, variance="heterogeneous", fixed=TRUE),
    "length q"
  )
})