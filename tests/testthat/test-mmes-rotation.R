rotation_fixture <- function(){
  set.seed(712)
  ids <- paste0("g", seq_len(6))
  data <- expand.grid(
    id=ids,
    env=paste0("e", seq_len(2)),
    KEEP.OUT.ATTRS=FALSE
  )
  A <- outer(
    seq_along(ids),
    seq_along(ids),
    function(i, j) 0.6^abs(i-j)
  )
  dimnames(A) <- list(ids, ids)
  Gu <- solve(A)
  attr(Gu, "inverse") <- TRUE
  genetic <- drop(chol(A) %*% rnorm(length(ids)))
  data$y <- 3 + rep(c(-0.35, 0.35), each=length(ids)) +
    rep(genetic, 2) + rnorm(nrow(data), sd=0.25)
  list(data=data, Gu=Gu)
}

test_that("rotation preserves a balanced Henderson fit", {
  fixture <- rotation_fixture()
  common <- list(
    fixed=y~env, rcov=~units, data=fixture$data,
    nIters=12, verbose=FALSE, dateWarning=FALSE, computeCi=0,
    henderson=TRUE
  )
  ordinary <- do.call(
    mmes,
    c(common, list(random=~vsm(ism(id), Gu=fixture$Gu)))
  )
  rotated <- do.call(
    mmes,
    c(common, list(random=~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE)))
  )

  expect_equal(tail(ordinary$llik, 1), tail(rotated$llik, 1), tolerance=1e-7)
  expect_equal(unname(ordinary$theta), unname(rotated$theta), tolerance=1e-7)
  expect_equal(as.numeric(fitted(ordinary)), as.numeric(fitted(rotated)), tolerance=1e-7)
  expect_equal(ordinary$uList[[1]], rotated$uList[[1]], tolerance=1e-7)
  expect_equal(rotate_back_mmes(rotated), rotated$uList[[1]], tolerance=1e-12)
  expect_equal(residuals(rotated), rotated$y - fitted(rotated), tolerance=1e-12)

  ordinaryPrediction <- predict(ordinary, D="id")
  rotatedPrediction <- predict(rotated, D="id")
  expect_equal(
    ordinaryPrediction$pvals$predicted.value,
    rotatedPrediction$pvals$predicted.value,
    tolerance=1e-7
  )
  expect_equal(
    ordinaryPrediction$vcov,
    rotatedPrediction$vcov,
    tolerance=1e-7
  )

  expect_error(postPEV(rotated, mode=1), "require mode=2")
  ordinaryPev <- postPEV(ordinary, mode=2)
  rotatedPev <- postPEV(rotated, mode=2)
  expect_equal(
    ordinaryPev$uPevList[[1]],
    rotatedPev$uPevList[[1]],
    tolerance=1e-7
  )
})

test_that("direct rotation agrees with the Henderson rotation", {
  fixture <- rotation_fixture()
  common <- list(
    fixed=y~env,
    random=~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE),
    rcov=~units, data=fixture$data, nIters=15,
    verbose=FALSE, dateWarning=FALSE, computeCi=0
  )
  henderson <- do.call(mmes, c(common, list(henderson=TRUE)))
  direct <- do.call(mmes, c(common, list(henderson=FALSE)))

  expect_equal(unname(henderson$theta), unname(direct$theta), tolerance=2e-3)
  expect_equal(as.numeric(fitted(henderson)), as.numeric(fitted(direct)), tolerance=2e-3)
  expect_equal(henderson$uList[[1]], direct$uList[[1]], tolerance=2e-3)
  expect_identical(
    direct$rotation$formulation,
    "direct-observation-covariance"
  )
})

test_that("rotation respects covariance-coordinate blocks and shuffled rows", {
  fixture <- rotation_fixture()
  set.seed(91)
  data <- fixture$data[sample(nrow(fixture$data)), , drop=FALSE]
  common <- list(
    fixed=y~env, rcov=~units, data=data, nIters=12,
    verbose=FALSE, dateWarning=FALSE, computeCi=0,
    henderson=TRUE
  )
  ordinary <- do.call(
    mmes,
    c(common, list(random=~vsm(dsm(env), ism(id), Gu=fixture$Gu)))
  )
  rotated <- do.call(
    mmes,
    c(common, list(
      random=~vsm(dsm(env), ism(id), Gu=fixture$Gu, rotation=TRUE)
    ))
  )

  expect_equal(unname(ordinary$theta), unname(rotated$theta), tolerance=1e-6)
  expect_equal(as.numeric(fitted(ordinary)), as.numeric(fitted(rotated)), tolerance=1e-6)
  expect_equal(ordinary$uList[[1]], rotated$uList[[1]], tolerance=1e-6)
})

test_that("rotation rejects unsupported or incomplete models", {
  fixture <- rotation_fixture()
  expect_error(
    with(fixture$data, vsm(ism(id), rotation=TRUE)),
    "requires a supplied Gu"
  )

  badGu <- fixture$Gu
  badGu[1,1] <- -1
  attr(badGu, "inverse") <- TRUE
  expect_error(
    with(fixture$data, vsm(ism(id), Gu=badGu, rotation=TRUE)),
    "positive-definite"
  )

  incomplete <- fixture$data[-1, , drop=FALSE]
  expect_error(
    mmes(
      y~env,
      random=~vsm(dsm(env), ism(id), Gu=fixture$Gu, rotation=TRUE),
      rcov=~units, data=incomplete, nIters=1,
      verbose=FALSE, dateWarning=FALSE
    ),
    "every Gu level equally often|complete balanced layout"
  )

  expect_error(
    mmes(
      y~env,
      random=~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE),
      rcov=~units, data=fixture$data, nIters=1,
      computeCi=1, verbose=FALSE, dateWarning=FALSE
    ),
    "requires computeCi=0"
  )
})