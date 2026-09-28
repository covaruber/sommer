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

expect_rotation_equivalent <- function(common, randomOrdinary, randomRotated,
                                       tolerance=1e-6){
  ordinary <- do.call(mmes, c(common, list(random=randomOrdinary)))
  rotated <- do.call(mmes, c(common, list(random=randomRotated)))
  expect_equal(tail(ordinary$llik, 1), tail(rotated$llik, 1), tolerance=tolerance)
  expect_equal(unname(ordinary$theta), unname(rotated$theta), tolerance=tolerance)
  expect_equal(as.numeric(fitted(ordinary)), as.numeric(fitted(rotated)),
               tolerance=tolerance)
  expect_equal(ordinary$uList[[1]], rotated$uList[[1]], tolerance=tolerance)
  invisible(rotated)
}

test_that("rotation supports multi-trait Kronecker residuals", {
  fixture <- rotation_fixture()
  set.seed(33)
  ids <- rownames(fixture$Gu)
  data <- expand.grid(id=ids, trait=c("t1", "t2"), KEEP.OUT.ATTRS=FALSE)
  genetic <- drop(chol(solve(as.matrix(fixture$Gu))) %*% rnorm(length(ids)))
  data$y <- 3 + (data$trait == "t2") * 0.5 +
    rep(genetic, 2) * rep(c(1, 0.6), each=length(ids)) +
    rnorm(nrow(data), sd=0.4)
  data <- data[sample(nrow(data)), , drop=FALSE]
  common <- list(fixed=y~trait, rcov=~vsm(usm(trait), ism(units)), data=data,
                 nIters=15, verbose=FALSE, dateWarning=FALSE, computeCi=0)
  for(henderson in c(TRUE, FALSE)){
    ordinary <- do.call(mmes, c(common, list(
      random=~vsm(usm(trait), ism(id), Gu=fixture$Gu), henderson=henderson
    )))
    rotated <- do.call(mmes, c(common, list(
      random=~vsm(usm(trait), ism(id), Gu=fixture$Gu, rotation=TRUE),
      henderson=henderson
    )))
    expect_equal(tail(ordinary$llik, 1), tail(rotated$llik, 1), tolerance=1e-6)
    expect_equal(unname(ordinary$theta), unname(rotated$theta), tolerance=1e-6)
    expect_equal(ordinary$uList[[1]], rotated$uList[[1]], tolerance=1e-6)
    if(henderson){
      expect_equal(as.numeric(fitted(ordinary)), as.numeric(fitted(rotated)),
                   tolerance=1e-6)
    }
  }
  expect_true(length(unique(rotated$rotation$residualBlock)) == 6L)
})

test_that("rotation supports heterogeneous residuals across rotation blocks", {
  fixture <- rotation_fixture()
  expect_rotation_equivalent(
    list(fixed=y~env, rcov=~vsm(dsm(env), ism(units)), data=fixture$data,
         nIters=15, verbose=FALSE, dateWarning=FALSE, computeCi=0,
         henderson=TRUE),
    ~vsm(ism(id), Gu=fixture$Gu),
    ~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE)
  )
})

test_that("rotation supports weights that are constant within rotation blocks", {
  fixture <- rotation_fixture()
  W <- diag(ifelse(fixture$data$env == "e1", 1, 2))
  expect_rotation_equivalent(
    list(fixed=y~env, rcov=~units, data=fixture$data, W=W,
         nIters=15, verbose=FALSE, dateWarning=FALSE, computeCi=0,
         henderson=TRUE),
    ~vsm(ism(id), Gu=fixture$Gu),
    ~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE)
  )
})

test_that("rotation rejects residual structures that are not invariant", {
  fixture <- rotation_fixture()
  data <- fixture$data
  data$grp <- factor(ifelse(data$id %in% c("g1", "g2", "g3"), "a", "b"))
  fit <- function(rcov, W=NULL){
    args <- list(y~env,
                 random=~vsm(ism(id), Gu=fixture$Gu, rotation=TRUE),
                 rcov=rcov, data=data, nIters=1,
                 verbose=FALSE, dateWarning=FALSE)
    if(!is.null(W)) args$W <- W
    do.call(mmes, args)
  }
  expect_error(fit(~vsm(ar1m(id), ism(units))), "different relationship levels")
  expect_error(fit(~vsm(dsm(grp), ism(units))), "same residual covariance coordinate")
  expect_error(fit(~units, W=diag(seq_len(nrow(data)))), "constant weight")
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