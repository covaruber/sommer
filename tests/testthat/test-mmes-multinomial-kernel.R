multinomialKernel <- function(response, predictor, direction=matrix(0,
    nrow(predictor), ncol(predictor)), memory=1000){
  sommer:::sommer_multinomial_family_eval(response, predictor, direction, memory)
}

test_that("joint multinomial likelihood has category-block scores and curvature", {
  eta <- rbind(c(-1, .4, 1.8), c(.2, -.3, .7), c(3, -2, .5), c(-2, -1, -3))
  y <- 1:4
  value <- multinomialKernel(y, eta)
  reference <- exp(cbind(eta, 0))
  reference <- reference / rowSums(reference)
  expect_equal(value$probability, reference, tolerance=1e-12)
  expect_equal(as.numeric(value$logDensityRows), log(reference[cbind(1:4, y)]),
    tolerance=1e-12)
  expect_equal(value$logDensity, sum(log(reference[cbind(1:4, y)])), tolerance=1e-12)
  indicator <- matrix(0, 4, 3)
  indicator[cbind(1:3, 1:3)] <- 1
  expect_equal(value$score, indicator-reference[, 1:3], tolerance=1e-12)
  for(i in 1:4){
    p <- reference[i, 1:3]
    expected <- diag(p)-tcrossprod(p)
    expect_equal(value$curvature[,,i], expected, tolerance=1e-12)
    expect_gt(min(eigen(value$curvature[,,i], symmetric=TRUE)$values), 0)
    expect_true(all(value$curvature[,,i][lower.tri(expected)] < 0))
  }
  # A joint category block is not three independent Bernoulli curvatures.
  expect_gt(max(abs(value$curvature[1, 2, ])), .01)
  expect_equal(value$curvatureD1, array(0, c(3, 3, 4)))
  expect_equal(value$curvatureD2, array(0, c(3, 3, 4)))
})

test_that("multinomial block derivatives agree with independent finite differences", {
  set.seed(918)
  for(d in c(1L, 2L, 4L)){
    eta <- matrix(rnorm(5*d), 5, d)
    direction <- matrix(rnorm(5*d), 5, d)
    y <- rep(seq_len(d+1L), length.out=5)
    value <- multinomialKernel(y, eta, direction)
    h <- 1e-4
    plus <- multinomialKernel(y, eta+h*direction, direction)
    minus <- multinomialKernel(y, eta-h*direction, direction)
    expect_lte(max(abs((plus$curvature-minus$curvature)/(2*h)-value$curvatureD1)), 1e-6)
    expect_lte(max(abs((plus$curvatureD1-minus$curvatureD1)/(2*h)-
      value$curvatureD2)), 1e-6)
    expect_lte(max(abs((plus$curvature-2*value$curvature+minus$curvature)/h^2-
      value$curvatureD2)), 1e-6)
    expect_lte(max(abs((plus$logDensityRows-minus$logDensityRows)/(2*h)-
      rowSums(value$score*direction))), 1e-6)
    other <- matrix(rnorm(5*d),5,d)
    sumDirection <- multinomialKernel(y,eta,direction+other)
    otherDirection <- multinomialKernel(y,eta,other)
    mixed <- (sumDirection$curvatureD2-value$curvatureD2-otherDirection$curvatureD2)/2
    upMixed <- multinomialKernel(y,eta+h*other,direction)
    downMixed <- multinomialKernel(y,eta-h*other,direction)
    expect_lte(max(abs((upMixed$curvatureD1-downMixed$curvatureD1)/(2*h)-mixed)),1e-6)
    for(j in seq_len(d)){
      upper <- lower <- eta
      upper[,j] <- upper[,j]+h
      lower[,j] <- lower[,j]-h
      up <- multinomialKernel(y, upper)
      down <- multinomialKernel(y, lower)
      expect_lte(max(abs((up$logDensityRows-down$logDensityRows)/(2*h)-
        value$score[,j])), 1e-6)
      for(i in 1:5){
        expect_lte(max(abs((up$score[i,]-down$score[i,])/(2*h)+
          value$curvature[,j,i])), 1e-6)
      }
    }
  }
})

test_that("binary multinomial limit matches the native binomial likelihood through order four", {
  for(eta in c(-30, -3, 0, 2, 30)){
    for(y in 1:2){
      value <- multinomialKernel(y, matrix(eta, 1, 1), matrix(1, 1, 1))
      binary <- sommer:::sommer_glmm_family_eval(as.numeric(y==1), eta,
        list(list(family="binomial", link="logit", logParameter=0, free=FALSE)),
        1L, 1)
      expect_lte(abs(value$logDensity-binary$logDensity), 1e-12)
      expect_lte(abs(value$score[1]-binary$score[1]), 1e-12)
      expect_lte(abs(value$curvature[1]-binary$curvature[1]), 1e-12)
      expect_lte(abs(value$curvatureD1[1]-binary$curvatureD1[1]), 1e-12)
      expect_lte(abs(value$curvatureD2[1]-binary$curvatureD2[1]), 1e-12)
    }
  }
})

test_that("multinomial category permutations and absent outcomes preserve the joint experiment", {
  eta <- rbind(c(-1, .5, 2), c(1, 3, -.3))
  y <- c(1, 4)
  value <- multinomialKernel(y, eta)
  permutation <- c(4, 2, 1, 3)
  full <- cbind(eta, 0)[,permutation]
  rebased <- full[,1:3]-full[,4]
  changed <- multinomialKernel(match(y, permutation), rebased)
  expect_equal(changed$probability, value$probability[,permutation], tolerance=1e-12)
  expect_equal(changed$logDensity, value$logDensity, tolerance=1e-12)
  # Linear map from old nonreference logits to the rebased nonreference logits.
  basis <- rbind(diag(3), rep(0, 3))
  map <- basis[permutation[1:3],,drop=FALSE] -
    matrix(basis[permutation[4],], 3, 3, byrow=TRUE)
  expect_equal(changed$score %*% map, value$score, tolerance=1e-12)
  for(i in 1:2){
    expect_equal(t(map) %*% changed$curvature[,,i] %*% map,
      value$curvature[,,i], tolerance=1e-12)
  }
  direction <- rbind(c(.2,-.4,.8),c(-.5,.7,.1))
  original <- multinomialKernel(y,eta,direction)
  transformed <- multinomialKernel(match(y,permutation),rebased,direction %*% t(map))
  for(i in 1:2){
    expect_equal(t(map) %*% transformed$curvatureD1[,,i] %*% map,
      original$curvatureD1[,,i],tolerance=1e-12)
    expect_equal(t(map) %*% transformed$curvatureD2[,,i] %*% map,
      original$curvatureD2[,,i],tolerance=1e-12)
  }
  absent <- multinomialKernel(rep(1, 2), eta)
  expect_equal(absent$probability, value$probability)
  expect_equal(dim(absent$curvature), c(3L, 3L, 2L))
})

test_that("multinomial stable tails, memory admission and invalid inputs are explicit", {
  eta <- rbind(c(1000, -1000, 0), c(-1000, -900, -800), c(1000, 1000, 1000))
  value <- multinomialKernel(c(2, 4, 1), eta, matrix(1, 3, 3))
  expect_equal(as.numeric(value$logDensityRows), c(-2000, 0, -log(3)), tolerance=1e-12)
  expect_lte(max(abs(rowSums(value$probability)-1)), 1e-12)
  expect_true(all(is.finite(value$score)))
  expect_true(all(is.finite(value$curvatureD2)))
  near <- multinomialKernel(1, matrix(c(35,0),1), matrix(c(1,0),1))
  expect_gt(near$curvature[1,1,1], 0)
  saturated <- multinomialKernel(1, matrix(60,1,1),matrix(1,1,1))
  expect_equal(saturated$logDensity/(-log1p(exp(-60))),1,tolerance=1e-12)
  expect_equal(saturated$score[1]/plogis(-60),1,tolerance=1e-12)
  expect_equal(saturated$curvature[1]/plogis(-60),1,tolerance=1e-12)
  for(y in list(0, 5, 1.2, NA_real_, Inf)){
    expect_error(multinomialKernel(y, matrix(0,1,3)), "codes|finite")
  }
  expect_error(multinomialKernel(1, matrix(Inf,1,3)), "finite")
  expect_error(multinomialKernel(1, matrix(0,1,3), matrix(NA_real_,1,3)), "finite")
  expect_error(multinomialKernel(1, matrix(0,2,3)), "dimensions")
  expect_error(multinomialKernel(1, matrix(numeric(),0,3)), "dimensions")
  expect_error(multinomialKernel(1, matrix(0,1,3), matrix(0,1,2)), "dimensions")
  for(budget in c(0, -1, Inf, NA_real_, 1e-8)){
    expect_error(multinomialKernel(1, matrix(0,1,3), memory=budget), "denseMemoryMB")
  }
  set.seed(91)
  before <- .Random.seed
  invisible(multinomialKernel(1, matrix(0,1,3)))
  expect_identical(.Random.seed, before)
})

test_that("multinomial fixed-design likelihood optimization agrees with an independent softmax reference", {
  set.seed(691)
  n <- 180L
  x <- cbind(1, rnorm(n))
  eta <- x %*% rbind(c(.2,-.4),c(.5,-.3))
  weights <- exp(cbind(eta,0))
  probability <- weights/rowSums(weights)
  uniform <- runif(n)
  y <- 1L + (uniform > probability[,1]) +
    (uniform > probability[,1]+probability[,2])
  evaluate <- function(beta) multinomialKernel(y, x %*% matrix(beta,2,2))
  objective <- function(beta) -evaluate(beta)$logDensity
  gradient <- function(beta) -as.numeric(crossprod(x,evaluate(beta)$score))
  independent <- function(beta){
    eta <- cbind(x %*% matrix(beta,2,2),0)
    maximum <- apply(eta,1,max)
    logDenominator <- maximum+log(rowSums(exp(eta-maximum)))
    sum(logDenominator-eta[cbind(seq_len(n),y)])
  }
  native <- optim(rep(0,4),objective,gradient,method="BFGS",
    control=list(reltol=1e-12))
  reference <- optim(rep(0,4),independent,method="BFGS",
    control=list(reltol=1e-12))
  expect_equal(native$convergence,0L)
  expect_equal(reference$convergence,0L)
  expect_lte(abs(native$value-reference$value),1e-8)
  expect_lte(max(abs(native$par-reference$par)),1e-5)
  expect_lte(max(abs(gradient(native$par))),1e-4)
  value <- evaluate(native$par)
  hessian <- matrix(0,4,4)
  for(i in seq_len(n)){
    design <- kronecker(diag(2), x[i,,drop=FALSE])
    hessian <- hessian + t(design) %*% value$curvature[,,i] %*% design
  }
  expect_equal(hessian,optimHess(native$par,objective,gradient),tolerance=1e-6)
})
