multinomialModeEval <- function(context,predictor,start,maxit=100L,tolerance=1e-10){
  sommer:::sommer_multinomial_mode_eval(context,predictor,start,maxit,tolerance)
}

multinomialModeReference <- function(fixture){
  n <- length(fixture$response)
  d <- ncol(fixture$predictor)
  Q <- as.matrix(fixture$precision)
  Z <- as.matrix(fixture$design)
  logDensity <- function(u){
    eta <- cbind(fixture$predictor + matrix(Z %*% u,n,d),0)
    maximum <- apply(eta,1,max)
    sum(eta[cbind(seq_len(n),fixture$response)] -
      maximum-log(rowSums(exp(eta-maximum)))) -
      drop(crossprod(u,Q %*% u))/2
  }
  list(logDensity=logDensity,
    mode=optim(rep(0,ncol(Z)),function(u)-logDensity(u),method="BFGS",
      control=list(reltol=1e-13,maxit=1000)))
}

test_that("multinomial block mode and Laplace likelihood match independent dense calculations", {
  f <- multinomialModeFixture()
  context <- multinomialModeContext(f)
  value <- multinomialModeEval(context,f$predictor,rep(0,3))
  reference <- multinomialModeReference(f)
  expect_true(value$valid)
  expect_equal(reference$mode$convergence,0L)
  expect_equal(as.numeric(value$mode),reference$mode$par,tolerance=1e-6)
  expect_lte(value$modeScore,1e-9)
  eta <- f$predictor + matrix(f$design %*% value$mode,5,2)
  probability <- exp(cbind(eta,0))
  probability <- probability/rowSums(probability)
  expect_equal(value$linearPredictors,eta,tolerance=1e-12)
  expect_equal(value$probability,probability,tolerance=1e-12)
  expected <- as.matrix(f$precision)
  diagonalOnly <- expected
  for(row in 1:5){
    p <- probability[row,1:2]
    block <- diag(p)-tcrossprod(p)
    z <- as.matrix(f$design[c(row,row+5),])
    expected <- expected + crossprod(z,block %*% z)
    diagonalOnly <- diagonalOnly + crossprod(z,diag(diag(block)) %*% z)
  }
  expect_equal(as.matrix(value$modeHessian),expected,tolerance=1e-12)
  expect_gt(max(abs(expected-diagonalOnly)),.01)
  numerical <- optimHess(as.numeric(value$mode),function(u)-reference$logDensity(u))
  expect_equal(numerical,expected,tolerance=2e-6)
  logDet <- as.numeric(determinant(expected,logarithm=TRUE)$modulus)
  expectedLogLik <- reference$logDensity(as.numeric(value$mode)) +
    as.numeric(determinant(f$precision,logarithm=TRUE)$modulus)/2 - logDet/2
  expect_equal(value$logDetHessian,logDet,tolerance=1e-12)
  expect_equal(value$logLik,expectedLogLik,tolerance=1e-12)
  expect_equal(value$work$symbolicAnalyses,1L)
  expect_gte(value$work$factorizations,2L)
})

test_that("multinomial modes preserve binary equivalence and nonzero offsets", {
  X <- Matrix::Matrix(cbind(1,c(-1,0,1,.5)),sparse=TRUE)
  Z <- Matrix::Matrix(matrix(c(1,1,0,0,0,1,1,-1),4),sparse=TRUE)
  Q <- Matrix::Matrix(matrix(c(2,.3,.3,1.5),2),sparse=TRUE)
  offset <- c(-.1,0,.2,-.2)
  beta <- c(.2,.4)
  y <- c(0,1,1,0)
  binary <- sommer:::sommer_laplace_glmm_eval(X,Z,Q,y,offset,beta,c(0,0),
    "binomial",log(2*1.5-.3^2),100L,1e-10)
  f <- list(response=ifelse(y==1,1,2),design=Z,precision=Q,
    predictor=matrix(as.numeric(X %*% beta)+offset,ncol=1))
  value <- multinomialModeEval(multinomialModeContext(f),f$predictor,c(0,0))
  expect_true(value$valid)
  expect_equal(value$logLik,binary$logLik,tolerance=1e-12)
  expect_equal(value$mode,binary$mode,tolerance=1e-10)
  expect_equal(value$modeHessian,binary$modeHessian,tolerance=1e-10)
  expect_equal(value$probability[,1],as.numeric(binary$fitted),tolerance=1e-12)
})

test_that("multinomial mode contexts reuse symbolic structure through saturation and evaluation order", {
  f <- multinomialModeFixture()
  context <- multinomialModeContext(f)
  base <- multinomialModeEval(context,f$predictor,rep(0,3))
  extreme <- matrix(rep(c(1000,-1000),each=5),5,2)
  saturated <- multinomialModeEval(context,extreme,rep(0,3))
  again <- multinomialModeEval(context,f$predictor,c(2,-3,1))
  reverse <- multinomialModeContext(f)
  invisible(multinomialModeEval(reverse,extreme,rep(0,3)))
  reversed <- multinomialModeEval(reverse,f$predictor,rep(0,3))
  for(value in list(base,saturated,again,reversed)){
    expect_true(value$valid)
    expect_equal(value$work$symbolicAnalyses,1L)
  }
  expect_equal(again$logLik,base$logLik,tolerance=1e-12)
  expect_equal(reversed$logLik,base$logLik,tolerance=1e-12)
  expect_equal(again$mode,base$mode,tolerance=1e-9)
  expect_equal(again$work$evaluations,3L)
  failed <- multinomialModeEval(context,f$predictor,rep(50,3),maxit=1L)
  expect_false(failed$valid)
  expect_match(failed$reason,"iteration limit")
  expect_false("logLik" %in% names(failed))
  overflow <- multinomialModeEval(context,f$predictor,rep(.Machine$double.xmax,3))
  expect_false(overflow$valid)
  expect_match(overflow$reason,"non-finite")
  expect_false("logLik" %in% names(overflow))
  expect_true(multinomialModeEval(context,f$predictor,rep(0,3))$valid)
  zero <- f
  zero$design <- Matrix::Matrix(0,nrow=10,ncol=3,sparse=TRUE)
  independent <- multinomialModeEval(multinomialModeContext(zero),f$predictor,rep(0,3))
  kernel <- sommer:::sommer_multinomial_family_eval(f$response,f$predictor,matrix(0,5,2))
  expect_true(independent$valid)
  expect_equal(independent$iterations,0L)
  expect_equal(independent$logLik,kernel$logDensity,tolerance=1e-12)
  expect_equal(as.matrix(independent$modeHessian),as.matrix(f$precision),tolerance=1e-12)
})

test_that("multinomial category reference changes preserve modes and integrated likelihood", {
  f <- multinomialModeFixture()
  original <- multinomialModeEval(multinomialModeContext(f),f$predictor,rep(0,3))
  # Swap reference (old category 3) with old category 1.
  map <- rbind(c(-1,0),c(-1,1))
  changed <- f
  changed$response <- match(f$response,c(3,2,1))
  changed$predictor <- f$predictor %*% t(map)
  changed$design <- Matrix::kronecker(Matrix::Matrix(map,sparse=TRUE),
    Matrix::Diagonal(5)) %*% f$design
  value <- multinomialModeEval(multinomialModeContext(changed),changed$predictor,rep(0,3))
  expect_true(value$valid)
  expect_equal(value$logLik,original$logLik,tolerance=1e-12)
  expect_equal(value$mode,original$mode,tolerance=1e-9)
  expect_equal(value$modeHessian,original$modeHessian,tolerance=1e-10)
  expect_equal(value$probability,original$probability[,c(3,2,1)],tolerance=1e-10)
})

test_that("multinomial mode guards reject bad contexts, precision, dimensions and budgets", {
  f <- multinomialModeFixture()
  context <- multinomialModeContext(f)
  for(pointer in list(NULL,1,new.env())){
    expect_error(multinomialModeEval(pointer,f$predictor,rep(0,3)),"Invalid multinomial")
  }
  serialized <- unserialize(serialize(context,NULL))
  expect_error(multinomialModeEval(serialized,f$predictor,rep(0,3)),"fresh context")
  expect_error(multinomialModeEval(context,matrix(0,4,2),rep(0,3)),"dimensions")
  expect_error(multinomialModeEval(context,f$predictor,rep(0,2)),"dimensions")
  expect_error(multinomialModeEval(context,f$predictor,rep(Inf,3)),"finite")
  expect_error(multinomialModeEval(context,f$predictor,rep(0,3),maxit=0),"positive")
  expect_error(multinomialModeEval(context,f$predictor,rep(0,3),tolerance=NA_real_),"positive")
  for(budget in c(0,NA_real_,Inf,1e-8)){
    expect_error(multinomialModeContext(f,budget),"denseMemoryMB")
  }
  baseBytes <- 8*(16*5*2^2 + 24*5*2 + 16*3 +
    8*sum(as.matrix(f$design) != 0) + 8*sum(as.matrix(f$precision) != 0))
  expect_error(multinomialModeContext(f,(baseBytes+1)/1048576),"structural Hessian")
  expect_error(sommer:::sommer_multinomial_mode_create(f$response,
    f$design,f$precision,0L),"dimensions")
  expect_error(sommer:::sommer_multinomial_mode_create(f$response,
    f$design,f$precision,.Machine$integer.max),"sparse index limits")
  expect_error(sommer:::sommer_multinomial_mode_create(f$response,
    f$design[-1,],f$precision,2L),"dimensions")
  bad <- f
  bad$precision[1,2] <- .5
  expect_error(multinomialModeContext(bad),"symmetric")
  bad$precision <- Matrix::Diagonal(3,c(1,-1,1))
  expect_warning(
    expect_error(multinomialModeContext(bad),"positive definite"),
    "CHOLMOD warning.*not positive definite")
  bad <- f
  bad$response[1] <- 1.5
  expect_error(multinomialModeContext(bad),"codes")
  bad <- f
  bad$design[1,1] <- NA_real_
  expect_error(multinomialModeContext(bad),"finite")
  set.seed(917)
  before <- .Random.seed
  invisible(multinomialModeEval(context,f$predictor,rep(0,3)))
  expect_identical(.Random.seed,before)
})

test_that("sparse multinomial modes agree across category and random-effect dimensions", {
  set.seed(913)
  for(d in c(1L,3L,5L)){
    n <- 9L
    q <- d+1L
    Z <- matrix(rnorm(n*d*q),n*d,q)
    Z[abs(Z)<.7] <- 0
    basis <- matrix(rnorm(q*q),q,q)
    f <- list(response=sample.int(d+1L,n,replace=TRUE),
      predictor=matrix(rnorm(n*d,sd=2),n,d),
      design=Matrix::Matrix(Z,sparse=TRUE),
      precision=Matrix::Matrix(crossprod(basis)+diag(.7,q),sparse=TRUE))
    value <- multinomialModeEval(multinomialModeContext(f),f$predictor,rep(0,q))
    reference <- multinomialModeReference(f)
    numerical <- optimHess(as.numeric(value$mode),function(u)-reference$logDensity(u))
    expect_true(value$valid)
    expect_equal(reference$mode$convergence,0L)
    expect_equal(as.numeric(value$mode),reference$mode$par,tolerance=2e-6)
    expect_equal(as.matrix(value$modeHessian),numerical,tolerance=1e-5)
    expect_equal(rowSums(value$probability),rep(1,n),tolerance=1e-14)
    expect_true(is.finite(value$logLik))
  }
})
