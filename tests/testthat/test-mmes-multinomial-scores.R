multinomialScoreFixture <- function(){
  f <- multinomialModeFixture()
  f$fixedDesign <- Matrix::Matrix(cbind(
    rep(c(1,0),each=5),rep(c(0,1),each=5),
    c(-1,0,.4,.7,1,.3,-.6,.1,.8,-.4)),sparse=TRUE)
  f
}

multinomialScorePrecision <- function(parameters=c(.2,-.1,.3,.25)){
  L <- diag(exp(parameters[1:3]))
  L[2,1] <- parameters[4]
  L[3,1] <- -.15
  derivatives <- lapply(1:4,function(k){
    derivative <- matrix(0,3,3)
    if(k<=3) derivative[k,k] <- L[k,k] else derivative[2,1] <- 1
    Matrix::Matrix(tcrossprod(derivative,L)+tcrossprod(L,derivative),sparse=TRUE)
  })
  list(value=Matrix::Matrix(tcrossprod(L),sparse=TRUE),derivatives=derivatives)
}

multinomialScoreEval <- function(context,f,precision=multinomialScorePrecision(),
    gradient=TRUE,start=rep(0,ncol(f$design)),maxit=100L){
  sommer:::sommer_multinomial_laplace_eval(context,f$predictor,precision$value,
    f$fixedDesign,precision$derivatives,start,maxit,1e-11,gradient)
}

multinomialScoreDense <- function(f,value,precision){
  n <- length(f$response)
  d <- ncol(f$predictor)
  Z <- as.matrix(f$design)
  mode <- as.numeric(value$mode)
  p <- value$probability
  H <- as.matrix(value$modeHessian)
  W <- matrix(0,n*d,n*d)
  score <- -p[,seq_len(d),drop=FALSE]
  for(row in seq_len(n)){
    if(f$response[row]<=d) score[row,f$response[row]] <- score[row,f$response[row]]+1
    index <- row+n*(seq_len(d)-1)
    W[index,index] <- diag(p[row,seq_len(d)],d)-tcrossprod(p[row,seq_len(d)])
  }
  directions <- c(lapply(seq_len(ncol(f$fixedDesign)),function(k)
    list(a=as.numeric(f$fixedDesign[,k]),Q=matrix(0,ncol(Z),ncol(Z)))),
    lapply(precision$derivatives,function(Q)list(a=numeric(n*d),Q=as.matrix(Q))))
  t(vapply(directions,function(direction){
    modeDerivative <- -solve(H,crossprod(Z,W %*% direction$a)+direction$Q %*% mode)
    etaDerivative <- direction$a+Z %*% modeDerivative
    derivativeW <- function(v){
      result <- matrix(0,n*d,n*d)
      for(row in seq_len(n)){
        index <- row+n*(seq_len(d)-1)
        third <- matrix(0,d,d)
        for(category in seq_len(d+1)){
          centered <- -p[row,seq_len(d)]
          if(category<=d) centered[category] <- centered[category]+1
          third <- third+p[row,category]*tcrossprod(centered)*
            sum(centered*v[index])
        }
        result[index,index] <- third
      }
      result
    }
    direct <- sum(score*direction$a)-drop(crossprod(mode,direction$Q %*% mode))/2+
      sum(diag(solve(as.matrix(precision$value),direction$Q)))/2
    logdet <- sum(diag(solve(H,direction$Q+crossprod(Z,derivativeW(etaDerivative) %*% Z))))/2
    frozenLogdet <- sum(diag(solve(H,direction$Q+crossprod(Z,derivativeW(direction$a) %*% Z))))/2
    c(full=direct-logdet,frozen=direct-frozenLogdet,direct=direct)
  },numeric(3)))
}

test_that("multinomial marginal scores include full mode and determinant dependence", {
  f <- multinomialScoreFixture()
  context <- multinomialModeContext(f)
  precision <- multinomialScorePrecision()
  value <- multinomialScoreEval(context,f,precision)
  expect_true(value$valid)
  expected <- multinomialScoreDense(f,value,precision)
  expect_equal(as.numeric(value$gradient),expected[,"full"],tolerance=1e-10)
  expect_gt(max(abs(expected[,"full"]-expected[,"frozen"])),1e-3)
  expect_gt(max(abs(expected[,"full"]-expected[,"direct"])),1e-3)
  expect_equal(value$work$scoreEvaluations,1L)
  expect_gt(value$work$derivativeSolves,0L)
  expect_equal(value$work$symbolicAnalyses,1L)
  expect_equal(value$work$precisionSymbolicAnalyses,1L)
  beta <- c(.1,-.2,.3)
  theta <- c(.2,-.1,.3,.25)
  evaluate <- function(parameters){
    candidate <- f
    candidate$predictor <- f$predictor+
      matrix(f$fixedDesign %*% parameters[1:3],5,2)
    out <- multinomialScoreEval(context,candidate,
      multinomialScorePrecision(parameters[4:7]),gradient=FALSE)
    expect_true(out$valid)
    expect_false("gradient" %in% names(out))
    out$logLik
  }
  f0 <- f
  f0$predictor <- f$predictor+matrix(f$fixedDesign %*% beta,5,2)
  analytic <- as.numeric(multinomialScoreEval(context,f0,precision)$gradient)
  for(step in c(1e-3,3e-4,1e-4)){
    parameters <- c(beta,theta)
    numerical <- vapply(seq_along(parameters),function(k){
      plus <- minus <- parameters
      plus[k] <- plus[k]+step
      minus[k] <- minus[k]-step
      (evaluate(plus)-evaluate(minus))/(2*step)
    },numeric(1))
    expect_equal(analytic,numerical,tolerance=2e-6)
  }
  direction <- c(.3,-.4,.7,-.2,.5,.1,-.6)
  step <- 1e-4
  numerical <- (evaluate(c(beta,theta)+step*direction)-
    evaluate(c(beta,theta)-step*direction))/(2*step)
  expect_equal(sum(analytic*direction),numerical,tolerance=1e-7)
})

test_that("multinomial precision score proposals reuse factors and recover from failures", {
  f <- multinomialScoreFixture()
  context <- multinomialModeContext(f)
  baseline <- multinomialScoreEval(context,f)
  for(theta in list(c(-.3,.4,.1,.2),c(.1,.2,-.1,-.3),c(.2,-.1,.3,.25))){
    precision <- multinomialScorePrecision(theta)
    shared <- multinomialScoreEval(context,f,precision)
    fresh <- multinomialScoreEval(multinomialModeContext(f),f,precision)
    expect_true(shared$valid)
    expect_equal(shared$logLik,fresh$logLik,tolerance=1e-12)
    expect_equal(shared$gradient,fresh$gradient,tolerance=1e-10)
    expect_equal(shared$work$symbolicAnalyses,1L)
    expect_equal(shared$work$precisionSymbolicAnalyses,1L)
  }
  failed <- multinomialScoreEval(context,f,start=rep(100,3),maxit=1L)
  expect_false(failed$valid)
  expect_false("gradient" %in% names(failed))
  expect_equal(failed$work$scoreEvaluations,4L)
  bad <- list(value=Matrix::Diagonal(3,c(1,-1,1)),derivatives=list())
  expect_warning(expect_error(multinomialScoreEval(context,f,bad),"positive definite"),
    "CHOLMOD warning.*not positive definite")
  recovered <- multinomialScoreEval(context,f)
  expect_equal(recovered$logLik,baseline$logLik,tolerance=1e-12)
  expect_equal(recovered$gradient,baseline$gradient,tolerance=1e-10)
  off <- multinomialScoreEval(context,f,gradient=FALSE)
  expect_equal(off$work$derivativeSolves,recovered$work$derivativeSolves)
  expect_equal(off$work$scoreEvaluations,recovered$work$scoreEvaluations)
})

test_that("multinomial marginal scores transform with category reference and parameter coordinates", {
  f <- multinomialScoreFixture()
  precision <- multinomialScorePrecision()
  original <- multinomialScoreEval(multinomialModeContext(f),f,precision)
  map <- rbind(c(-1,0),c(-1,1))
  changed <- f
  changed$response <- match(f$response,c(3,2,1))
  changed$predictor <- f$predictor %*% t(map)
  stacked <- Matrix::kronecker(Matrix::Matrix(map,sparse=TRUE),Matrix::Diagonal(5))
  changed$design <- stacked %*% f$design
  changed$fixedDesign <- stacked %*% f$fixedDesign
  rebased <- multinomialScoreEval(multinomialModeContext(changed),changed,precision)
  expect_equal(rebased$logLik,original$logLik,tolerance=1e-12)
  expect_equal(rebased$gradient,original$gradient,tolerance=1e-10)
  transform <- matrix(c(1,.2,-.1,0,1,.3),3,2)
  changed <- f
  changed$fixedDesign <- Matrix::Matrix(f$fixedDesign %*% transform,sparse=TRUE)
  transformed <- multinomialScoreEval(multinomialModeContext(changed),changed,precision)
  expected <- c(as.numeric(crossprod(transform,original$gradient[1:3])),
    original$gradient[4:7])
  expect_equal(as.numeric(transformed$gradient),expected,tolerance=1e-10)
})

test_that("multinomial scores preserve binary, independent and saturated limits", {
  f <- multinomialScoreFixture()
  f$predictor <- f$predictor[,1,drop=FALSE]
  f$response <- ifelse(f$response==1,1,2)
  f$design <- f$design[1:5,]
  f$fixedDesign <- f$fixedDesign[1:5,c(1,3)]
  precision <- multinomialScorePrecision()
  value <- multinomialScoreEval(multinomialModeContext(f),f,precision)
  expected <- multinomialScoreDense(f,value,precision)
  expect_true(value$valid)
  expect_equal(as.numeric(value$gradient),expected[,"full"],tolerance=1e-10)
  for(k in 1:2){
    h <- 1e-4
    plus <- minus <- f
    plus$predictor <- f$predictor+matrix(f$fixedDesign[,k]*h,5,1)
    minus$predictor <- f$predictor-matrix(f$fixedDesign[,k]*h,5,1)
    scalar <- function(candidate){
      sommer:::sommer_laplace_glmm_eval(
        Matrix::Matrix(matrix(0,5,1),sparse=TRUE),candidate$design,precision$value,
        as.numeric(candidate$response==1),as.numeric(candidate$predictor),0,rep(0,3),
        "binomial",as.numeric(determinant(precision$value,logarithm=TRUE)$modulus),
        100L,1e-11)$logLik
    }
    expect_equal(as.numeric(value$gradient)[k],(scalar(plus)-scalar(minus))/(2*h),
      tolerance=1e-7)
  }
  f <- multinomialScoreFixture()
  f$design <- Matrix::Matrix(0,10,3,sparse=TRUE)
  independent <- multinomialScoreEval(multinomialModeContext(f),f,precision)
  kernel <- sommer:::sommer_multinomial_family_eval(f$response,f$predictor,matrix(0,5,2))
  expect_equal(as.numeric(independent$gradient),
    c(as.numeric(crossprod(f$fixedDesign,as.numeric(kernel$score))),rep(0,4)),
    tolerance=1e-12)
  f <- multinomialScoreFixture()
  f$predictor <- matrix(rep(c(1000,-1000),each=5),5,2)
  context <- multinomialModeContext(f)
  saturated <- multinomialScoreEval(context,f,precision)
  expect_true(saturated$valid)
  expect_true(all(is.finite(saturated$gradient)))
  expected <- multinomialScoreDense(f,saturated,precision)
  expect_equal(as.numeric(saturated$gradient),expected[,"full"],tolerance=1e-10)
  expect_equal(saturated$work$symbolicAnalyses,1L)
})

test_that("multinomial score derivatives and bounded solve workspaces are validated", {
  f <- multinomialScoreFixture()
  context <- multinomialModeContext(f)
  precision <- multinomialScorePrecision()
  bad <- f
  bad$fixedDesign <- f$fixedDesign[-1,]
  expect_error(multinomialScoreEval(context,bad),"category-major")
  bad <- f
  bad$fixedDesign[1,1] <- Inf
  expect_error(multinomialScoreEval(context,bad),"finite")
  for(derivative in list(Matrix::Diagonal(2),
    Matrix::Matrix(matrix(c(1,2,0,1),2),sparse=TRUE),
    Matrix::Diagonal(3,c(1,Inf,1)))){
    badPrecision <- precision
    badPrecision$derivatives <- list(derivative)
    expect_error(multinomialScoreEval(context,f,badPrecision),"precision derivatives")
  }
  badPrecision <- precision
  badPrecision$value <- Matrix::Diagonal(2)
  expect_error(multinomialScoreEval(context,f,badPrecision),"dimensions")
  zeroPrecision <- precision
  zeroPrecision$derivatives <- list(Matrix::Matrix(0,3,3,sparse=TRUE))
  value <- multinomialScoreEval(context,f,zeroPrecision)
  expect_equal(as.numeric(value$gradient)[4],0)
  noFixed <- f
  noFixed$fixedDesign <- Matrix::Matrix(0,10,0,sparse=TRUE)
  onlyPrecision <- multinomialScoreEval(context,noFixed,precision)
  expect_length(onlyPrecision$gradient,4)
  none <- precision
  none$derivatives <- list()
  noParameters <- multinomialScoreEval(context,noFixed,none)
  expect_length(noParameters$gradient,0)
  small <- multinomialModeContext(f,memory=.011)
  expect_true(multinomialScoreEval(small,f,gradient=FALSE)$valid)
  expect_error(multinomialScoreEval(small,f),"score workspace")
  set.seed(8123)
  before <- .Random.seed
  invisible(multinomialScoreEval(context,f))
  expect_identical(before,.Random.seed)
})

test_that("multinomial solve panels preserve scores beyond one panel", {
  set.seed(3192)
  f <- multinomialScoreFixture()
  q <- 39L
  Z <- matrix(rnorm(10*q,sd=.2),10,q)
  Z[abs(Z)<.12] <- 0
  f$design <- Matrix::Matrix(Z,sparse=TRUE)
  basis <- matrix(rnorm(q*q,sd=.1),q,q)
  Q <- crossprod(basis)+diag(1,q)
  direction <- Matrix::Matrix(crossprod(matrix(rnorm(q*2,sd=.1),2,q)),sparse=TRUE)
  precision <- list(value=Matrix::Matrix(Q,sparse=TRUE),derivatives=list(direction))
  f$precision <- precision$value
  full <- multinomialScoreEval(multinomialModeContext(f),f,precision)
  expected <- multinomialScoreDense(f,full,precision)
  expect_equal(as.numeric(full$gradient),expected[,"full"],tolerance=1e-10)
  structuralPairs <- sum(vapply(1:5,function(row)
    sum(colSums(abs(Z[c(row,row+5),,drop=FALSE]))>0)^2,numeric(1)))
  nnzZ <- sum(Z != 0)
  nnzQ <- sum(Q != 0)
  workspace <- 8*(16*5*2^2+24*5*2+16*q+8*nnzZ+8*nnzQ)+
    64*(structuralPairs+nnzQ)
  derivativeBytes <- 8*(12*10+8*2^2+12*q+4*nnzZ+3+1)+
    32*sum(as.matrix(direction) != 0)
  columnBytes <- 8*(4*q+2)
  memory <- (workspace+derivativeBytes+1.5*columnBytes)/1048576
  panelled <- multinomialScoreEval(multinomialModeContext(f,memory),f,precision)
  expect_equal(panelled$gradient,full$gradient,tolerance=1e-10)
  expect_gt(panelled$work$derivativeSolves,full$work$derivativeSolves)
  h <- 1e-4
  plus <- minus <- precision
  plus$value <- precision$value+h*direction
  minus$value <- precision$value-h*direction
  numeric <- (multinomialScoreEval(multinomialModeContext(f),f,plus,FALSE)$logLik-
    multinomialScoreEval(multinomialModeContext(f),f,minus,FALSE)$logLik)/(2*h)
  expect_equal(as.numeric(full$gradient)[4],numeric,tolerance=1e-7)

  n <- 4L
  d <- 35L
  q <- 3L
  f <- list(response=c(1,12,30,36),predictor=matrix(rnorm(n*d),n,d),
    design=Matrix::Matrix(matrix(rnorm(n*d*q,sd=.2),n*d,q),sparse=TRUE),
    fixedDesign=Matrix::Matrix(matrix(rnorm(n*d*2),n*d,2),sparse=TRUE),
    precision=Matrix::Diagonal(q))
  precision <- list(value=f$precision,derivatives=list())
  value <- multinomialScoreEval(multinomialModeContext(f),f,precision)
  expect_true(value$valid)
  expected <- multinomialScoreDense(f,value,precision)
  expect_equal(as.numeric(value$gradient),expected[,"full"],tolerance=1e-10)
})

test_that("analytic multinomial scores drive the same fixed-effect optimum as numerical scores", {
  f <- multinomialScoreFixture()
  f$fixedDesign <- f$fixedDesign[,1:2]
  precision <- multinomialScorePrecision()
  precision$derivatives <- list()
  context <- multinomialModeContext(f)
  evaluate <- function(beta,gradient){
    candidate <- f
    candidate$predictor <- f$predictor+matrix(f$fixedDesign %*% beta,5,2)
    multinomialScoreEval(context,candidate,precision,gradient)
  }
  numerical <- optim(c(0,0),function(beta)-evaluate(beta,FALSE)$logLik,
    method="BFGS",control=list(reltol=1e-11))
  analytic <- optim(c(0,0),function(beta)-evaluate(beta,FALSE)$logLik,
    gr=function(beta)-as.numeric(evaluate(beta,TRUE)$gradient),
    method="BFGS",control=list(reltol=1e-11))
  expect_equal(analytic$convergence,0L)
  expect_equal(numerical$convergence,0L)
  expect_equal(analytic$par,numerical$par,tolerance=1e-5)
  expect_equal(analytic$value,numerical$value,tolerance=1e-10)
  expect_lt(max(abs(evaluate(analytic$par,TRUE)$gradient)),1e-5)
})

test_that("precision sparsity changes refresh symbolic analysis without stale scores", {
  f <- multinomialScoreFixture()
  f$design <- Matrix::Matrix(rbind(diag(3),matrix(0,7,3)),sparse=TRUE)
  f$precision <- Matrix::Diagonal(3)
  context <- multinomialModeContext(f)
  derivative <- Matrix::Matrix(matrix(c(0,.3,0,.3,0,.2,0,.2,0),3),sparse=TRUE)
  values <- list()
  for(theta in c(0,.4,0)){
    precision <- list(value=Matrix::drop0(f$precision+theta*derivative),
      derivatives=list(derivative))
    out <- multinomialScoreEval(context,f,precision)
    fresh <- multinomialScoreEval(multinomialModeContext(f),f,precision)
    expect_true(out$valid)
    expect_equal(out$gradient,fresh$gradient,tolerance=1e-12)
    expect_equal(out$logLik,fresh$logLik,tolerance=1e-12)
    values[[length(values)+1L]] <- out
    h <- 1e-4
    plus <- minus <- precision
    plus$value <- precision$value+h*derivative
    minus$value <- precision$value-h*derivative
    numerical <- (multinomialScoreEval(context,f,plus,FALSE)$logLik-
      multinomialScoreEval(context,f,minus,FALSE)$logLik)/(2*h)
    expect_equal(as.numeric(out$gradient)[4],numerical,tolerance=1e-7)
  }
  expect_gt(values[[2]]$work$symbolicAnalyses,values[[1]]$work$symbolicAnalyses)
  expect_gt(values[[3]]$work$symbolicAnalyses,values[[2]]$work$symbolicAnalyses)
  expect_gt(values[[2]]$work$precisionSymbolicAnalyses,
    values[[1]]$work$precisionSymbolicAnalyses)
  expect_equal(values[[1]]$gradient,values[[3]]$gradient,tolerance=1e-12)
})
