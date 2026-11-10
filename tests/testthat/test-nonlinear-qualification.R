nonlinearQualificationData <- function(groups=16L,points=12L,curve="exponential",
    randomCurve=FALSE,seed=652L){
  set.seed(seed)
  data <- expand.grid(x=seq(0,5,length.out=points),group=factor(seq_len(groups)))
  u <- rnorm(groups,sd=.2)[data$group]
  mean <- switch(curve,exponential=3*exp(-.65*data$x),
    logistic=3*plogis(1.2*(data$x-2.5)),saturation=3*data$x/(.8+data$x))
  data$y <- if(randomCurve) 3*exp(-.65*exp(u)*data$x)+rnorm(nrow(data),sd=.1) else
    mean+u+rnorm(nrow(data),sd=.1)
  data
}

nonlinearQualificationControl <- function(curve="exponential",profile="off"){
  start <- switch(curve,exponential=c(asymptote=3,rate=.65),
    logistic=c(asymptote=3,rate=1.2,midpoint=2.5),saturation=c(asymptote=3,half=.8))
  mmes_control(task=list(objective="ml",nonlinear=list(start=start)),
    convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE,profile=profile))
}

test_that("sparse admission scales by connected components before matrix assembly",{
  independent <- Matrix::Diagonal(10000L)
  bytes <- sommer:::sommer_nonlinear_sparse_admission(independent,independent,independent,50)
  expect_lt(bytes,50*1048576)
  expect_gt(bytes,0)
  expect_error(sommer:::sommer_nonlinear_sparse_admission(independent,independent,
    independent,1),"before Hessian assembly")
  chain <- Matrix::bandSparse(200L,k=c(-1,0,1),
    diagonals=list(rep(-.1,199),rep(1,200),rep(-.1,199)))
  diagonal <- Matrix::Diagonal(200L)
  expect_error(sommer:::sommer_nonlinear_sparse_admission(diagonal,chain,diagonal,1),
    "factor-fill bound")
  coupled <- Matrix::sparseMatrix(i=rep(1L,200),j=1:200,x=1,dims=c(1,200))
  expect_error(sommer:::sommer_nonlinear_sparse_admission(coupled,diagonal,
    Matrix::Diagonal(1L),1),"factor-fill bound")
  expect_error(sommer:::sommer_nonlinear_sparse_admission(diagonal,diagonal,chain,1),
    "factor-fill bound")
  expect_error(sommer:::sommer_nonlinear_sparse_admission(diagonal,diagonal,diagonal,0),
    "Invalid")
  expect_error(sommer:::sommer_nonlinear_sparse_admission(diagonal,diagonal,
    Matrix::Diagonal(2L)),"Invalid")
})

test_that("profiling is opt-in, deterministic, serializable and inference-free",{
  data <- nonlinearQualificationData()
  fit <- mmes(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearQualificationControl())
  expect_null(fit$nonlinearDiagnostics$profile)
  timed <- mmes(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearQualificationControl(profile="force"))
  expect_equal(coef(timed),coef(fit),tolerance=1e-12)
  expect_equal(timed$llik,fit$llik,tolerance=1e-12)
  diagnostics <- timed$nonlinearDiagnostics
  expect_equal(sum(diagnostics$profile$seconds),diagnostics$profile$totalSeconds,
    tolerance=1e-8)
  expect_true(all(diagnostics$profile$seconds>=0))
  expect_true(is.finite(diagnostics$profile$marginalEvaluationSeconds))
  expect_equal(diagnostics$optimizer$iterations,fit$optimizer$iterations)
  expect_identical(timed$executionPlan$optimizer$gradient,"analytic")
  expect_identical(timed$executionPlan$starts$curve,"supplied")
  expect_equal(timed$executionPlan$dimensions["randomFull"],16,ignore_attr=TRUE)
  before <- .Random.seed
  profile <- profile_nonlinear(timed,replications=1L)
  expect_identical(before,.Random.seed)
  expect_named(profile,c("table","contextSeconds","fitProfile","dimensions","work",
    "optimizer","starts","memory","likelihood","converged","boundaryParameters",
    "environment","note"))
  expect_equal(nrow(profile$table),6L)
  expect_true(all(profile$table$seconds>=0))
  expect_true(profile$memory$fitBytes>0)
  expect_null(timed$nonlinearInference)
  expect_false(timed$nonlinearDiagnostics$uncertaintyComputed)
  expect_identical(unserialize(serialize(profile,NULL)),profile)
  expect_error(profile_nonlinear(timed,replications=0),"positive")
  expect_error(profile_nonlinear(list()),"nonlinear mmes")
  expect_error(profile_nonlinear(timed,inference=NA),"inference")
  expect_error(profile_nonlinear(timed,
    control=mmes_control(execution=list(denseMemoryMB=.0001))),"workspace")
  withInference <- profile_nonlinear(timed,replications=1L,inference=TRUE)
  expect_true("postFitInference" %in% withInference$table$phase)
  expect_null(timed$nonlinearInference)
})

for(curve in c("exponential","logistic","saturation")){
  test_that(paste(curve,"Gaussian curve matches independent dense maximum likelihood"),{
    starts <- list(exponential=c(logA=log(3),logRate=log(.65),location=0),
      logistic=c(logA=log(3),logRate=log(1.2),midpoint=2.5,location=0),
      saturation=c(logA=log(3),logHalf=log(.8),location=0))
    data <- nonlinearQualificationData(curve=curve)
    fit <- mmes(y~1,random=~group,data=data,nonlinear=nonlinear_curve("x",curve),
      control=nonlinearQualificationControl(curve))
    k <- length(starts[[curve]])-1L
    meanValue <- function(at) switch(curve,
      exponential=exp(at[1])*exp(-exp(at[2])*data$x),
      logistic=exp(at[1])*plogis(exp(at[2])*(data$x-at[3])),
      saturation=exp(at[1])*data$x/(exp(at[2])+data$x))+at[k+1L]
    sameGroup <- outer(data$group,data$group,"==")*1
    objective <- function(at){
      variance <- exp(at[k+2L])*sameGroup+diag(exp(at[k+3L]),nrow(data))
      factor <- chol(variance)
      residual <- data$y-meanValue(at)
      whitened <- forwardsolve(t(factor),residual)
      .5*(nrow(data)*log(2*pi)+2*sum(log(diag(factor)))+sum(whitened^2))
    }
    reference <- optim(c(starts[[curve]],log(.04),log(.01)),objective,
      method="BFGS",control=list(maxit=500L,reltol=1e-12))
    variance <- unname(exp(tail(reference$par,2L)))
    covariance <- variance[1]*sameGroup+diag(variance[2],nrow(data))
    prediction <- meanValue(reference$par)+as.numeric(
      variance[1]*sameGroup%*%solve(covariance,data$y-meanValue(reference$par)))
    expect_equal(reference$convergence,0L)
    expect_true(fit$convergence)
    expect_equal(fit$llik,-reference$value,tolerance=2e-6)
    expect_equal(as.numeric(fitted(fit)),prediction,tolerance=2e-5)
    referenceCurve <- reference$par[seq_len(k)]
    referenceCurve[seq_len(2)] <- exp(referenceCurve[seq_len(2)])
    expect_equal(unname(fit$nonlinearParameters),unname(referenceCurve),tolerance=2e-4)
    expect_equal(unname(fit$covPar[[1L]]),variance[1],tolerance=2e-4)
    expect_equal(unname(fit$dispersion),variance[2],tolerance=2e-4)
  })
}

test_that("random curve recovery and predictions agree with generating truth",{
  data <- nonlinearQualificationData(groups=30L,points=16L,randomCurve=TRUE)
  fit <- mmes(y~0,data=data,nonlinear=nonlinear_curve("x",random=list(rate=~group)),
    control=nonlinearQualificationControl())
  expect_true(fit$convergence)
  expect_lt(sqrt(mean((fitted(fit)-data$y)^2)),.15)
  expect_equal(unname(fit$nonlinearParameters),c(3,.65),tolerance=.08)
  expect_match(fit$likelihoodMethod,"Laplace")
  expect_true(all(is.finite(profile_nonlinear(fit,replications=1L)$table$seconds)))
})

test_that("native family profiling preserves likelihood and no inference work",{
  data <- nonlinearQualificationData()
  set.seed(738)
  data$y <- rpois(nrow(data),exp(1.2*exp(-.65*data$x)))
  control <- nonlinearQualificationControl(profile="force")
  control$task$nonlinear$start <- c(asymptote=1.2,rate=.65)
  fit <- mmes(y~0,data=data,family=poisson(),nonlinear=nonlinear_curve("x",
    random=list(rate=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE))),
    control=control)
  expect_true(fit$convergence)
  profile <- profile_nonlinear(fit,replications=1L)
  expect_identical(profile$likelihood,"Laplace family ML")
  expect_equal(nrow(profile$table),6L)
  expect_gt(fit$executionPlan$mode$iterations,0)
  expect_true(all(is.finite(profile$table$seconds)))
  expect_null(fit$nonlinearInference)
})
