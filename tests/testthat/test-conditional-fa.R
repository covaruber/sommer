.conditional_fixture <- function(){
  set.seed(162)
  n <- 7L
  data <- expand.grid(id=factor(paste0("g", seq_len(n))), rep=1:2,
    trait=factor(paste0("t", 1:3)))
  data$record <- interaction(data$id, data$rep)
  data$x <- rep(c(-1, 1), each=n, length.out=nrow(data))
  data$y <- rnorm(nrow(data)) + as.integer(data$id)/4 +
    as.integer(data$trait)/3
  A <- outer(seq_len(n), seq_len(n), function(i,j) 0.4^abs(i-j))
  precision <- solve(A)
  dimnames(precision) <- list(levels(data$id), levels(data$id))
  attr(precision, "inverse") <- TRUE
  list(data=data, precision=precision, A=A)
}

.conditional_control <- function(iterations=4L, panel=8L){
  mmes_control(update=list(method="conditional-fa"),
    convergence=list(maxIters=iterations),
    output=list(covariance="none", coefficientMatrix=FALSE, panelSize=panel),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

test_that("conditional FA handles incomplete replicated records with exact REML and BLUP", {
  fixture <- .conditional_fixture()
  data <- fixture$data[-c(2, 15, 31, 38),]
  data <- data[rev(seq_len(nrow(data))),]
  precision <- fixture$precision
  Lg <- matrix(c(.7, -.3, .5), 3, 1)
  Le <- matrix(c(.4, .6, -.2), 3, 1)
  pg <- c(.3, .5, .4); pe <- c(.7, .6, .8)
  G <- tcrossprod(Lg) + diag(pg); R <- tcrossprod(Le) + diag(pe)
  sg <- fam(data$trait, 1, loadings=Lg, specific=pg, fixed=rep(TRUE, 5L))
  se <- fam(data$trait, 1, loadings=Le, specific=pe, fixed=rep(TRUE, 5L))
  fit <- mmes(y~trait+trait:x, random=~vsm(sg, ism(id), Gu=precision,
    sigma2=G[1,1], fixedSigma2=TRUE),
    rcov=~vsm(se, ism(record), sigma2=R[1,1], fixedSigma2=TRUE),
    data=data, control=.conditional_control())
  X <- model.matrix(~trait+trait:x, data)
  ti <- as.integer(data$trait); ii <- as.integer(data$id)
  Vg <- G[ti,ti] * fixture$A[ii,ii]
  Vr <- R[ti,ti] * outer(data$record, data$record, `==`)
  V <- Vg + Vr
  ViX <- solve(V, X)
  b <- solve(crossprod(X, ViX), crossprod(ViX, data$y))
  Py <- solve(V, data$y - X %*% b)
  ll <- -.5 * (as.numeric(determinant(V, logarithm=TRUE)$modulus) +
    as.numeric(determinant(crossprod(X, ViX), logarithm=TRUE)$modulus) +
    sum(data$y * Py))
  expect_true(fit$convergence)
  expect_equal(as.numeric(fit$b), as.numeric(b), tolerance=1e-9)
  expect_equal(as.numeric(fitted(fit)), as.numeric(X%*%b + Vg%*%Py), tolerance=1e-9)
  expect_equal(tail(fit$conditionalFA$rawLikelihood,1L), ll, tolerance=1e-9)
  expect_equal(dim(fit$llik),c(1L,1L))
  expect_equal(tail(as.numeric(fit$llik),1L), ll+(nrow(data)-ncol(X))*log(sd(data$y)),
    tolerance=1e-9)
  expect_equal(unname(covmatrix_mmes(fit,1L,se=FALSE)$covariance), G, tolerance=1e-10)
  expect_equal(unname(covmatrix_mmes(fit,2L,se=FALSE)$covariance), R, tolerance=1e-10)
  expect_null(fit$C)
  expect_null(fit$theta_se)
  expect_true(fit$executionPlan$admission$conditionalFA$active)
  expect_true(all(is.finite(scores_mmes(fit, term=names(fit$uList)[1L]))))
  expect_output(print(summary(fit)), "Covariance information omitted")
  expect_error(covmatrix_mmes(fit,1L), "Covariance information was omitted")
  expect_error(postPEV(fit), "Coefficient matrix C was not retained")
  ordinary <- mmes(y~trait+trait:x, random=~vsm(sg,ism(id),Gu=precision,
    sigma2=G[1,1],fixedSigma2=TRUE),
    rcov=~vsm(se,ism(record),sigma2=R[1,1],fixedSigma2=TRUE),data=data,
    control=mmes_control(diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_equal(as.numeric(fit$b),as.numeric(ordinary$b),tolerance=1e-9)
  expect_equal(as.numeric(fit$u),as.numeric(ordinary$u),tolerance=1e-9)
  expect_equal(tail(as.numeric(fit$llik),1L),
    tail(as.numeric(ordinary$llik),1L),tolerance=1e-9)
})

test_that("independent and shared conditional updates increase observed REML", {
  fixture <- .conditional_fixture()
  data <- fixture$data
  data$y[c(3, 18, 39)] <- NA
  precision <- fixture$precision
  for(shared in list(NULL, "axes")){
    sg <- fam(data$trait, 1, shared=shared)
    se <- fam(data$trait, 1, shared=shared)
    expect_warning(fit <- mmes(y~trait, random=~vsm(sg, ism(id), Gu=precision),
      rcov=~vsm(se, ism(record)), data=data, control=.conditional_control()),
      "did not converge")
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-8))
    expect_true(all(is.finite(fitted(fit))))
    expect_equal(fit$engineDiagnostics$conditionalFA$symbolicAnalyses, 1L)
    expect_equal(fit$engineDiagnostics$conditionalFA$numericFactorizations, 5L)
    expect_lte(fit$engineDiagnostics$conditionalFA$panelColumns, 8L)
    state <- fit$conditionalFA$state
    if(!is.null(shared)){
      expect_equal(state$weightG+state$weightE, 1)
      expectedG <- tcrossprod(sweep(state$loading,2,sqrt(state$weightG),`*`)) +
        diag(state$specificG)
      expectedE <- tcrossprod(sweep(state$loading,2,sqrt(state$weightE),`*`)) +
        diag(state$specificE)
      expect_equal(unname(covmatrix_mmes(fit,1L,se=FALSE)$covariance), expectedG,
        tolerance=1e-9)
      expect_equal(unname(covmatrix_mmes(fit,2L,se=FALSE)$covariance), expectedE,
        tolerance=1e-9)
    }
  }
})

test_that("diagonal and homogeneous residuals use the same conditional blocks", {
  fixture <- .conditional_fixture()
  precision <- fixture$precision
  for(residual in list(~units, ~vsm(dsm(trait),ism(record)),
                       ~vsm(ism(trait),ism(record)))){
    expect_warning(fit <- mmes(y~trait,
      random=~vsm(fam(trait,1),ism(id),Gu=precision),
      rcov=residual, data=fixture$data, control=.conditional_control(2L)),
      "did not converge")
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-8))
    expect_true(all(is.finite(fitted(fit))))
  }
})

test_that("factor score regression uses the low-rank contraction on both scales", {
  fixture <- .conditional_fixture()
  fit <- mmes(y~trait, random=~vsm(fam(trait,1),ism(id)), data=fixture$data,
    control=mmes_control(convergence=list(maxIters=2L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  term <- names(fit$uList)[1L]
  for(scale in c(FALSE, TRUE)) for(rotate in c(FALSE, TRUE)){
    loadings <- loadings_mmes(fit, term, varianceScale=scale, rotation=rotate)$loadings
    expected <- fit$uList[[term]] %*% solve(fit$theta[[term]], loadings)
    expect_equal(scores_mmes(fit, term, varianceScale=scale, rotation=rotate),
      expected, tolerance=1e-9)
  }
})

test_that("shared loadings and conditional controls fail explicitly outside their scope", {
  fixture <- .conditional_fixture()
  data <- fixture$data
  precision <- fixture$precision
  expect_error(fam(data$trait, shared=TRUE), "loading-group")
  expect_error(mmes_control(update=list(method="unknown")), "update\\$method")
  expect_error(mmes(y~trait, random=~vsm(fam(trait,1,shared="a"),ism(id)),
    data=data, control=mmes_control(diagnostics=list(verbose=FALSE))), "requires update")
  expect_error(mmes(y~trait, random=~vsm(fam(trait,1),ism(id)),
    data=data, control=mmes_control(update=list(method="conditional-fa"))),
    "output=list")
  expect_error(mmes(y~trait, random=~id, data=data,
    control=.conditional_control()), "one fam")
  expect_error(mmes(y~trait, random=~vsm(fam(trait,1,shared="a"),ism(id)),
    rcov=~vsm(fam(trait,1,shared="b"),ism(record)), data=data,
    control=.conditional_control()), "matching")
  expect_error(mmes(y~trait, random=~vsm(fam(trait,1),ism(id),Gu=precision),
    rcov=~vsm(fam(trait,1),ism(record)), data=data, W=diag(nrow(data))*2,
    control=.conditional_control()), "W/weights")
})

test_that("conditional EM converges to an independently optimized restricted likelihood", {
  fixture <- .conditional_fixture()
  data <- fixture$data[-c(2,15,31,38),]
  precision <- fixture$precision
  L <- matrix(c(.7,-.3,.5),3,1)
  G <- tcrossprod(L)+diag(c(.3,.5,.4))
  sg <- fam(data$trait,1,loadings=L,specific=c(.3,.5,.4),
    fixed=rep(TRUE,5L))
  control <- .conditional_control(300L)
  control$convergence$likelihoodTol <- 1e-9
  control$convergence$parameterTol <- 1e-7
  fit <- mmes(y~trait+trait:x,
    random=~vsm(sg,ism(id),Gu=precision,sigma2=G[1,1],fixedSigma2=TRUE),
    rcov=~units,data=data,control=control)
  X <- model.matrix(~trait+trait:x,data)
  Vg <- G[as.integer(data$trait),as.integer(data$trait)] *
    fixture$A[as.integer(data$id),as.integer(data$id)]
  evaluate <- function(logVariance, coefficients=FALSE){
    V <- Vg+diag(exp(logVariance),nrow(data))
    ViX <- solve(V,X)
    b <- solve(crossprod(X,ViX),crossprod(ViX,data$y))
    if(coefficients) return(as.numeric(b))
    -.5*(as.numeric(determinant(V,logarithm=TRUE)$modulus) +
      as.numeric(determinant(crossprod(X,ViX),logarithm=TRUE)$modulus) +
      sum(data$y*solve(V,data$y-X%*%b)))
  }
  optimum <- optimize(evaluate,c(-8,5),maximum=TRUE,tol=1e-9)
  expect_true(fit$convergence)
  expect_equal(tail(fit$conditionalFA$rawLikelihood,1L),optimum$objective,
    tolerance=1e-9)
  expect_equal(unname(covmatrix_mmes(fit,2L,se=FALSE)$covariance)[1,1],
    exp(optimum$maximum),tolerance=1e-6)
  expect_equal(as.numeric(fit$b),evaluate(optimum$maximum,TRUE),tolerance=1e-6)
  expect_equal(fit$AIC,-2*tail(as.numeric(fit$llik),1L)+2*ncol(X))
  expect_equal(fit$BIC,-2*tail(as.numeric(fit$llik),1L)+log(nrow(data))*ncol(X))
  controlAI <- mmes_control(representation=list(equations="conditional"),
    convergence=list(maxIters=100L,likelihoodTol=1e-9,parameterTol=1e-7),
    output=list(covariance="none",coefficientMatrix=FALSE),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  fitAI <- mmes(y~trait+trait:x,
    random=~vsm(sg,ism(id),Gu=precision,sigma2=G[1,1],fixedSigma2=TRUE),
    rcov=~units,data=data,control=controlAI)
  expect_true(fitAI$convergence)
  expect_equal(tail(fitAI$conditionalFA$rawLikelihood,1L),optimum$objective,
    tolerance=1e-9)
  expect_equal(as.numeric(fitted(fitAI)),as.numeric(fitted(fit)),tolerance=1e-6)
  expect_equal(fitAI$conditionalFA$optimizer,"ai")
  expect_lt(length(fitAI$llik),length(fit$llik))
  expect_lte(fitAI$conditionalFA$ai$scoreNorm,1e-7)
})

test_that("unequal FA ranks and homogeneous specifics preserve marginal descriptors", {
  set.seed(271)
  data <- expand.grid(id=factor(seq_len(6)),rep=1:2,trait=factor(seq_len(4)))
  data$record <- interaction(data$id,data$rep)
  data$y <- rnorm(nrow(data))
  data$y[c(3,17)] <- NA
  expect_warning(fit <- mmes(y~trait,
    random=~vsm(fam(trait,2,variance="homogeneous"),ism(id)),
    rcov=~vsm(fam(trait,1,variance="homogeneous"),ism(record)),
    data=data,control=.conditional_control(3L)), "did not converge")
  state <- fit$conditionalFA$state
  expect_equal(dim(state$genetic),c(4L,2L))
  expect_equal(dim(state$residual),c(4L,1L))
  expect_equal(state$specificG,rep(state$specificG[1L],4))
  expect_equal(state$specificE,rep(state$specificE[1L],4))
  expect_equal(unname(covmatrix_mmes(fit,1L,se=FALSE)$covariance),
    tcrossprod(state$genetic)+diag(state$specificG),tolerance=1e-9)
  expect_true(all(diff(as.numeric(fit$llik))>=-1e-8))
})
