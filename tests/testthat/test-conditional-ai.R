.conditional_ai_fixture <- function(shared=NULL, homogeneous=FALSE){
  set.seed(613)
  data <- expand.grid(id=factor(seq_len(6)),rep=seq_len(3),
    trait=factor(paste0("t",seq_len(4))))
  data$record <- interaction(data$id,data$rep)
  data$x <- rnorm(nrow(data))
  data$y <- rnorm(nrow(data))+.3*data$x
  data <- data[-c(3,21,50),]
  precision <- solve(outer(seq_len(6),seq_len(6),function(i,j) .3^abs(i-j)))
  dimnames(precision) <- list(levels(data$id),levels(data$id))
  attr(precision,"inverse") <- TRUE
  L <- cbind(c(.8,.4,-.5,.6),c(0,.7,.2,-.3))
  variance <- if(homogeneous) "homogeneous" else "heterogeneous"
  sg <- fam(data$trait,2,loadings=L,specific=rep(.5,4),
    shared=shared,variance=variance)
  se <- fam(data$trait,2,loadings=L*.7,specific=rep(.8,4),
    shared=shared,variance=variance)
  control <- mmes_control(representation=list(equations="conditional"),
    output=list(inverse=0,covariance="none",coefficientMatrix=FALSE,setup=TRUE),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  setup <- mmes(y~trait+x,random=~vsm(sg,ism(id),Gu=precision),
    rcov=~vsm(se,ism(record)),data=data,control=control)
  model <- sommer:::.conditional_fa_prepare(setup$X,setup$Z,setup$Ai,setup$yvar,
    setup$residualBlock,setup$residualIndex,setup$covStruct)
  list(model=model,engine=sommer:::.mmes_compile_control(control),data=data)
}

.conditional_ai_covariance <- function(model,state){
  G <- if(is.null(model$shared)) tcrossprod(state$genetic) else
    tcrossprod(sweep(state$loading,2,sqrt(state$weightG),`*`))
  E <- if(is.null(model$shared)) tcrossprod(state$residual) else
    tcrossprod(sweep(state$loading,2,sqrt(state$weightE),`*`))
  G <- G+diag(state$specificG)
  E <- E+diag(state$specificE)
  A <- solve(as.matrix(model$Ai))
  G[model$trait,model$trait]*A[model$id,model$id]+
    E[model$trait,model$trait]*outer(model$block,model$block,`==`)
}

test_that("conditional REML scores and AI match independent marginal derivatives",{
  for(shared in list(NULL,"axes")) for(homogeneous in c(FALSE,TRUE)){
    fixture <- .conditional_ai_fixture(shared,homogeneous)
    model <- fixture$model; state <- model$state; engine <- fixture$engine
    posterior <- sommer:::.conditional_fa_posterior(model,state,NULL,engine,TRUE)
    ai <- sommer:::.conditional_fa_ai(model,state,posterior,engine)
    V <- .conditional_ai_covariance(model,state)
    Vi <- solve(V)
    X <- as.matrix(model$X)
    P <- Vi-Vi%*%X%*%solve(crossprod(X,Vi%*%X),t(X)%*%Vi)
    Py <- as.numeric(P%*%model$y)
    coordinates <- ai$coordinates
    derivatives <- lapply(seq_along(coordinates$values),function(j){
      plus <- minus <- coordinates$values
      plus[j] <- plus[j]+1e-5
      minus[j] <- minus[j]-1e-5
      (.conditional_ai_covariance(model,sommer:::.conditional_fa_decode(
        model,state,coordinates,plus))-
        .conditional_ai_covariance(model,sommer:::.conditional_fa_decode(
          model,state,coordinates,minus)))/2e-5
    })
    expectedScore <- vapply(derivatives,function(dV)
      .5*(sum(Py*as.numeric(dV%*%Py))-sum(P*t(dV))),numeric(1))
    B <- vapply(derivatives,function(dV) as.numeric(dV%*%Py),numeric(length(Py)))
    expect_equal(ai$score,expectedScore,tolerance=1e-7)
    expect_equal(ai$information,.5*crossprod(B,P%*%B),tolerance=1e-7)
    expect_equal(ai$coordinates$values,
      sommer:::.conditional_fa_coordinates(model,
        sommer:::.conditional_fa_decode(model,state,coordinates,coordinates$values))$values)
    limited <- engine
    limited$denseMemoryMB <- 1e-6
    expect_error(sommer:::.conditional_fa_ai(model,state,posterior,limited),
      "workspace exceeds")
  }
})

test_that("conditional factor solves reuse the current numeric factorization",{
  fixture <- .conditional_ai_fixture()
  model <- fixture$model
  posterior <- sommer:::.conditional_fa_posterior(model,model$state,NULL,
    fixture$engine,TRUE)
  H <- sommer:::.conditional_fa_system(model,model$state)$H
  rhs <- matrix(seq_len(model$d*2)/model$d,model$d,2)
  solved <- sommer:::mme_conditional_solve_cpp(posterior$workspace,rhs,10)
  expect_equal(solved,unname(as.matrix(solve(H,rhs))),tolerance=1e-9)
  expect_error(sommer:::mme_conditional_solve_cpp(posterior$workspace,rhs,1e-6),
    "budget exceeded")
})

test_that("conditional representation selects AI independently of legacy EM",{
  fixture <- .conditional_ai_fixture("axes")
  control <- mmes_control(representation=list(equations="conditional"),
    convergence=list(maxIters=8L),
    output=list(covariance="none",coefficientMatrix=FALSE),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  expect_warning(fit <- mmes(y~trait+x,
    random=~vsm(fam(trait,1,shared="axes"),ism(id)),
    rcov=~vsm(fam(trait,1,shared="axes"),ism(record)),
    data=fixture$data,control=control),"did not converge")
  expect_identical(fit$conditionalFA$optimizer,"ai")
  expect_true(any(fit$conditionalFA$updates=="ai"))
  expect_true(all(diff(as.numeric(fit$llik))>=-1e-8))
  expect_true(all(is.finite(fitted(fit))))
  expect_true(fit$executionPlan$admission$conditionalFA$active)
  expect_identical(fit$executionPlan$representation$equations,"conditional")
  expect_equal(fit$conditionalFA$state$weightG+fit$conditionalFA$state$weightE,
    rep(1,1))
  expect_error(mmes_control(update=list(method="em")),"requires representation")
  expect_error(mmes_control(representation=list(equations="conditional"),
    proposal=list(extrapolation="aitken")),"Incompatible")
})

test_that("conditional AI enforces numerical bounds with projected-score convergence",{
  fixture <- .conditional_ai_fixture()
  data <- fixture$data
  data$y <- as.integer(data$trait)+.2*data$x
  sg <- fam(data$trait,1,fixed=rep(TRUE,7L))
  fit <- mmes(y~trait+x,
    random=~vsm(sg,ism(id),sigma2=1,fixedSigma2=TRUE),
    rcov=~vsm(ism(units),sigma2=.5),data=data,
    control=mmes_control(representation=list(equations="conditional"),
      convergence=list(maxIters=80L,parameterTol=1e-8),
      output=list(covariance="none",coefficientMatrix=FALSE),
      execution=list(inverseTol=1e-4),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_true(fit$convergence)
  expect_equal(unname(fit$conditionalFA$state$specificE),
    rep(1e-4*var(data$y),4),tolerance=1e-10)
  expect_equal(fit$conditionalFA$ai$scoreNorm,0)
  expect_length(fit$conditionalFA$ai$activeBounds,1L)
  expect_true(all(diff(as.numeric(fit$llik))>=-1e-8))
  model <- fixture$model; engine <- fixture$engine
  state <- model$state
  state$specificE[1L] <- engine$conditionalInverseTol*model$responseVariance/10
  posterior <- sommer:::.conditional_fa_posterior(model,state,NULL,engine,TRUE)
  expect_error(sommer:::.conditional_fa_ai_step(model,state,posterior,engine,1e-4),
    "starting variances/allocations")
})

test_that("explicit conditional EM preserves the earlier selector exactly",{
  fixture <- .conditional_ai_fixture("axes")
  run <- function(control) suppressWarnings(mmes(y~trait+x,
    random=~vsm(fam(trait,1,shared="axes"),ism(id)),
    rcov=~vsm(fam(trait,1,shared="axes"),ism(record)),
    data=fixture$data,control=control))
  modules <- list(convergence=list(maxIters=2L),
    output=list(covariance="none",coefficientMatrix=FALSE),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  old <- run(do.call(mmes_control,c(list(update=list(method="conditional-fa")),modules)))
  new <- run(do.call(mmes_control,c(list(representation=list(equations="conditional"),
    update=list(method="em")),modules)))
  expect_identical(old$llik,new$llik)
  expect_identical(old$covPar,new$covPar)
  expect_identical(old$bu,new$bu)
  expect_identical(new$conditionalFA$optimizer,"em")
})

test_that("conditional AI records EM recovery when proposals are unavailable",{
  fixture <- .conditional_ai_fixture()
  local_mocked_bindings(.conditional_fa_ai_step=function(model,state,posterior,
    engine,parameterTol) list(accepted=FALSE,stationary=FALSE,scoreNorm=1,
      rank=0L,evaluations=0L,workspaceBytes=0,activeBounds=integer(),
      varianceFloor=1e-6,allocationFloor=1e-6),.package="sommer")
  expect_warning(fit <- mmes(y~trait+x,
    random=~vsm(fam(trait,1),ism(id)),
    rcov=~vsm(dsm(trait),ism(record)),data=fixture$data,
    control=mmes_control(representation=list(equations="conditional"),
      convergence=list(maxIters=2L),
      output=list(covariance="none",coefficientMatrix=FALSE),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE))),"did not converge")
  expect_identical(fit$conditionalFA$updates,rep("em-recovery",2L))
  expect_true(all(diff(as.numeric(fit$llik))>=-1e-8))
})
