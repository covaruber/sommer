nonlinearJointControl <- function(curve="exponential",...){
  start <- switch(curve,exponential=c(asymptote=2,rate=.6),
    saturation=c(asymptote=2,half=.6),
    logistic=c(asymptote=2,rate=1.1,midpoint=2.5))
  sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(start=start,
    population=list(order=8L,maxOrder=64L,tol=1e-6))),
    convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE),...)
}

nonlinearJointData <- function(curve="exponential",groups=10L,points=12L,seed=912L){
  set.seed(seed)
  data <- expand.grid(x=seq(0,5,length.out=points),group=factor(seq_len(groups)))
  u <- matrix(rnorm(2L*groups,sd=.08),groups,2L)
  u[,2L] <- .4*u[,1L]+u[,2L]
  theta <- nonlinearJointControl(curve)$task$nonlinear$start
  theta[names(theta)!="midpoint"] <- log(theta[names(theta)!="midpoint"])
  deviation <- matrix(0,nrow(data),length(theta))
  deviation[,1:2] <- u[data$group,]
  data$eta <- as.numeric(sommer:::sommer_nonlinear_curve_joint_eval(curve,data$x,
    theta,deviation))
  data$y <- data$eta+rnorm(nrow(data),sd=.1)
  data
}

nonlinearJointModel <- function(curve="exponential",structure="usm"){
  formula <- switch(structure,
    usm=~vsm(usm(.parameter,theta=matrix(c(1,.3,.3,1),2)),
      ism(group),sigma2=.01),
    fixed=~vsm(usm(.parameter,theta=matrix(c(1,.3,.3,1),2),fixed=rep(TRUE,2)),
      ism(group),sigma2=.01,fixedSigma2=TRUE),
    rr=~vsm(rrm(.parameter,k=1),ism(group),sigma2=.01))
  nonlinear_curve("x",curve,joint=list(parameters=c("asymptote",
    if(curve=="saturation") "half" else "rate"),random=formula))
}

nonlinearJointFixture <- function(curve="exponential",family=gaussian(),structure="usm"){
  data <- nonlinearJointData(curve,groups=3L,points=8L)
  data$z <- rep(c(-1,1),length.out=nrow(data))
  mixture <- if(inherits(family,"sommer_familym")) family else NULL
  if(!is.null(mixture)) family <- gaussian()
  if(family$family=="poisson") data$y <- round(exp(data$eta))
  if(family$family=="binomial") data$y <- rep(c(0,1,1),length.out=nrow(data))
  if(family$family=="Gamma") data$y <- exp(data$eta)
  if(grepl("^Negative Binomial",family$family)) data$y <- round(exp(data$eta))
  if(family$family=="ziPoisson") data$y <- rep(c(0,2,3),length.out=nrow(data))
  if(family$family=="ordinal")
    data$y <- ordered(rep(c("low","middle","high"),length.out=nrow(data)),
      levels=c("low","middle","high"))
  if(family$family=="multinomial")
    data$y <- factor(rep(c("red","green","blue"),length.out=nrow(data)),
      levels=c("red","green","blue"))
  if(!is.null(mixture)){
    data$trait <- rep(c("count","normal"),each=4L,length.out=nrow(data))
    selected <- data$trait=="count"
    data$y[selected] <- if(mixture$families$count$family=="ordinal")
      rep(1:3,length.out=sum(selected)) else round(exp(data$eta[selected]))
    family <- mixture
  }
  model <- nonlinearJointModel(curve,structure)
  control <- nonlinearJointControl(curve)
  random <- sommer:::.nonlinear_joint_formula(NULL,model,data)
  setup <- sommer:::.nonlinear_family_prepare(y~0,random,data,family,"x",control,
    "fail","fail","fail","fail",NULL)
  prepared <- setup$prepared
  layout <- sommer:::.nonlinear_random_layout(model,prepared)
  structures <- prepared$covStruct[seq_along(prepared$rtermss)]
  z <- do.call(cbind,prepared$Z)
  context <- sommer:::sommer_nonlinear_family_create(prepared$X,z,
    setup$response$nativeResponse,setup$offset,setup$input,curve,
    layout$coordinate[1L],layout$start,structures,prepared$Ai,
    setup$specifications,setup$index-1L,setup$response$trials,
    coordinateMap=layout$coordinateMap*1.0)
  theta <- sommer:::.nonlinear_coordinates(model,control$task$nonlinear$start)
  parameters <- c(theta,unlist(lapply(structures,function(s) s$par[s$free])),
    unlist(lapply(setup$specifications,function(s){
      if(!s$free) return(numeric())
      if(s$family=="ordinal") sommer:::.laplace_ordinal_coordinates(s$thresholds)
      else if(is.null(s$zeroDesign)) s$logParameter else s$zeroCoefficients
    })))
  list(data=data,model=model,setup=setup,layout=layout,structures=structures,
    z=z,parameters=parameters,context=context,
    evaluate=function(at,gradient=TRUE)
      sommer:::sommer_nonlinear_family_eval(context,at,gradient))
}

test_that("joint syntax preserves old formulas and rejects ambiguous parameter axes",{
  old <- nonlinear_curve("x",random=list(rate=~vsm(ism(group))))
  expect_identical(old$version,2L)
  expect_identical(sommer:::.nonlinear_validate(old),old)
  independent <- nonlinear_curve("x",random=list(asymptote=~vsm(ism(group)),
    rate=~vsm(ism(group))))
  expect_identical(independent$version,3L)
  expect_identical(sommer:::.nonlinear_validate(independent),independent)
  joint <- nonlinearJointModel()
  expect_identical(sommer:::.nonlinear_validate(joint),joint)
  expect_error(nonlinear_curve("x",random=list(rate=~group),
    joint=list(parameters=c("asymptote","rate"),random=~vsm(usm(.parameter),ism(group)))),
    "disjoint")
  expect_error(nonlinear_curve("x",joint=list(parameters=c("rate","rate"),
    random=~group)),"joint")
  expect_error(nonlinear_curve("x",random=list(rate=~group,rate=~group)),"uniquely")
})

for(curve in c("exponential","logistic","saturation")){
  test_that(paste(curve,"joint conditional curvature matches independent dense algebra"),{
    fixture <- nonlinearJointFixture(curve)
    at <- fixture$parameters
    result <- fixture$evaluate(at)
    expect_true(result$valid,info=result$reason)
    conditional <- sommer:::sommer_nonlinear_family_conditional(
      fixture$context,at,result$latentMode)
    k <- length(fixture$model$parameters)
    descriptor <- fixture$structures[[1L]]
    descriptor$par[descriptor$free] <- at[k+seq_len(sum(descriptor$free))]
    covariance <- kronecker(sommer:::evaluate_covstruct_cpp(descriptor,descriptor$par),
      solve(as.matrix(fixture$setup$prepared$Ai[[1L]])))
    precision <- solve(covariance)
    bases <- lapply(seq_len(k),function(j)
      as.matrix(fixture$z)%*%diag(fixture$layout$coordinateMap[,j]))
    dense <- function(mode){
      deviations <- do.call(cbind,lapply(bases,function(b) as.numeric(b%*%mode)))
      mean <- as.numeric(sommer:::sommer_nonlinear_curve_joint_eval(curve,
        fixture$data$x,at[seq_len(k)],deviations))
      variance <- unname(exp(tail(at,1)))
      .5*(nrow(fixture$data)*log(2*pi*variance)+sum((fixture$data$y-mean)^2)/variance+
        as.numeric(crossprod(mode,precision%*%mode)))
    }
    expect_equal(conditional$objective,dense(result$latentMode),tolerance=1e-9)
    expected <- optimHess(as.numeric(result$latentMode),dense,
      control=list(ndeps=rep(1e-4,length(result$latentMode))))
    expect_equal(as.matrix(conditional$hessian),expected,tolerance=2e-5,
      ignore_attr=TRUE)
    expectedLogLik <- -dense(result$latentMode)-.5*as.numeric(determinant(covariance,
      logarithm=TRUE)$modulus)-.5*as.numeric(determinant(expected,logarithm=TRUE)$modulus)
    expect_equal(result$logLik,expectedLogLik,tolerance=2e-5)
    numerical <- vapply(seq_along(at),function(j){
      up <- down <- at; up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
      (fixture$evaluate(up,FALSE)$logLik-fixture$evaluate(down,FALSE)$logLik)/2e-5
    },numeric(1))
    expect_equal(as.numeric(result$gradient),numerical,tolerance=5e-5)
  })
}

for(family in list(binomial(),binomial("probit"),binomial("cloglog"),poisson(),
    nbinom2(3),Gamma("log"),zipoisson(.25),ordinal(),ordinal("probit"),
    multinomial(reference="red"),zipoisson(zero=~z))){
  test_that(paste(family$family,family$link,"joint Laplace scores include cross curvature"),{
    fixture <- nonlinearJointFixture(family=family)
    at <- fixture$parameters
    result <- fixture$evaluate(at)
    expect_true(result$valid,info=result$reason)
    numerical <- vapply(seq_along(at),function(j){
      up <- down <- at; up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
      a <- fixture$evaluate(up,FALSE); b <- fixture$evaluate(down,FALSE)
      expect_true(a$valid,info=a$reason); expect_true(b$valid,info=b$reason)
      (a$logLik-b$logLik)/2e-5
    },numeric(1))
    expect_equal(as.numeric(result$gradient),numerical,tolerance=8e-5)
  })
}

for(categorical in c(FALSE,TRUE)){
  test_that(paste("joint mixed",if(categorical) "categorical" else "scalar",
      "scores preserve family coordinates"),{
    family <- familym("trait",count=if(categorical) ordinal() else poisson(),
      normal=gaussian())
    fixture <- nonlinearJointFixture(family=family)
    at <- fixture$parameters
    result <- fixture$evaluate(at)
    expect_true(result$valid,info=result$reason)
    numerical <- vapply(seq_along(at),function(j){
      up <- down <- at; up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
      a <- fixture$evaluate(up,FALSE); b <- fixture$evaluate(down,FALSE)
      expect_true(a$valid,info=a$reason); expect_true(b$valid,info=b$reason)
      (a$logLik-b$logLik)/2e-5
    },numeric(1))
    expect_equal(as.numeric(result$gradient),numerical,tolerance=8e-5)
  })
}

test_that("joint reduced-rank scores include moving parameter-axis loadings",{
  fixture <- nonlinearJointFixture(structure="rr")
  at <- fixture$parameters
  result <- fixture$evaluate(at)
  expect_true(result$valid,info=result$reason)
  expect_lt(length(result$latentMode),length(result$mode))
  numerical <- vapply(seq_along(at),function(j){
    up <- down <- at; up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
    (fixture$evaluate(up,FALSE)$logLik-fixture$evaluate(down,FALSE)$logLik)/2e-5
  },numeric(1))
  expect_equal(as.numeric(result$gradient),numerical,tolerance=8e-5)
})

test_that("joint fits predict simulate serialize and defer uncertainty",{
  skip_on_cran()
  data <- nonlinearJointData()
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=nonlinearJointModel(structure="fixed"),
    control=nonlinearJointControl())
  expect_true(fit$convergence)
  expect_null(fit$nonlinearInference)
  expect_equal(unname(fit$nonlinearParameters),c(2,.6),tolerance=.08)
  expect_equal(predict(fit),as.numeric(fitted(fit)),ignore_attr=TRUE)
  expect_equal(predict(fit,data[1:3,]),as.numeric(fitted(fit)[1:3]),tolerance=1e-10)
  expect_equal(dim(fit$groupParameters),c(10L,2L))
  expect_gt(max(abs(predict(fit,data[1:3,],re.form="population")-
    predict(fit,data[1:3,],re.form=NA))),1e-5)
  expect_identical(simulate(fit,2,seed=16,re.form="population"),
    simulate(fit,2,seed=16,re.form="population"))
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(predict(restored,data[1:3,]),predict(fit,data[1:3,]))
  expect_true(all(is.finite(profile_nonlinear(restored,replications=1L)$table$seconds)))
  post <- post_nonlinear(restored,mode=2)
  expect_true(all(is.finite(vcov(post,type="all"))))
  expect_equal(dim(post$nonlinearInference$groupParameterSE),dim(fit$groupParameters))
  expect_true(all(post$nonlinearInference$groupParameterSE>0))
  expect_equal(diag(vcov(post,type="group")),
    as.numeric(post$nonlinearInference$groupParameterSE)^2,ignore_attr=TRUE)
  expect_true(isSymmetric(vcov(post,type="group")))
  natural <- as.numeric(fit$groupParameters)
  randomJacobian <- matrix(0,length(natural),length(fit$randomMode))
  randomJacobian[cbind(seq_along(natural),
    as.numeric(fit$nonlinearRandom$groupMap))] <- natural
  direct <- matrix(0,length(natural),length(fit$optimizer$par))
  direct[cbind(seq_along(natural),rep(1:2,each=nrow(fit$groupParameters)))] <- natural
  total <- direct+randomJacobian%*%post$nonlinearInference$modeSensitivity
  expected <- total%*%vcov(post,type="working")%*%t(total)+
    randomJacobian%*%vcov(post,type="conditional")%*%t(randomJacobian)
  expect_equal(vcov(post,type="group"),expected,tolerance=1e-6,ignore_attr=TRUE)
  prediction <- predict(post,data[1:2,],se.fit=TRUE,cov.fit=TRUE)
  expect_equal(prediction$se.fit^2,diag(prediction$cov.fit),ignore_attr=TRUE)
  expect_null(fit$nonlinearInference)
})

test_that("independent named curves and ordinary additive vsm terms coexist",{
  data <- nonlinearJointData(groups=6L)
  model <- nonlinear_curve("x",random=list(
    asymptote=~vsm(ism(group),sigma2=.01,fixedSigma2=TRUE),
    rate=~vsm(ism(group),sigma2=.01,fixedSigma2=TRUE)))
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.002,fixedSigma2=TRUE),
    data=data,nonlinear=model,control=nonlinearJointControl())
  expect_true(fit$convergence)
  expect_equal(predict(fit,data[1:3,]),as.numeric(fitted(fit)[1:3]),tolerance=1e-10)
  expect_equal(length(fit$additiveRandomMode),6L)
  expect_equal(dim(fit$groupParameters),c(6L,2L))
})

test_that("all logistic coordinates vary jointly on their declared scales",{
  data <- nonlinearJointData("logistic",groups=6L)
  model <- nonlinear_curve("x","logistic",joint=list(
    parameters=c("asymptote","rate","midpoint"),
    random=~vsm(usm(.parameter,theta=diag(3),fixed=rep(TRUE,5)),
      ism(group),sigma2=.005,fixedSigma2=TRUE)))
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=nonlinearJointControl("logistic"))
  expect_true(fit$convergence)
  expect_equal(dim(fit$groupParameters),c(6L,3L))
  expect_true(all(fit$groupParameters[,1:2]>0))
  expect_equal(predict(fit,data[1:2,]),as.numeric(fitted(fit)[1:2]),tolerance=1e-10)
  expect_true(all(is.finite(predict(fit,data[1:2,],re.form="population"))))
  expect_null(fit$nonlinearInference)
})

test_that("joint covariance and independent midpoint coexist without changing vsm",{
  skip_on_cran()
  data <- nonlinearJointData("logistic",groups=6L)
  model <- nonlinearJointModel("logistic","fixed")
  model <- nonlinear_curve("x","logistic",random=list(
    midpoint=~vsm(ism(group),sigma2=.002,fixedSigma2=TRUE)),joint=model$joint)
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=nonlinearJointControl("logistic"))
  expect_true(fit$convergence)
  expect_equal(dim(fit$groupParameters),c(6L,3L))
  expect_equal(predict(fit,data[1:2,]),as.numeric(fitted(fit)[1:2]),tolerance=1e-10)
  expect_equal(length(fit$randomMode),18L)
})

test_that("joint known relationships preserve full covariance and prediction basis",{
  data <- nonlinearJointData(groups=6L)
  precision <- Matrix::Matrix(diag(6L)+matrix(.1,6L,6L),sparse=TRUE)
  dimnames(precision) <- list(levels(data$group),levels(data$group))
  attr(precision,"inverse") <- TRUE
  model <- nonlinear_curve("x",joint=list(parameters=c("asymptote","rate"),
    random=~vsm(usm(.parameter,theta=diag(2),fixed=rep(TRUE,2)),
      ism(group),Gu=precision,sigma2=.01,fixedSigma2=TRUE)))
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=nonlinearJointControl())
  expect_true(fit$convergence)
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(predict(restored,data[1:2,]),as.numeric(fitted(fit)[1:2]),tolerance=1e-10)
  expect_error(predict(fit,transform(data[1:2,],group="unseen")),"unseen")
})

for(family in list(poisson(),ordinal(),multinomial(reference="red"))){
  test_that(paste(family$family,"joint fits retain response and uncertainty shapes"),{
    skip_on_cran()
    data <- nonlinearJointData(groups=8L,points=14L)
    set.seed(351)
    if(family$family=="poisson") data$y <- rpois(nrow(data),exp(data$eta))
    if(family$family=="ordinal"){
      p <- cbind(plogis(-.6-data$eta),
        plogis(.8-data$eta)-plogis(-.6-data$eta),1-plogis(.8-data$eta))
      code <- 1L+rowSums(runif(nrow(data))>t(apply(p,1,cumsum))[,-3,drop=FALSE])
      data$y <- ordered(c("low","middle","high")[code],
        levels=c("low","middle","high"))
    }
    if(family$family=="multinomial"){
      p <- cbind(1,exp(data$eta),exp(data$eta))
      p <- p/rowSums(p)
      code <- 1L+rowSums(runif(nrow(data))>t(apply(p,1,cumsum))[,-3,drop=FALSE])
      data$y <- factor(c("red","green","blue")[code],levels=c("red","green","blue"))
    }
    fit <- sommer:::.mmes_dispatch(y~0,data=data,family=family,
      nonlinear=nonlinearJointModel(structure="fixed"),control=nonlinearJointControl())
    expect_true(fit$convergence)
    expect_null(fit$nonlinearInference)
    prediction <- predict(fit,data[1:2,])
    expected <- fitted(fit)
    expected <- if(is.matrix(expected)) expected[1:2,,drop=FALSE] else expected[1:2]
    expect_equal(prediction,expected,ignore_attr=TRUE)
    expect_true(all(is.finite(unlist(predict(fit,data[1:2,],re.form="population")))))
    post <- post_nonlinear(fit,mode=2)
    uncertainty <- predict(post,data[1:2,],se.fit=TRUE,cov.fit=TRUE)
    expect_equal(as.numeric(uncertainty$se.fit)^2,diag(uncertainty$cov.fit),
      ignore_attr=TRUE)
    expect_equal(dim(uncertainty$fit),dim(prediction))
    expect_true(all(is.finite(vcov(post,type="group"))))
    expect_identical(simulate(post,1,seed=12),simulate(post,1,seed=12))
  })
}

test_that("joint reduced-rank fits restore post-fit full group covariance",{
  data <- nonlinearJointData(groups=6L)
  model <- nonlinear_curve("x",joint=list(parameters=c("asymptote","rate"),
    random=~vsm(rrm(.parameter,k=1,loadings=matrix(c(1,.5),2,1),
      fixed=c(TRUE,TRUE)),ism(group),sigma2=.01,fixedSigma2=TRUE)))
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=nonlinearJointControl())
  expect_true(fit$convergence)
  expect_lt(length(fit$latentRandomMode),length(fit$randomMode))
  post <- post_nonlinear(fit,mode=2)
  expect_equal(dim(vcov(post,type="group")),c(12L,12L))
  expect_equal(diag(vcov(post,type="group")),
    as.numeric(post$nonlinearInference$groupParameterSE)^2,ignore_attr=TRUE)
  expect_equal(predict(post,data[1:2,]),as.numeric(fitted(fit)[1:2]),tolerance=1e-10)
})

test_that("joint specifications reject repeated and missing parameter axes explicitly",{
  data <- nonlinearJointData(groups=3L)
  expect_error(nonlinear_curve("x",joint=list(parameters=c("asymptote","rate"),
    random=~vsm(usm(.parameter),usm(.parameter),ism(group)))),
    "exactly one")
  expect_error(nonlinear_curve("x",random=list(asymptote=~group+other,rate=~group)),
    "one-term")
  expect_error(nonlinear_curve("x",joint=list(parameters=c("asymptote","rate"),
    random=~vsm(ism(group)))),"containing .parameter")
})

test_that("covariance constraints retain independent parameter formula semantics",{
  data <- nonlinearJointData(groups=8L)
  model <- nonlinear_curve("x",random=list(
    asymptote=~vsm(ism(group),sigma2=.01),
    rate=~vsm(ism(group),sigma2=.01)))
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=nonlinearJointControl(),
    vcc=data.frame(parameter=1:2,group=1L,scale=c(1,1)))
  expect_true(fit$convergence)
  expect_equal(length(coef(fit,type="working")),4L)
  expect_false(anyDuplicated(names(coef(fit,type="all")))>0L)
  expect_equal(fit$covPar[[1L]],fit$covPar[[2L]],ignore_attr=TRUE)
  expect_equal(predict(fit,data[1:2,]),as.numeric(fitted(fit)[1:2]),tolerance=1e-10)
})

test_that("diagonal joint covariance agrees with independent named formulas",{
  data <- nonlinearJointData(groups=6L)
  independent <- nonlinear_curve("x",random=list(
    asymptote=~vsm(ism(group),sigma2=.01,fixedSigma2=TRUE),
    rate=~vsm(ism(group),sigma2=.01,fixedSigma2=TRUE)))
  joint <- nonlinear_curve("x",joint=list(parameters=c("rate","asymptote"),
    random=~vsm(usm(.parameter,theta=diag(2),fixed=rep(TRUE,2)),
      ism(group),sigma2=.01,fixedSigma2=TRUE)))
  a <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=independent,control=nonlinearJointControl())
  b <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=joint,control=nonlinearJointControl())
  expect_true(a$convergence); expect_true(b$convergence)
  expect_equal(a$llik,b$llik,tolerance=1e-9)
  expect_equal(coef(a),coef(b),tolerance=1e-8)
  expect_equal(fitted(a),fitted(b),tolerance=1e-8)
  expect_equal(a$groupParameters,b$groupParameters,tolerance=1e-8)
  expect_false(anyDuplicated(names(coef(a,type="all")))>0L)
  expect_equal(predict(a,data[1:2,],re.form="population"),
    predict(b,data[1:2,],re.form="population"),tolerance=1e-7)
})

test_that("joint mixed-family fits retain scalar trait output and post-fit covariance",{
  skip_on_cran()
  data <- nonlinearJointData(groups=6L,points=10L)
  normal <- data
  data$trait <- "count"; normal$trait <- "normal"
  set.seed(351)
  data$y <- rpois(nrow(data),exp(data$eta))
  data <- rbind(data,normal)
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=nonlinearJointModel(structure="fixed"),
    family=familym("trait",count=poisson(),normal=gaussian()),
    control=nonlinearJointControl())
  expect_true(fit$convergence)
  selected <- c(1L,61L)
  expect_equal(predict(fit,data[selected,]),fitted(fit)[selected],ignore_attr=TRUE)
  expect_true(all(is.finite(predict(fit,data[selected,],re.form="population"))))
  post <- post_nonlinear(fit,mode=2)
  expect_true(all(is.finite(vcov(post,type="group"))))
  uncertainty <- predict(post,data[selected,],se.fit=TRUE,cov.fit=TRUE)
  expect_equal(uncertainty$se.fit^2,diag(uncertainty$cov.fit),ignore_attr=TRUE)
  expect_null(fit$nonlinearInference)
})
