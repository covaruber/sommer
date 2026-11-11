test_that("method defaults preserve REML and PQL working fits", {
  method <- mmes_method()
  control <- mmes_control()
  resolved <- sommer:::.mmes_resolve_controls(control,method,FALSE,NULL,list())
  expect_identical(resolved$task$objective,"reml")
  expect_identical(resolved$task$glmm$method,"pql")
  expect_identical(resolved$task$pql,sommer:::.mme_pql_defaults)
  expect_identical(resolved$task$glmm,sommer:::.mme_glmm_defaults)
  expect_identical(control$solve,sommer:::.mmes_control_internal()$solve)
  expect_error(mmes_control(task=list(objective="ml")),"unused argument")
  expect_error(mmes_method(integration="unknown"),"integration")
  expect_error(mmes_control(optimizer=list(unknown=1)),"Unknown optimizer")
  expect_error(mmes_control(mode=list(tol=0)),"mode\\$tol")
  expect_null(control$task)
  expect_identical(control$optimizer$gradient,"auto")
  control$optimizer$maxit <- 73L
  resolved <- sommer:::.mmes_resolve_controls(control,method,FALSE,NULL,list())
  expect_identical(resolved$convergence$maxIters,73L)
  control$optimizer$unknown <- TRUE
  expect_error(sommer:::.mmes_resolve_controls(control,method,FALSE,NULL,list()),
    "Unknown optimizer")
})

test_that("outer and conditional mode controls have distinct ownership", {
  control <- mmes_control(optimizer=list(maxit=123L),
    mode=list(maxit=17L,tol=1e-7))
  method <- mmes_method(objective="ml",integration="laplace")
  resolved <- sommer:::.mmes_resolve_controls(control,method,FALSE,NULL,list())
  expect_identical(resolved$task$glmm$maxit,123L)
  expect_identical(resolved$task$glmm$modeMaxit,17L)
  expect_equal(resolved$task$glmm$modeTol,1e-7)
  resolved <- sommer:::.mmes_resolve_controls(control,method,TRUE,
    list(parameters=c(asymptote=3,rate=1)),list())
  expect_identical(resolved$convergence$maxIters,123L)
  expect_identical(resolved$task$nonlinear$mode$maxIters,17L)
  expect_equal(resolved$task$nonlinear$start,c(asymptote=3,rate=1))
  expect_identical(resolved$task$glmm$method,"pql")
  expect_error(sommer:::.mmes_resolve_controls(mmes_control(),
    mmes_method(approximation="pql",integration="laplace"),FALSE,NULL,list()),
    "working-model")
})

test_that("Gaussian mixed-model defaults are unchanged at the public boundary", {
  set.seed(91)
  dat <- data.frame(group=factor(rep(1:8,each=5)),x=rep(1:5,8))
  dat$y <- 2+0.3*dat$x+rep(rnorm(8),each=5)+rnorm(40,sd=0.3)
  ctrl <- mmes_control(diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  actual <- mmes(y~x,random=~group,data=dat,control=ctrl)
  reference <- sommer:::.mmes_dispatch(y~x,random=~group,data=dat,
    control=sommer:::.mmes_control_internal(
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_equal(actual$b,reference$b,tolerance=1e-10)
  expect_equal(actual$llik,reference$llik,tolerance=1e-10)
  expect_identical(actual$REML,TRUE)
  expect_identical(actual$method$objective,"reml")
  expect_null(actual$control$task)
})

test_that("default binomial mixed fits retain the REML PQL numerical path", {
  set.seed(94)
  dat <- data.frame(group=factor(rep(1:10,each=12)),x=rep(seq(-1,1,length.out=12),10))
  dat$y <- rbinom(nrow(dat),1,plogis(.2+.7*dat$x+
    rep(rnorm(10,sd=.7),each=12)))
  actual <- mmes(y~x,random=~group,data=dat,family=binomial(),
    control=mmes_control(diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  reference <- sommer:::.mmes_dispatch(y~x,random=~group,data=dat,
    family=binomial(),control=sommer:::.mmes_control_internal(
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_equal(actual$b,reference$b,tolerance=1e-10)
  expect_equal(fitted(actual),fitted(reference),tolerance=1e-10)
  expect_true(actual$REML)
  expect_identical(actual$method$objective,"reml")
})

test_that("native curves use named predictor components and natural starts", {
  set.seed(92)
  dat <- data.frame(group=factor(rep(1:6,each=8)),x=rep(seq(0,4,length.out=8),6))
  dat$y <- 4*exp(-0.5*dat$x)+rep(rnorm(6,sd=0.3),each=8)+rnorm(48,sd=0.1)
  actual <- mmes(y~0+graph_term("growth"),random=~group,data=dat,
    predictors=list(growth=nonlinear_curve("x")),
    method=mmes_method(objective="ml",integration="exact"),
    start=list(parameters=c(asymptote=4,rate=0.5)),
    control=mmes_control(optimizer=list(maxit=200L),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_s3_class(actual,"mmes.nonlinear")
  expect_true(actual$convergence)
  expect_equal(unname(actual$nonlinearParameters),c(4,0.5),tolerance=0.3)
  expect_error(mmes(y~x,random=~group,data=dat,
    predictors=list(growth=nonlinear_curve("x")),
    method=mmes_method(objective="ml")),"graph_term")
  expect_error(mmes(y~x,random=~group,data=dat,nonlinear=nonlinear_curve("x")),
    "unused argument")
  expect_equal(predict(actual,target="conditional"),fitted(actual),tolerance=1e-9)
  expect_equal(predict(actual,target="zero"),predict(actual,re.form=NA))
  expect_error(predict(actual,target="zero",re.form=NA),"not both")
  post <- post_nonlinear(actual)
  expect_true(all(is.finite(vcov(post))))
})

test_that("fixed-only and parameter-regression calls use the public specification", {
  set.seed(93)
  dat <- expand.grid(x=seq(-2,2,length.out=9),Env=factor(c("A","B")))
  dat$y <- 1+.2*(dat$Env=="B")+2*tanh(.8*dat$x)+rnorm(nrow(dat),sd=.1)
  ml <- mmes_method(objective="ml",integration="exact")
  ctrl <- mmes_control(optimizer=list(maxit=300L),
    convergence=list(likelihoodTol=1e-10,parameterTol=1e-9),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  fit <- mmes(y~Env+a*tanh(k*x),parameters=list(
    a=parameter(),k=parameter(link="log")),data=dat,
    start=list(parameters=c(a=2,k=.8)),method=ml,control=ctrl,rcov=~units)
  expect_s3_class(fit,"mmes.specification")
  expect_true(fit$convergence)
  expect_equal(predict(fit,dat[1:3,]),fitted(fit)[1:3],tolerance=1e-9)
  expect_null(fit$control$task)
  expect_identical(fit$control$optimizer$maxit,300L)
  expect_match(fit$executionPlan$derivatives,"numerical")
  expect_true(all(is.finite(vcov(postPEV(fit)))))
  expect_error(anova(fit),"Supply object2")
  dat$h <- as.integer(dat$x>0)
  fixed <- mmes(y~Env+x,data=dat,method=ml)
  expect_equal(as.numeric(fixed$b),unname(coef(lm(y~Env+x,dat))),tolerance=1e-8)
  expect_error(mmes(h~x,data=dat,family=binomial()),"require.*ml")
})

test_that("fixed-only Gaussian models retain rank handling and covariance output", {
  set.seed(96)
  dat <- expand.grid(locale=factor(c("north","south")),
    cultivar=factor(c("early","late")),replicate=1:6)
  dat <- subset(dat,!(locale=="south" & cultivar=="late"))
  dat$y <- 3+as.integer(dat$locale)+as.integer(dat$cultivar)+rnorm(nrow(dat))
  formula <- y~locale+locale:cultivar
  diagnostics <- list(verbose=FALSE,dateWarning=FALSE)
  expect_error(mmes(formula,data=dat,
    control=mmes_control(diagnostics=diagnostics)),"rank deficient")
  expect_warning(fit <- mmes(formula,data=dat,control=mmes_control(
    representation=list(fixedRank="drop"),diagnostics=diagnostics)),
    "rank reduction")
  reference <- lm(formula,data=dat)
  expect_true(fit$convergence)
  expect_false(inherits(fit,"mmes.specification"))
  expect_equal(as.numeric(fitted(fit)),as.numeric(fitted(reference)),tolerance=1e-8)
  expect_equal(as.numeric(fit$covPar[[1L]]),summary(reference)$sigma^2,
    tolerance=1e-5)
  expect_true(nrow(summary(fit)$varcomp)>0L)
  expect_identical(fit$control$representation$fixedRank,"drop")
  expect_true(fit$REML)
})

test_that("correlation imputation retains its fixed-only Gaussian coefficient path", {
  set.seed(97)
  wide <- matrix(rnorm(40)+rep(1:5,each=8),8,5,
    dimnames=list(paste0("g",1:8),paste0("e",1:5)))
  wide[2,3] <- NA
  relationship <- matrix(.2,8,8,dimnames=list(rownames(wide),rownames(wide)))
  diag(relationship) <- 1
  result <- corImputation(wide,Gu=relationship,nearest=3)
  expect_true(all(is.finite(result$imputed)))
  expect_equal(result$imputed[!is.na(wide)],wide[!is.na(wide)])
})

test_that("named native graphs retain natural starts and public controls", {
  set.seed(95)
  dat <- data.frame(group=factor(rep(1:8,each=10)),
    x=rep(seq(-2,2,length.out=10),8))
  dat$y <- 2*tanh(.8*dat$x)+.1+rep(rnorm(8,sd=.2),each=10)+
    rnorm(nrow(dat),sd=.08)
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    scaled=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    hidden=list(op="tanh",input="scaled"),
    output=list(op="affine",inputs="hidden",weights=matrix(3L),bias=4L),
    mean=list(op="identity",input="output")),
    output="mean",parameters=c("rate","hiddenBias","amplitude","offset"),
    parameterMap=list(rate=list(type="fixed",transform="log"),
      hiddenBias=list(type="fixed"),amplitude=list(type="fixed"),
      offset=list(type="fixed")))
  fit <- mmes(y~0+graph_term("net"),random=~group,data=dat,
    predictors=list(net=graph),
    start=list(parameters=c("net::rate"=.8,"net::hiddenBias"=0,
      "net::amplitude"=2,"net::offset"=.1)),
    method=mmes_method(objective="ml",integration="exact"),
    control=mmes_control(optimizer=list(maxit=200L),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_true(fit$convergence)
  expect_null(fit$control$task)
  expect_equal(predict(fit,target="conditional"),fitted(fit),tolerance=1e-8)
  expect_equal(predict(fit,target="zero"),predict(fit,re.form=NA),tolerance=1e-8)
  expect_true(all(is.finite(vcov(post_nonlinear(fit)))))
})
