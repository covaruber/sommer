specificationProfileControl <- function(){
  sommer:::.mmes_control_build(task=list(objective="ml"),
    convergence=list(maxIters=120L,likelihoodTol=1e-9,parameterTol=1e-9),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

specificationProfileData <- function(){
  set.seed(5291)
  data <- expand.grid(x=seq(-1.5,1.5,length.out=8L),
    group=factor(seq_len(4L)))
  data$z <- rep(c(-1,1),length.out=nrow(data))
  data$y <- .35+.75*data$x+.3*data$z+
    rep(c(-.16,.04,.12,0),each=8L)+rnorm(nrow(data),sd=.18)
  data
}

test_that("specification performance profile times saved operations without refitting",{
  data <- specificationProfileData()
  control <- specificationProfileControl()
  fit <- sommer:::.mmes_specification_fit(y~a*x,data=data,
    random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
    familyParameters=list(dispersion=.0324,estimate=FALSE),
    parameters=list(a=parameter()),start=c(a=.6),control=control)
  parameters <- fit$optimizer$par
  fitted <- fit$fitted.values
  profile <- profile_nonlinear(fit,replications=1L)
  expect_named(profile$table,
    c("phase","seconds","minimumSeconds","replications"))
  expect_setequal(profile$table$phase,
    c("predictor","jacobian","conditionalCurvature","objective",
      "score","prediction"))
  expect_true(all(is.finite(profile$table$seconds)))
  expect_true(all(profile$table$seconds>=profile$table$minimumSeconds))
  expect_true(all(profile$table$replications==1L))
  expect_equal(profile$contextSeconds,0)
  expect_equal(profile$likelihood,fit$likelihoodMethod)
  expect_equal(fit$optimizer$par,parameters)
  expect_equal(fit$fitted.values,fitted)
  expect_match(profile$note,"outer model is not refit")
})

test_that("specification ANOVA allows only regular nested fixed-term tests",{
  data <- specificationProfileData()
  data$y <- .35+.75*data$x+.3*data$z+rnorm(nrow(data),sd=.18)
  control <- specificationProfileControl()
  familyParameters <- list(dispersion=.0324,estimate=FALSE)
  parameters <- list(a=parameter())
  small <- sommer:::.mmes_specification_fit(y~a*x,data=data,
    familyParameters=familyParameters,parameters=parameters,start=c(a=.6),
    control=control)
  large <- sommer:::.mmes_specification_fit(y~z+a*x,data=data,
    familyParameters=familyParameters,parameters=parameters,start=c(a=.6),
    control=control)
  comparison <- stats::anova(small,large)
  expect_s3_class(comparison,"anova")
  expect_equal(comparison["larger","Df"],1)
  expect_true(is.finite(comparison["larger","Pr(>Chisq)"]))
  parm <- names(coef(small))
  null <- stats::setNames(numeric(length(parm)),parm)
  wald <- stats::anova(small,parm=parm,null=null)
  expect_s3_class(wald,"anova")
  expect_true(is.finite(wald$Chisq))
  expect_error(stats::anova(large,small),"strictly nested|supply the smaller")
  expect_error(stats::anova(small,parm=parm,null=unname(null)),"exactly named null")
  expect_error(stats::anova(small,parm="family:response:dispersion.1",
    null=c("family:response:dispersion.1"=0)),"mean-coordinate")
  changed <- data
  changed$x[1L] <- changed$x[1L]+.1
  other <- sommer:::.mmes_specification_fit(y~a*x,data=changed,
    familyParameters=familyParameters,parameters=parameters,start=c(a=.6),
    control=control)
  expect_error(stats::anova(small,other),"identical data")
})
