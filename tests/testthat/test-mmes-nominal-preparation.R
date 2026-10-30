nominalPreparationData <- function(){
  data.frame(y=factor(rep(c("red","green","blue"),12),
      levels=c("red","green","blue")),
    x=seq(-1,1,length.out=36),
    treatment=factor(rep(c("a","b"),each=18)),
    group=factor(rep(letters[1:6],each=6)),
    block=factor(rep(1:3,each=12)), exposure=seq(.6,1.3,length.out=36))
}

nominalPreparationControl <- function(...){
  mmes_control(task=list(objective="ml",glmm=list(method="laplace",modeTol=1e-10)),
    execution=list(blasThreads=1L),diagnostics=list(verbose=FALSE,dateWarning=FALSE),...)
}

prepareNominalFixture <- function(data=nominalPreparationData(),fixed=y~x,
    random=~vsm(usm(.nominalCategory),ism(group),sigma2=.6),
    reference=NULL,control=nominalPreparationControl(),...){
  sommer:::.laplace_prepare_nominal(fixed,random,data,reference,control,...)
}

test_that("nominal preparation labels categories and gives each contrast its own fixed coefficients", {
  data <- nominalPreparationData()
  original <- data
  prepared <- prepareNominalFixture(data,y~x+treatment+offset(log(exposure)),
    reference="green",contrasts=list(treatment="contr.sum"))
  expect_identical(data,original)
  expect_s3_class(prepared,"sommer_nominal_preparation")
  expect_false(inherits(prepared,"mmes"))
  expect_identical(prepared$nominal$levels,c("red","green","blue"))
  expect_identical(prepared$nominal$nativeLevels,c("red","blue","green"))
  expect_identical(prepared$nominal$predictorLevels,c("red","blue"))
  expect_identical(prepared$nominal$reference,"green")
  expect_equal(prepared$nominal$response,rep(c(1,3,2),12))
  expect_equal(prepared$nominal$counts,c(red=12L,blue=12L,green=12L))
  base <- model.matrix(y~x+treatment+offset(log(exposure)),data,
    contrasts.arg=list(treatment="contr.sum"))
  expect_equal(unname(as.matrix(prepared$fixedDesign)),
    unname(kronecker(diag(2),base)))
  expect_identical(colnames(prepared$fixedDesign),
    c("red:(Intercept)","red:x","red:treatment1",
      "blue:(Intercept)","blue:x","blue:treatment1"))
  expect_identical(prepared$nominal$fixed$contrasts,list(treatment="contr.sum"))
  expect_identical(prepared$nominal$fixed$xlevels$treatment,levels(data$treatment))
  expect_equal(prepared$nominal$fixed$coefficientMap$assignment,rep(c(0L,1L,2L),2))
  expect_equal(prepared$offset,matrix(rep(log(data$exposure),2),36,2,
    dimnames=list(NULL,c("red","blue"))))
  expect_equal(prepared$nominal$rowMap$originalRow,rep(1:36,2))
  expect_equal(prepared$nominal$rowMap$category,rep(c("red","blue"),each=36))
  expect_identical(prepared$obsInfo$included,rep(TRUE,36))
  expect_equal(dim(prepared$randomDesign),c(72L,12L))
  expect_equal(length(prepared$covStruct),1L)
  expect_equal(length(prepared$initialParameters),9L)
  expect_false(any(grepl("units",names(prepared$initialParameters))))
})

test_that("nominal reference and response representations cannot be silently reordered", {
  data <- nominalPreparationData()
  data$y <- factor(data$y,levels=c("blue","red","green"))
  prepared <- prepareNominalFixture(data)
  expect_identical(prepared$nominal$reference,"green")
  expect_identical(prepared$nominal$nativeLevels,levels(data$y))
  data$y <- as.integer(data$y)
  numeric <- prepareNominalFixture(data,reference="1")
  expect_identical(numeric$nominal$nativeLevels,c("2","3","1"))
  expect_equal(numeric$nominal$response,match(as.character(data$y),c("2","3","1")))
  for(value in list(as.character(data$y),ordered(data$y),matrix(data$y),
    rep(NA_real_,36),rep(1,36),rep(c(1,3),18),rep(c(1,2.2,3),12),
    rep(c(0,1,2),12),rep(c(1,Inf,3),12))){
    data$y <- value
    expect_error(prepareNominalFixture(data),"Nominal")
  }
  data <- nominalPreparationData()
  for(reference in list("missing",1,c("red","blue"),NA_character_)){
    expect_error(prepareNominalFixture(data,reference=reference),"reference")
  }
  data$y <- factor(data$y,levels=c(levels(data$y),"absent"))
  expect_error(prepareNominalFixture(data),"every declared category")
  data <- nominalPreparationData()
  data$x[data$y=="red"] <- NA_real_
  expect_error(prepareNominalFixture(data),"every declared category")
})

test_that("nominal filtering keeps all category rows of an observation together", {
  data <- nominalPreparationData()
  data$y[1] <- NA
  data$x[2] <- NA_real_
  data$group[3] <- NA
  prepared <- prepareNominalFixture(data)
  included <- 4:36
  expect_equal(which(prepared$obsInfo$included),included)
  expect_equal(prepared$nominal$rowMap$originalRow,rep(included,2))
  expect_equal(prepared$nominal$response,as.numeric(data$y[included]))
  expect_equal(nrow(prepared$fixedDesign),66L)
  expect_equal(nrow(prepared$randomDesign),66L)
  expect_equal(nrow(prepared$offset),33L)
  expect_equal(prepared$nominal$counts,stats::setNames(tabulate(data$y[included],3),
    levels(data$y)))
  omitted <- prepareNominalFixture(data,naMethodY="omit",naMethodX="omit",
    naMethodRandom="omit")
  expect_equal(omitted$fixedDesign,prepared$fixedDesign)
  expect_equal(omitted$randomDesign,prepared$randomDesign)
  expect_error(prepareNominalFixture(data,naMethodY="fail"),"response")
  expect_error(prepareNominalFixture(data,naMethodX="fail"),"fixed-effect")
  expect_error(prepareNominalFixture(data,naMethodRandom="fail"),"random-effect")
  expect_error(prepareNominalFixture(data,naMethodY="include"),"missing-data methods")
  data <- nominalPreparationData()
  data$x[2] <- Inf
  expect_error(prepareNominalFixture(data),"finite")
  data <- nominalPreparationData()
  data$duplicate <- data$x
  expect_error(prepareNominalFixture(data,y~x+duplicate),"full-rank")
})

test_that("nominal fixed bases preserve intercept choices, transformations and dot expansion", {
  data <- nominalPreparationData()
  fixed <- y~0+poly(x,2)+treatment
  prepared <- prepareNominalFixture(data,fixed)
  base <- model.matrix(fixed,data)
  expect_equal(unname(as.matrix(prepared$fixedDesign)),
    unname(kronecker(diag(2),base)),tolerance=1e-12)
  expect_identical(prepared$nominal$fixed$columns,colnames(base))
  noFixed <- prepareNominalFixture(data,y~0)
  expect_equal(dim(noFixed$fixedDesign),c(72L,0L))
  dotData <- data[,c("y","x")]
  dotData$group <- rep(c(0,1),each=18)
  dot <- prepareNominalFixture(dotData,y~.)
  expect_identical(dot$nominal$fixed$columns,c("(Intercept)","x","group"))
  expect_false(any(grepl(".nominal",colnames(dot$fixedDesign),fixed=TRUE)))
})

test_that("nominal random designs explicitly distinguish shared and category-specific effects", {
  data <- nominalPreparationData()
  shared <- prepareNominalFixture(data,random=~group)
  independent <- prepareNominalFixture(data,
    random=~vsm(ism(.nominalCategory),ism(group)))
  correlated <- prepareNominalFixture(data,
    random=~vsm(usm(.nominalCategory),ism(group)))
  base <- model.matrix(~0+group,data)
  expect_equal(unname(as.matrix(shared$randomDesign)),unname(rbind(base,base)))
  expect_equal(unname(as.matrix(independent$randomDesign)),
    unname(kronecker(diag(2),base)))
  expect_equal(correlated$randomDesign,independent$randomDesign)
  expect_length(shared$initialParameters,5)
  expect_length(independent$initialParameters,5)
  expect_length(correlated$initialParameters,7)
  binary <- data
  binary$y <- factor(ifelse(data$y=="red","red","other"))
  binary <- prepareNominalFixture(binary,random=~group)
  expect_equal(dim(binary$fixedDesign),c(36L,2L))
  expect_equal(dim(binary$offset),c(36L,1L))
})

test_that("nominal preparation uses existing controls and protects unsupported public surfaces", {
  data <- nominalPreparationData()
  for(control in list(mmes_control(),
    mmes_control(formulation="direct"),
    mmes_control(task=list(objective="ml",glmm=list(method="laplace")),
      output=list(inverse=1L)))){
    expect_error(prepareNominalFixture(data,control=control),"Nominal|nominal|Laplace")
  }
  expect_error(prepareNominalFixture(data,
    control=mmes_control(task=list(objective="ml",glmm=list(method="laplace")),
      score=list(denseCrossproduct="force"))),"Gaussian engine controls")
  expect_error(prepareNominalFixture(data,
    control=mmes_control(task=list(objective="ml",glmm=list(method="laplace",
      familyParameters=list(estimate=FALSE))))),"familyParameters")
  expect_error(prepareNominalFixture(data,
    control=mmes_control(task=list(objective="ml",glmm=list(method="laplace")),
      execution=list(denseMemoryMB=1e-6))),"denseMemoryMB")
  data$.nominalCategory <- "reserved"
  expect_error(prepareNominalFixture(data),"reserves")
  expect_error(prepareNominalFixture(data=NULL),"data frame")
  expect_error(prepareNominalFixture(fixed=y~.nominalCategory),"reserves")
  expect_false("nominal" %in% getNamespaceExports("sommer"))
  expect_false(".laplace_prepare_nominal" %in% getNamespaceExports("sommer"))
  expect_error(sommer:::.laplace_nominal_objective(list()),"preparation")
})

test_that("nominal descriptor objective matches independent correlated precision and marginal scores", {
  data <- nominalPreparationData()
  K <- matrix(c(1,.3,.3,1.4),2)
  Gu <- diag(6)
  Gu[1,2] <- Gu[2,1] <- .15
  dimnames(Gu) <- list(levels(data$group),levels(data$group))
  attr(Gu,"inverse") <- TRUE
  prepared <- prepareNominalFixture(data,y~x+offset(log(exposure)),
    random=~vsm(usm(.nominalCategory,theta=K),ism(group),sigma2=.6,Gu=Gu))
  objective <- sommer:::.laplace_nominal_objective(prepared)
  parameters <- objective$initialParameters
  parameters[1:4] <- c(.1,.2,-.2,.3)
  evaluation <- objective$evaluate(parameters)
  expect_true(evaluation$valid,info=evaluation$reason)
  expect_identical(names(evaluation$gradient),names(parameters))
  precision <- Matrix::Matrix(kronecker(solve(.6*K),Gu),sparse=TRUE)
  context <- sommer:::sommer_multinomial_mode_create(prepared$nominal$response,
    prepared$randomDesign,precision,2L)
  eta <- prepared$offset+matrix(prepared$fixedDesign %*% parameters[1:4],36,2)
  reference <- sommer:::sommer_multinomial_laplace_eval(context,eta,precision,
    prepared$fixedDesign,list(),rep(0,12),100L,1e-10,TRUE)
  expect_equal(evaluation$logLik,reference$logLik,tolerance=1e-10)
  expect_equal(evaluation$mode,reference$mode,tolerance=1e-9)
  expect_equal(evaluation$modeHessian,reference$modeHessian,tolerance=1e-9)
  expect_equal(unname(evaluation$gradient[1:4]),as.numeric(reference$gradient),
    tolerance=1e-9)
  expect_equal(evaluation$work$relationshipFactorizations,1L)
  for(step in c(3e-4,1e-4)){
    numerical <- vapply(seq_along(parameters),function(k){
      plus <- minus <- parameters
      plus[k] <- plus[k]+step
      minus[k] <- minus[k]-step
      a <- objective$evaluate(plus,FALSE)
      b <- objective$evaluate(minus,FALSE)
      expect_true(a$valid,info=a$reason)
      expect_true(b$valid,info=b$reason)
      (a$logLik-b$logLik)/(2*step)
    },numeric(1))
    expect_equal(unname(evaluation$gradient),numerical,tolerance=2e-6)
  }
  repeatValue <- objective$evaluate(parameters)
  expect_equal(repeatValue$gradient,evaluation$gradient,tolerance=1e-10)
  expect_equal(repeatValue$work$relationshipFactorizations,1L)
  expect_gt(repeatValue$work$descriptorReuses,0L)
  expect_equal(repeatValue$work$symbolicAnalyses,1L)
})

test_that("nominal coordinate plumbing handles multiple terms, fixed flags and invalid proposals", {
  data <- nominalPreparationData()
  prepared <- prepareNominalFixture(data,
    random=~vsm(usm(.nominalCategory,fixed=c(TRUE,FALSE)),ism(group),
      sigma2=.4,fixedSigma2=TRUE)+
      vsm(ism(.nominalCategory),ism(block),sigma2=.2))
  original <- prepared
  objective <- sommer:::.laplace_nominal_objective(prepared)
  p <- objective$initialParameters
  expect_length(p,6)
  value <- objective$evaluate(p)
  expect_true(value$valid,info=value$reason)
  expect_equal(value$work$relationshipFactorizations,2L)
  for(k in seq_along(p)){
    h <- 1e-4
    plus <- minus <- p
    plus[k] <- plus[k]+h
    minus[k] <- minus[k]-h
    numeric <- (objective$evaluate(plus,FALSE)$logLik-
      objective$evaluate(minus,FALSE)$logLik)/(2*h)
    expect_equal(unname(value$gradient[k]),numeric,tolerance=2e-6)
  }
  bad <- p
  bad[5] <- 1000
  failed <- objective$evaluate(bad)
  expect_false(failed$valid)
  expect_match(failed$reason,"exponential")
  expect_false("gradient" %in% names(failed))
  failed <- objective$evaluate(rep(NA_real_,length(p)))
  expect_false(failed$valid)
  expect_match(failed$reason,"non-finite")
  recovered <- objective$evaluate(p)
  expect_equal(recovered$logLik,value$logLik,tolerance=1e-10)
  expect_equal(recovered$gradient,value$gradient,tolerance=1e-9)
  expect_identical(prepared,original)
  off <- objective$evaluate(p,FALSE)
  expect_equal(off$work$scoreEvaluations,recovered$work$scoreEvaluations)
  expect_equal(off$work$evaluations,recovered$work$evaluations)
  expect_false("gradient" %in% names(off))
  repeated <- objective$evaluate(p)
  expect_equal(repeated$work$evaluations,off$work$evaluations)
  expect_equal(repeated$work$derivativeSolves,off$work$derivativeSolves)
  repeated$probability[,] <- 0
  expect_equal(objective$evaluate(p)$probability,recovered$probability)
  changedTolerance <- sommer:::sommer_multinomial_objective_eval(objective$context,p,
    tolerance=1e-11,gradient=FALSE)
  expect_gt(changedTolerance$work$evaluations,off$work$evaluations)
  expect_error(objective$evaluate(p[-1]),"one vector")
  expect_error(objective$evaluate(matrix(p,ncol=1)),"one vector")
  expect_error(objective$evaluate(p,NA),"TRUE/FALSE")
  serialized <- unserialize(serialize(objective$context,NULL))
  expect_error(sommer:::sommer_multinomial_objective_eval(serialized,p),"fresh context")
  expect_error(sommer:::sommer_multinomial_objective_eval(NULL,p),"Invalid multinomial")
  expect_error(sommer:::sommer_multinomial_objective_eval(objective$context,p[-1]),
    "incompatible length")
  expect_error(sommer:::sommer_multinomial_objective_eval(objective$context,p,0L),
    "positive")
  set.seed(9172)
  before <- .Random.seed
  invisible(objective$evaluate(p))
  expect_identical(before,.Random.seed)
})

test_that("nominal mode limits and failed proposals do not cache partial results", {
  prepared <- prepareNominalFixture()
  prepared$control$task$glmm$modeMaxit <- 1L
  objective <- sommer:::.laplace_nominal_objective(prepared)
  p <- objective$initialParameters
  p[1:4] <- c(.7,.5,-.4,.3)
  failed <- objective$evaluate(p)
  expect_false(failed$valid)
  expect_match(failed$reason,"iteration limit")
  expect_false("gradient" %in% names(failed))
  retry <- objective$evaluate(p)
  expect_gt(retry$work$evaluations,failed$work$evaluations)
  recovered <- sommer:::sommer_multinomial_objective_eval(objective$context,p,
    maxModeIterations=100L,tolerance=1e-10)
  expect_true(recovered$valid,info=recovered$reason)
  expect_true(all(is.finite(recovered$gradient)))
})

test_that("prepared binary nominal objectives agree with the scalar descriptor objective", {
  data <- nominalPreparationData()
  data$y <- factor(ifelse(data$y=="red","red","other"),levels=c("red","other"))
  prepared <- prepareNominalFixture(data,y~x+offset(log(exposure)),random=~group)
  objective <- sommer:::.laplace_nominal_objective(prepared)
  parameters <- objective$initialParameters
  parameters[1:2] <- c(.3,-.4)
  binary <- sommer:::sommer_laplace_objective_create(prepared$fixedDesign,
    prepared$randomDesign,as.numeric(data$y=="red"),as.numeric(prepared$offset),
    "binomial",prepared$covStruct,prepared$relationships)
  expected <- sommer:::sommer_laplace_objective_eval(binary,parameters,tolerance=1e-10)
  value <- objective$evaluate(parameters)
  expect_true(value$valid)
  expect_equal(value$logLik,expected$logLik,tolerance=1e-10)
  expect_equal(value$mode,expected$mode,tolerance=1e-9)
  expect_equal(unname(value$gradient),as.numeric(expected$gradient),tolerance=1e-9)

  noFixed <- prepareNominalFixture(data,y~0,
    random=~vsm(ism(group),sigma2=.3,fixedSigma2=TRUE))
  fixedObjective <- sommer:::.laplace_nominal_objective(noFixed)
  expect_length(fixedObjective$initialParameters,0)
  empty <- fixedObjective$evaluate(numeric())
  expect_true(empty$valid)
  expect_length(empty$gradient,0)
  expect_equal(empty$work$derivativeSolves,0L)
  expect_equal(empty$work$relationshipFactorizations,1L)
})

test_that("the internal nominal objective can jointly optimize fixed and covariance coordinates", {
  set.seed(7329)
  groups <- 24L
  perGroup <- 12L
  data <- data.frame(group=factor(rep(seq_len(groups),each=perGroup)),
    x=rep(seq(-1,1,length.out=perGroup),groups))
  effects <- matrix(rnorm(groups*2,sd=.7),groups,2)
  eta <- cbind(.3+.5*data$x,-.4-.3*data$x)+effects[as.integer(data$group),]
  weights <- exp(cbind(eta,0))
  probabilities <- weights/rowSums(weights)
  outcomes <- vapply(seq_len(nrow(data)),function(i)
    sample.int(3L,1L,prob=probabilities[i,]),integer(1))
  data$y <- factor(outcomes,levels=1:3)
  prepared <- prepareNominalFixture(data,
    random=~vsm(ism(.nominalCategory),ism(group),sigma2=.5))
  objective <- sommer:::.laplace_nominal_objective(prepared)
  fit <- nlminb(objective$initialParameters,
    function(parameters) -objective$evaluate(parameters,FALSE)$logLik,
    function(parameters) -objective$evaluate(parameters,TRUE)$gradient,
    control=list(iter.max=100,eval.max=400,rel.tol=1e-12,x.tol=1e-10))
  value <- objective$evaluate(fit$par)
  expect_equal(fit$convergence,0L,info=fit$message)
  expect_true(value$valid,info=value$reason)
  expect_lt(max(abs(value$gradient)),1e-3)
  expect_true(is.finite(exp(tail(fit$par,1))))
  expect_gt(exp(tail(fit$par,1)),.01)
  expect_equal(rowSums(value$probability),rep(1,nrow(data)),tolerance=1e-12)
  expect_false(inherits(value,"mmes"))
})
