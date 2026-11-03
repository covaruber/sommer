categoricalMixedFixture <- function(composition=c("ord", "nom", "count"), slopes=TRUE){
  set.seed(940)
  n <- 120L
  data <- data.frame(trait=factor(rep(composition, each=n), levels=composition),
    x=rnorm(n*length(composition)),
    group=factor(rep(rep(1:20, each=6), length(composition))))
  effect <- rnorm(20, sd=sqrt(0.15))[data$group]
  eta <- 0.35*data$x + effect
  data$y <- numeric(nrow(data))
  for(trait in composition){
    rows <- which(data$trait == trait)
    data$y[rows] <- switch(trait,
      ord=as.integer(cut(eta[rows]+rlogis(n), c(-Inf,-0.8,0.9,Inf))),
      nom={
        logits <- cbind(0.3+eta[rows], -0.2-0.2*data$x[rows]+effect[rows], 0)
        p <- exp(logits)/rowSums(exp(logits))
        1L+rowSums(runif(n) > t(apply(p, 1, cumsum))[,1:2])
      },
      count=rpois(n, exp(0.5+eta[rows])))
  }
  families <- list(ord=ordinal(), nom=multinomial(reference="2"), count=poisson())[composition]
  family <- do.call(familym, c(list(by="trait"), families))
  control <- mmes_control(task=list(objective="ml",glmm=list(method="laplace")),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  fixed <- if(slopes) y~trait+trait:x-1 else y~trait-1
  fit <- mmes(fixed, random=~vsm(ism(group),sigma2=0.15,fixedSigma2=TRUE),
    data=data, family=family, control=control)
  list(fit=fit,data=data,family=family,control=control)
}

test_that("familym fits ordinal, nominal and both with scalar likelihoods", {
  for(composition in list(c("ord","count"),c("nom","count"),c("ord","nom","count"))){
    fixture <- categoricalMixedFixture(composition)
    fit <- fixture$fit
    expect_true(fit$convergence)
    expect_identical(nobs(fit), nrow(fixture$data))
    expect_null(fit$Ci)
    expect_null(fit$fixedCovariance)
    expect_null(fit$uPevList)
    expect_false("hessian" %in% names(fit$optimizer))
    expect_true(all(is.na(fit$betaSummary[,"Std. Error"])))
    expect_identical(names(fitted(fit)), composition)
    expect_identical(fit$executionPlan$glmm$likelihood,"mixed scalar/category blocks")
    expect_match(fit$executionPlan$glmm$randomExpansion,"explicit formula")
    path <- fit$executionPlan$glmm$path
    expectedPath <- if("nom" %in% composition)
      "native-laplace-category-block" else "native-laplace-scalar"
    expect_true(path$eligible)
    expect_identical(path$requested,expectedPath)
    expect_identical(path$selected,expectedPath)
    expect_identical(path$reason,"eligible")
    expect_identical(path$likelihood,if("nom" %in% composition) "category-block" else "scalar")
    expect_identical(fit$executionPlan$glmm$selectedPath,path$selected)
    expect_identical(path,fit$glmmPathPlan)
    expect_true(is.numeric(fitted(fit)$count))
    for(trait in intersect(composition,c("ord","nom"))){
      p <- fitted(fit)[[trait]]
      expect_identical(dim(p),c(120L,3L))
      expect_identical(colnames(p),as.character(1:3))
      expect_equal(rowSums(p),rep(1,120),ignore_attr=TRUE,tolerance=1e-12)
    }
    expect_equal(predict(fit,type="response"),fitted(fit))
    expect_equal(predict(fit,type="link"),fitted(fit,type="link"))
    expect_equal(AIC(fit),-2*as.numeric(logLik(fit))+2*fit$npar)
    expect_equal(BIC(fit),-2*as.numeric(logLik(fit))+log(nobs(fit))*fit$npar)
    # The original shared effect stays one effect per group, not per logit.
    expect_identical(ncol(fit$randomDesign),20L)
    expect_equal(fit$randomBlockWidths,20L)
    expect_true(all(is.finite(unlist(residuals(fit,type="pearson")))))
    expect_true(all(is.finite(unlist(residuals(fit,type="deviance")))))
    if("nom" %in% composition){
      expect_error(residuals(fit,type="quantile"),"unordered categories")
      expect_error(residuals(fit,type="working"),"unordered categories")
    }else expect_true(all(is.finite(unlist(residuals(fit,type="quantile")))))
  }
})

test_that("nominal familym uses one block likelihood with off-diagonal curvature", {
  fixture <- categoricalMixedFixture(c("nom","count"))
  fit <- fixture$fit
  specification <- fit$familySpecifications$nom
  groups <- specification$nominalGroups
  evaluation <- sommer:::sommer_glmm_family_eval(fit$laplaceNativeResponse,
    fit$linear.predictors,fit$familySpecifications,fit$laplaceFamilyIndex,fit$priorWeights)
  logits <- matrix(fit$linear.predictors[as.vector(t(groups))], ncol=2,byrow=TRUE)
  p <- sommer:::.laplace_nominal_probabilities(logits,specification)
  response <- fit$y[groups[,1]]
  manual <- sum(log(p[cbind(seq_len(nrow(p)),
    match(specification$nativeLevels[response], specification$levels))])) +
      sum(dpois(fit$y[fit$laplaceFamilyIndex==2],
        fitted(fit)$count,log=TRUE))
  expect_equal(evaluation$logDensity,manual,tolerance=1e-10)
  native <- match(specification$predictorLevels,specification$levels)
  expect_equal(as.numeric(evaluation$blockCurvature[cbind(groups[,1],groups[,2])]),
    -p[,native[1]]*p[,native[2]],ignore_attr=TRUE,tolerance=1e-12)
  expect_equal(fit$laplaceObservationCurvature,evaluation$blockCurvature)
  expect_equal(fit$deviance,-2*manual)
})

test_that("mixed prediction preserves traits, row order, shared modes and offsets", {
  fixture <- categoricalMixedFixture()
  fit <- fixture$fit
  rows <- c(244,122,3,241,120,121)
  newdata <- fixture$data[rows,]
  predicted <- predict(fit,newdata=newdata,type="response")
  expect_equal(unname(predicted$ord),unname(fitted(fit)$ord[c(3,120),]))
  expect_equal(unname(predicted$nom),unname(fitted(fit)$nom[c(2,1),]))
  expect_equal(predicted$count,fitted(fit)$count[c(4,1)])
  expect_equal(predict(fit,newdata=newdata,re.form=NA,type="response"),
    sommer:::.laplace_mixed_output(fit,
      sommer:::.laplace_mixed_prediction_design(fit,newdata,FALSE)$offset +
      as.numeric(sommer:::.laplace_mixed_prediction_design(fit,newdata,FALSE)$fixed %*% fit$b),
      index=match(newdata$trait,names(fit$familySpecifications))[
        sommer:::.laplace_mixed_expand(newdata,fit$family,fit$familySpecifications)$rowMap],
      rowMap=sommer:::.laplace_mixed_expand(newdata,fit$family,fit$familySpecifications)$rowMap,
      coordinate=sommer:::.laplace_mixed_expand(newdata,fit$family,fit$familySpecifications)$coordinate))
  expect_identical(dim(predict(fit,newdata=fixture$data[1:3,])$nom),c(0L,3L))
  unseen <- newdata
  unseen$group <- "new"
  expect_error(predict(fit,newdata=unseen),"unseen")
  expect_true(all(is.finite(unlist(predict(fit,newdata=unseen,re.form=NA)))))
  population <- predict(fit,newdata=unseen,re.form="population")
  expect_equal(rowSums(population$nom),rep(1,2),ignore_attr=TRUE,tolerance=1e-10)
  expect_true(all(is.finite(unlist(population))))
})

test_that("mixed inference, simulation and bootstrap stay post-fit", {
  fixture <- categoricalMixedFixture()
  fit <- fixture$fit
  expect_true(all(is.finite(vcov(fit))))
  summary <- summary(fit)
  expect_true(all(is.finite(summary$betas[, "Std.Error"])))
  uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
  expect_true(uncertainty$valid)
  expect_true(all(is.finite(uncertainty$thresholds$StdError)))
  pev <- postPEV(fit,mode=2L)
  expect_true(all(is.finite(pev$Ci)))
  expect_true(all(unlist(pev$uPevList)>=0))
  newdata <- fixture$data[c(2,121,241),]
  for(mode in c("conditional","full")){
    predicted <- predict(fit,newdata=newdata,se.fit=TRUE,cov.fit=TRUE,uncertainty=mode)
    expect_equal(predicted$fit,predict(fit,newdata=newdata))
    expect_true(all(is.finite(unlist(predicted$se.fit))))
    expect_true(all(diag(predicted$cov.fit)>=0))
    expect_equal(predicted$cov.fit,t(predicted$cov.fit))
  }
  for(re.form in c("conditional","population")){
    simulations <- simulate(fit,nsim=2L,seed=3,re.form=re.form)
    expect_identical(dim(simulations),c(360L,2L))
    expect_true(all(as.matrix(simulations[fixture$data$trait!="count",]) %in% 1:3))
    expect_equal(simulations,simulate(fit,nsim=2L,seed=3,re.form=re.form))
  }
  boot <- boot_mmes(fit,nsim=2L,seed=41,maxit=300L)
  expect_true(all(boot$diagnostics$valid))
})

test_that("mixed comparisons and explicit random coordinates are supported", {
  fixture <- categoricalMixedFixture()
  small <- mmes(y~trait-1,random=~vsm(ism(group),sigma2=.15,fixedSigma2=TRUE),
    data=fixture$data,family=fixture$family,control=fixture$control)
  comparison <- anova(small,fixture$fit)
  expect_true(all(is.finite(comparison[,"logLik"])))
  bootstrap <- boot_mmes(small,fixture$fit,nsim=2L,seed=12,maxit=300L)
  expect_true(all(bootstrap$diagnostics$valid))
  explicit <- mmes(y~trait+trait:x-1,
    random=~vsm(ism(.familyCoordinate),ism(group),sigma2=.15,fixedSigma2=TRUE),
    data=fixture$data,family=fixture$family,control=fixture$control)
  expect_true(explicit$convergence)
  expect_identical(ncol(explicit$randomDesign),80L)
  predicted <- predict(explicit,newdata=fixture$data[c(2,121,241),])
  expect_equal(unname(predicted$ord),unname(fitted(explicit)$ord[2,,drop=FALSE]))
  expect_equal(unname(predicted$nom),unname(fitted(explicit)$nom[1,,drop=FALSE]))
  expect_equal(predicted$count,fitted(explicit)$count[1])
  population <- predict(explicit,newdata=fixture$data[c(2,121,241),],
    re.form="population")
  expect_true(all(is.finite(unlist(population))))
  expect_equal(rowSums(population$nom),1,ignore_attr=TRUE,tolerance=1e-10)
  expect_error(mmes(y~trait,random=~vsm(ism(.nominalCategory),ism(group)),
    data=fixture$data,family=fixture$family,control=fixture$control),"ambiguous")
})

test_that("categorical familym enforces ML, response coding and memory policy", {
  fixture <- categoricalMixedFixture(c("ord","count"))
  data <- fixture$data
  family <- fixture$family
  pql <- mmes_control(task=list(objective="ml",glmm=list(method="pql")),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  expect_error(mmes(y~trait,data=data,family=family,control=pql),"not supported by PQL")
  data$y[1] <- 5
  expect_error(mmes(y~trait,random=~vsm(ism(group)),data=data,
    family=family,control=fixture$control),"contiguous")
  limited <- fixture$fit
  limited$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(fitted(limited),"denseMemoryMB")
  expect_error(predict(limited,newdata=fixture$data[1:3,]),"denseMemoryMB")
  expect_error(vcov(limited),"denseMemoryMB")
  expect_error(postPEV(limited,1L),"denseMemoryMB")
})

test_that("categorical composition retains scalar shapes, trials and missingness", {
  fixture <- categoricalMixedFixture(c("ord","count"))
  scalarFamilies <- list(binomial=binomial("probit"), gaussian=gaussian(),
    Gamma=Gamma("log"), nbinom2=nbinom2(2),
    ziPoisson=zipoisson(zero=~x))
  for(code in names(scalarFamilies)){
    data <- fixture$data
    rows <- which(data$trait=="count")
    data$trials <- 1
    settings <- list()
    if(code == "binomial"){
      data$trials[rows] <- 5
      data$y[rows] <- rep(0:5,length.out=length(rows))
    }else if(code == "gaussian"){
      data$y[rows] <- log1p(data$y[rows])-1
      settings <- list(count=list(dispersion=0.7,estimate=FALSE))
    }else if(code == "Gamma"){
      data$y[rows] <- data$y[rows]+0.5
      settings <- list(count=list(dispersion=0.7,estimate=FALSE))
    }else if(code == "nbinom2"){
      settings <- list(count=list(theta=2,estimate=FALSE))
    }else if(code == "ziPoisson"){
      data$y[rows[seq.int(1L,length(rows),by=3L)]] <- 0
    }
    family <- familym("trait",ord=ordinal(),count=scalarFamilies[[code]],
      trials=if(code=="binomial") "trials" else NULL)
    control <- mmes_control(task=list(objective="ml",glmm=list(method="laplace",
      familyParameters=settings)),diagnostics=list(verbose=FALSE,dateWarning=FALSE))
    fit <- mmes(y~trait+trait:x-1,random=~vsm(ism(group),sigma2=.15,fixedSigma2=TRUE),
      data=data,family=family,control=control)
    expect_true(fit$convergence)
    expect_true(all(is.finite(unlist(fitted(fit)))))
    expect_true(all(is.finite(unlist(predict(fit,newdata=data[c(1,121),])))))
    expect_true(all(is.finite(unlist(residuals(fit,type="pearson")))))
    expect_identical(dim(simulate(fit,2,seed=2)),c(240L,2L))
  }
  fixture <- categoricalMixedFixture()
  data <- fixture$data
  data$x[c(1,121,241)] <- NA
  fit <- mmes(y~trait+trait:x-1,random=~vsm(ism(group),sigma2=.15,fixedSigma2=TRUE),
    data=data,family=fixture$family,control=fixture$control)
  expect_identical(nobs(fit),357L)
  expect_equal(vapply(fitted(fit),NROW,integer(1)),c(ord=119L,nom=119L,count=119L))
  expect_identical(dim(simulate(fit,2,seed=2)),c(357L,2L))
})

test_that("binary nominal traits and fitted contrasts retain their coordinate maps", {
  fixture <- categoricalMixedFixture()
  data <- fixture$data
  data$y[data$trait=="nom"] <- rep(1:2,60)
  data$treatment <- factor(rep(c("a","b","c"),120))
  fit <- mmes(y~trait+trait:x+trait:treatment-1,
    random=~vsm(ism(group),sigma2=.15,fixedSigma2=TRUE),data=data,
    family=fixture$family,control=fixture$control,contrasts=list(treatment="contr.sum"))
  expect_true(fit$convergence)
  expect_identical(dim(fitted(fit)$nom),c(120L,2L))
  expect_identical(dim(fit$familySpecifications$nom$nominalGroups),c(120L,1L))
  predicted <- predict(fit,newdata=data[c(1,122,241),])
  expect_equal(unname(predicted$nom),unname(fitted(fit)$nom[2,,drop=FALSE]))
  expect_equal(unname(predicted$ord),unname(fitted(fit)$ord[1,,drop=FALSE]))
})

test_that("mixed fitting never computes post-fit inference and filters ZIP-only covariates", {
  local_mocked_bindings(
    .mmes_laplace_uncertainty=function(...) stop("Fit-time uncertainty"),
    .mmes_laplace_postPEV=function(...) stop("Fit-time PEV"),
    .package="sommer")
  local_mocked_bindings(optimHess=function(...) stop("Fit-time Hessian"),.package="stats")
  fixture <- categoricalMixedFixture()
  expect_true(fixture$fit$convergence)
  data <- fixture$data
  rows <- which(data$trait=="count")
  data$z <- NA_real_
  data$z[rows] <- data$x[rows]
  data$z[rows[1]] <- NA
  data$y[rows[seq.int(1L,length(rows),by=3L)]] <- 0
  family <- familym("trait",ord=ordinal(),nom=multinomial(reference="2"),
    count=zipoisson(zero=~z))
  fit <- mmes(y~trait+trait:x-1,random=~vsm(ism(group),sigma2=.15,fixedSigma2=TRUE),
    data=data,family=family,control=fixture$control)
  expect_true(fit$convergence)
  expect_identical(nobs(fit),359L)
  expect_equal(vapply(fitted(fit),NROW,integer(1)),c(ord=120L,nom=120L,count=119L))
  expect_true(all(is.finite(unlist(predict(fit,newdata=data[c(1,121,242),])))))
})
