nbinom2 <- function(theta=1, link="log"){
  theta <- .mme_positive(theta, "theta")
  if(!identical(link, "log")){
    stop("nbinom2() currently supports link='log' only.", call.=FALSE)
  }
  linkFunctions <- stats::make.link(link)
  shapeValue <- theta
  result <- list(family="Negative Binomial NB2", link=link,
    theta=shapeValue, linkfun=linkFunctions$linkfun, linkinv=linkFunctions$linkinv,
    mu.eta=linkFunctions$mu.eta,
    variance=function(meanValue) meanValue + meanValue^2 / shapeValue,
    dev.resids=function(responseValue, meanValue, weightValue){
      responseTerm <- numeric(length(responseValue))
      positiveResponse <- responseValue > 0
      responseTerm[positiveResponse] <- responseValue[positiveResponse] *
        log(responseValue[positiveResponse] / meanValue[positiveResponse])
      2 * weightValue * (responseTerm - (responseValue + shapeValue) *
        log((responseValue + shapeValue) / (meanValue + shapeValue)))
    },
    aic=function(responseValue, trialsValue, meanValue, weightValue, devianceValue)
      -2 * sum(weightValue * stats::dnbinom(responseValue, size=shapeValue,
        mu=meanValue, log=TRUE)),
    initialize=expression({
      if(any(y < 0)) stop("Negative-binomial responses must be non-negative.")
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    }),
    validmu=function(meanValue) all(is.finite(meanValue)) && all(meanValue > 0),
    valideta=linkFunctions$valideta)
  class(result) <- "family"
  result
}

zipoisson <- function(pi=0.1, link="log", zero=~1){
  if(!is.numeric(pi) || length(pi) != 1L || is.na(pi) ||
     !is.finite(pi) || pi <= 0 || pi >= 1){
    stop("zipoisson() requires one finite starting probability strictly between 0 and 1.",
         call.=FALSE)
  }
  if(!inherits(zero, "formula") || length(zero) != 2L ||
     any(c("|", "||", "vsm", "ism") %in% all.names(zero))){
    stop("zipoisson(zero=) requires a one-sided fixed-effects formula without random effects.",
      call.=FALSE)
  }
  result <- .zipoisson_family(pi, link)
  result$zero <- zero
  result
}

.laplace_zip_zero_log <- function(probability, countMean){
  inflation <- log(probability)
  poisson <- log1p(-probability) - countMean
  larger <- pmax(inflation, poisson)
  larger + log1p(exp(pmin(inflation, poisson) - larger))
}

.zipoisson_family <- function(pi, link){
  if(!is.numeric(pi) || !length(pi) || anyNA(pi) ||
     any(!is.finite(pi)) || any(pi < 0 | pi >= 1)){
    stop("zipoisson() requires finite zero-inflation probabilities in [0, 1).",
         call.=FALSE)
  }
  if(!identical(link, "log")){
    stop("zipoisson() currently supports link='log' only.", call.=FALSE)
  }
  linkFunctions <- stats::make.link(link)
  countProbability <- 1 - pi
  result <- list(family="ziPoisson", link=link, pi=pi,
    linkfun=function(meanValue) log(meanValue / countProbability),
    linkinv=function(predictor) countProbability * linkFunctions$linkinv(predictor),
    mu.eta=function(predictor) countProbability * linkFunctions$mu.eta(predictor),
    variance=function(meanValue)
      meanValue + (pi / countProbability) * meanValue^2,
    dev.resids=function(responseValue, meanValue, weightValue){
      countMean <- meanValue / (1 - pi)
      logLikelihood <- numeric(length(responseValue))
      zeroRows <- responseValue == 0
      probability <- rep_len(pi, length(responseValue))
      logLikelihood[zeroRows] <- .laplace_zip_zero_log(probability[zeroRows],
        countMean[zeroRows])
      positiveRows <- !zeroRows
      logLikelihood[positiveRows] <- log1p(-probability[positiveRows]) +
        stats::dpois(responseValue[positiveRows],
          lambda=countMean[positiveRows], log=TRUE)
      saturated <- numeric(length(responseValue))
      saturated[positiveRows] <- stats::dpois(responseValue[positiveRows],
        lambda=responseValue[positiveRows], log=TRUE)
      2 * weightValue * (saturated - logLikelihood)
    },
    aic=function(responseValue, trialsValue, meanValue, weightValue, devianceValue){
      countMean <- meanValue / (1 - pi)
      zeroRows <- responseValue == 0
      logLikelihood <- numeric(length(responseValue))
      probability <- rep_len(pi, length(responseValue))
      logLikelihood[zeroRows] <- .laplace_zip_zero_log(probability[zeroRows],
        countMean[zeroRows])
      positiveRows <- !zeroRows
      logLikelihood[positiveRows] <- log1p(-probability[positiveRows]) +
        stats::dpois(responseValue[positiveRows],
          lambda=countMean[positiveRows], log=TRUE)
      -2 * sum(weightValue * logLikelihood)
    },
    initialize=expression({
      if(any(y < 0)) stop("Poisson responses must be non-negative.")
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    }),
    validmu=function(meanValue) all(is.finite(meanValue)) && all(meanValue > 0),
    valideta=linkFunctions$valideta)
  class(result) <- "family"
  result
}

.laplace_family_specs <- function(family, settings){
  mixed <- inherits(family, "sommer_familym")
  families <- if(mixed) family$families else list(response=family)
  if(!is.list(settings) || (length(settings) &&
     (is.null(names(settings)) || anyNA(names(settings)) ||
      any(!nzchar(names(settings))) || anyDuplicated(names(settings))))){
    stop("task$glmm$familyParameters must be a uniquely named list.", call.=FALSE)
  }
  if(mixed){
    unknown <- setdiff(names(settings), names(families))
    if(length(unknown)){
      stop("Unknown familyParameters trait(s): ", paste(unknown, collapse=", "), ".", call.=FALSE)
    }
    familySettings <- lapply(names(families), function(trait) settings[[trait]])
  }else{
    familySettings <- list(settings)
  }
  specifications <- lapply(seq_along(families), function(familyPosition){
    familyObject <- families[[familyPosition]]
    code <- if(grepl("^Negative Binomial", familyObject$family)) "nbinom2" else familyObject$family
    allowed <- switch(code, binomial=c("logit", "probit", "cloglog"),
      poisson="log", ziPoisson="log", nbinom2="log", gaussian="identity",
      Gamma="log", character())
    if(!familyObject$link %in% allowed){
      stop("Unsupported native Laplace family/link: ", familyObject$family,
        "(", familyObject$link, ").", call.=FALSE)
    }
    supplied <- familySettings[[familyPosition]]
    if(is.null(supplied)) supplied <- list()
    parameterName <- switch(code, nbinom2="theta", ziPoisson="zeroInflation",
      "dispersion")
    hasParameter <- code %in% c("nbinom2", "ziPoisson", "gaussian", "Gamma")
    defaultValue <- if(code == "nbinom2"){
      if(!is.null(familyObject$theta)) familyObject$theta else
        1 / (familyObject$variance(1) - 1)
    }else if(code == "ziPoisson"){
      familyObject$pi
    }else 1
    configuration <- .mme_module(supplied,
      if(hasParameter) stats::setNames(list(defaultValue, TRUE), c(parameterName, "estimate"))
      else list(), paste0("familyParameters$", names(families)[familyPosition]))
    value <- if(hasParameter) .mme_positive(configuration[[parameterName]], parameterName) else 1
    if(code == "ziPoisson" && (length(value) != 1L || value >= 1)){
      stop("familyParameters$zeroInflation must be strictly between 0 and 1.",
           call.=FALSE)
    }
    free <- hasParameter && .mme_flag(configuration$estimate, "familyParameters$estimate")
    coordinate <- if(code == "ziPoisson") stats::qlogis(value) else log(value)
    list(family=code, link=familyObject$link, logParameter=coordinate,
      free=free, parameter=if(hasParameter) parameterName else NULL,
      supplied=parameterName %in% names(supplied))
  })
  names(specifications) <- names(families)
  list(families=families, specifications=specifications, mixed=mixed)
}

.laplace_family_index <- function(family, data, observations){
  if(!inherits(family, "sommer_familym")) return(rep(1L, observations))
  if(is.null(data) || !family$by %in% names(data) || nrow(data) != observations){
    stop("familym() requires its trait column in data with one value per observation.", call.=FALSE)
  }
  traitValues <- as.character(data[[family$by]])
  familyIndex <- match(traitValues, names(family$families))
  unknown <- unique(traitValues[is.na(familyIndex) & !is.na(traitValues)])
  if(length(unknown)){
    stop("Traits without a family in familym(): ", paste(unknown, collapse=", "), ".", call.=FALSE)
  }
  familyIndex
}

.laplace_response <- function(response, included, specifications, familyIndex, declaredTrials=NULL){
  trials <- attr(response, "trials")
  if(!is.null(declaredTrials) && (!is.null(trials) || is.matrix(response))){
    stop("Declare grouped trials only once: familym(trials=...), binm(), or a two-column response.", call.=FALSE)
  }
  if(is.matrix(response)){
    if(ncol(response) != 2L ||
       any(vapply(specifications, function(specification)
         specification$family != "binomial", logical(1)))){
      stop("Matrix responses must be two-column binomial successes/failures.", call.=FALSE)
    }
    retainedCounts <- response[included, , drop=FALSE]
    if(any(!is.finite(retainedCounts)) || any(retainedCounts < 0) ||
       any(retainedCounts != floor(retainedCounts))){
      stop("Grouped-binomial successes and failures must be finite non-negative integers.", call.=FALSE)
    }
    trials <- rowSums(retainedCounts)
    response <- retainedCounts[, 1L] / trials
  }else{
    if(!is.numeric(response) || !is.null(dim(response))){
      stop("Native Laplace fitting requires one numeric response or two-column binomial counts.", call.=FALSE)
    }
    response <- as.numeric(response[included])
    trials <- if(!is.null(declaredTrials)) declaredTrials else
      if(is.null(trials)) rep(1, length(response)) else trials[included]
  }
  if(!is.numeric(trials) || any(!is.finite(response)) || any(!is.finite(trials)) || any(trials <= 0) ||
     any(trials != floor(trials))){
    stop("Retained responses must be finite; binomial trials must be positive integers.", call.=FALSE)
  }
  nativeResponse <- response
  for(familyPosition in seq_along(specifications)){
    rows <- familyIndex == familyPosition
    code <- specifications[[familyPosition]]$family
    if(code == "binomial"){
      if(!is.null(declaredTrials)){
        if(any(response[rows] != floor(response[rows]))){
          stop("Grouped-binomial successes must be integers.", call.=FALSE)
        }
        response[rows] <- response[rows] / trials[rows]
      }
      successes <- response[rows] * trials[rows]
      if(any(response[rows] < 0 | response[rows] > 1) ||
         any(trials[rows] == 1 & !response[rows] %in% c(0, 1)) ||
         any(abs(successes - round(successes)) >
           8 * .Machine$double.eps * pmax(1, trials[rows]))){
        stop("Binomial responses must be 0/1 or integer successes with declared trials.", call.=FALSE)
      }
      nativeResponse[rows] <- round(successes)
    }else{
      if(any(trials[rows] != 1)){
        stop("Trials may only be supplied for binomial observations.", call.=FALSE)
      }
      if(code %in% c("poisson", "ziPoisson", "nbinom2") &&
         any(response[rows] < 0 | response[rows] != floor(response[rows]))){
        stop("Native Laplace count responses must be non-negative integer counts.", call.=FALSE)
      }
      if(code == "Gamma" && any(response[rows] <= 0)){
        stop("Native Laplace Gamma responses must be strictly positive.", call.=FALSE)
      }
    }
  }
  list(response=response, nativeResponse=nativeResponse, trials=trials)
}

.laplace_zero_formula <- function(family){
  if(!identical(family$family, "ziPoisson") || is.null(family$zero)) return(NULL)
  terms <- stats::terms(family$zero)
  if(!length(attr(terms, "term.labels")) && attr(terms, "intercept") == 1L &&
     !length(attr(terms, "offset"))) return(NULL)
  family$zero
}

.laplace_zero_complete <- function(family, data, observations){
  families <- if(inherits(family, "sommer_familym")) family$families else list(family)
  index <- .laplace_family_index(family, data, observations)
  complete <- rep(TRUE, observations)
  for(position in seq_along(families)){
    formula <- .laplace_zero_formula(families[[position]])
    if(is.null(formula)) next
    rows <- which(index == position)
    frame <- stats::model.frame(formula, data=data[rows, , drop=FALSE],
      na.action=stats::na.pass)
    if(nrow(frame) != length(rows)){
      stop("Zero-component formula must describe one row per observation.", call.=FALSE)
    }
    complete[rows] <- stats::complete.cases(frame)
  }
  complete
}

.laplace_zero_designs <- function(families, specifications, data, included, index){
  for(position in seq_along(families)){
    formula <- .laplace_zero_formula(families[[position]])
    if(is.null(formula)) next
    if(!specifications[[position]]$free){
      stop("A zero-component formula requires estimate=TRUE; a fixed scalar probability uses zero=~1.",
        call.=FALSE)
    }
    if(is.null(data)) stop("A zero-component formula requires data.", call.=FALSE)
    rows <- which(index == position)
    frame <- stats::model.frame(formula,
      data=data[included, , drop=FALSE][rows, , drop=FALSE],
      na.action=stats::na.fail, drop.unused.levels=TRUE)
    terms <- attr(frame, "terms")
    design <- stats::model.matrix(terms, frame)
    if(!ncol(design) || any(!is.finite(design)) ||
       qr(design)$rank != ncol(design)){
      stop("The zero-component fixed design must be finite, nonempty and full rank.",
        call.=FALSE)
    }
    offset <- stats::model.offset(frame)
    if(is.null(offset)) offset <- numeric(length(rows))
    if(any(!is.finite(offset))) stop("Zero-component offsets must be finite.", call.=FALSE)
    fullDesign <- matrix(0, length(index), ncol(design),
      dimnames=list(NULL, colnames(design)))
    fullDesign[rows, ] <- design
    fullOffset <- numeric(length(index))
    fullOffset[rows] <- offset
    coefficients <- stats::setNames(numeric(ncol(design)), colnames(design))
    intercept <- match("(Intercept)", colnames(design))
    if(!is.na(intercept)) coefficients[intercept] <- specifications[[position]]$logParameter
    specifications[[position]]$zeroDesign <- fullDesign
    specifications[[position]]$zeroOffset <- fullOffset
    specifications[[position]]$zeroCoefficients <- coefficients
    specifications[[position]]$zeroTerms <- terms
    specifications[[position]]$zeroXlevels <- stats::.getXlevels(terms, frame)
    specifications[[position]]$zeroContrasts <- attr(design, "contrasts")
  }
  specifications
}

.laplace_family_count <- function(specification){
  if(!specification$free) return(0L)
  if(is.null(specification$zeroDesign)) 1L else ncol(specification$zeroDesign)
}

.laplace_zero_probability <- function(specification, rows=NULL){
  if(is.null(specification$zeroDesign)) return(stats::plogis(specification$logParameter))
  if(is.null(rows)) rows <- seq_len(nrow(specification$zeroDesign))
  stats::plogis(as.numeric(specification$zeroDesign[rows, , drop=FALSE] %*%
    specification$zeroCoefficients) + specification$zeroOffset[rows])
}

.laplace_prediction_specs <- function(object, newdata, index){
  specifications <- object$familySpecifications
  if(is.null(newdata)) return(specifications)
  for(position in seq_along(specifications)){
    specification <- specifications[[position]]
    if(is.null(specification$zeroDesign)) next
    rows <- which(index == position)
    fullDesign <- matrix(0, length(index), ncol(specification$zeroDesign),
      dimnames=list(NULL, colnames(specification$zeroDesign)))
    offset <- numeric(length(index))
    if(length(rows)){
      frame <- stats::model.frame(specification$zeroTerms,
        data=newdata[rows, , drop=FALSE], na.action=stats::na.fail,
        xlev=specification$zeroXlevels)
      design <- stats::model.matrix(specification$zeroTerms, frame,
        contrasts.arg=specification$zeroContrasts)
      if(!identical(colnames(design), colnames(fullDesign)) ||
         any(!is.finite(design))){
        stop("Zero-component prediction design differs from the fitted basis or is nonfinite.",
          call.=FALSE)
      }
      fullDesign[rows, ] <- design
      newOffset <- stats::model.offset(frame)
      if(!is.null(newOffset)) offset[rows] <- newOffset
    }
    if(any(!is.finite(offset))) stop("Zero-component prediction offsets must be finite.", call.=FALSE)
    specifications[[position]]$zeroDesign <- fullDesign
    specifications[[position]]$zeroOffset <- offset
  }
  specifications
}

.laplace_proposal_specs <- function(specifications, parameters, layout){
  counts <- vapply(specifications, .laplace_family_count, integer(1))
  starts <- cumsum(c(0L, head(counts, -1L)))
  for(position in which(counts > 0L)){
    values <- parameters[layout$fixedCount + layout$covarianceCount +
      starts[position] + seq_len(counts[position])]
    if(is.null(specifications[[position]]$zeroDesign)){
      specifications[[position]]$logParameter <- values[[1L]]
    }else specifications[[position]]$zeroCoefficients[] <- values
  }
  specifications
}

.laplace_zero_table <- function(specifications){
  do.call(rbind, lapply(seq_along(specifications), function(position){
    specification <- specifications[[position]]
    if(is.null(specification$zeroDesign)) return(NULL)
    data.frame(trait=names(specifications)[position],
      term=names(specification$zeroCoefficients),
      estimate=unname(specification$zeroCoefficients), StdError=NA_real_)
  }))
}

.laplace_update_families <- function(families, specifications, familyIndex=NULL){
  lapply(seq_along(families), function(familyPosition){
    if(specifications[[familyPosition]]$family == "nbinom2")
      nbinom2(exp(specifications[[familyPosition]]$logParameter))
    else if(specifications[[familyPosition]]$family == "ziPoisson"){
      rows <- if(is.null(familyIndex)) NULL else which(familyIndex == familyPosition)
      probability <- .laplace_zero_probability(specifications[[familyPosition]], rows)
      if(!length(probability)) probability <- stats::plogis(specifications[[familyPosition]]$logParameter)
      if(any(probability >= 1)){
        stop("The fitted zero-inflation probability reached the numerical boundary at one; a finite Poisson count mean is unavailable.",
             call.=FALSE)
      }
      result <- .zipoisson_family(probability, "log")
      result$zero <- families[[familyPosition]]$zero
      result
    }
    else families[[familyPosition]]
  })
}

.laplace_family_parameter_value <- function(specification){
  if(specification$family == "ziPoisson")
    stats::plogis(specification$logParameter) else
    exp(specification$logParameter)
}

.laplace_family_parameter_derivative <- function(specification){
  value <- .laplace_family_parameter_value(specification)
  if(specification$family == "ziPoisson") value * (1 - value) else value
}

.laplace_deviance <- function(response, mean, families, familyIndex, trials, dispersion){
  .pql_apply(families, familyIndex, "dev.resids", response, mean, trials) /
    dispersion[familyIndex]
}

.laplace_random_effect_sampling_factors <- function(object, maxMemoryMB){
  covarianceFactors <- vector("list", length(object$covStruct))
  for(termIndex in seq_along(object$covStruct)){
    blockWidth <- object$randomBlockWidths[termIndex]
    descriptorCovariance <- evaluate_covstruct_cpp(
      object$covStruct[[termIndex]], object$covStruct[[termIndex]]$par)
    relationshipPrecision <- Matrix::Matrix(
      object$randomRelationships[[termIndex]], sparse=TRUE)
    relationshipCount <- nrow(relationshipPrecision)
    descriptorCount <- nrow(descriptorCovariance)
    expectedWidth <- descriptorCount * relationshipCount
    if(expectedWidth != blockWidth){
      stop("Population random-effect simulation expected ", expectedWidth,
        " coefficients from descriptor/relationship dimensions, but term ",
        termIndex, " has ", blockWidth, " design columns.", call.=FALSE)
    }
    if(!all(dim(descriptorCovariance) == c(descriptorCount, descriptorCount)) ||
       !all(dim(relationshipPrecision) == c(relationshipCount,
         relationshipCount)) || any(!is.finite(descriptorCovariance))){
      stop("Population random-effect covariance metadata are invalid for term ",
        termIndex, ".", call.=FALSE)
    }
    workspaceBytes <- 24 * (as.double(relationshipCount) * descriptorCount +
      as.double(descriptorCount)^2) +
      24 * Matrix::nnzero(relationshipPrecision)
    workspaceMB <- workspaceBytes / 1048576
    if(!is.finite(workspaceMB) || workspaceMB > maxMemoryMB){
      stop("Population random-effect simulation needs approximately ",
        format(round(workspaceMB, 1), scientific=FALSE),
        " MB for sparse precision factors and draw workspaces, exceeding maxMemoryMB=",
        maxMemoryMB, ".", call.=FALSE)
    }
    relationshipFactor <- Matrix::Cholesky(relationshipPrecision,
      LDL=FALSE, perm=TRUE)
    covarianceFactors[[termIndex]] <- list(
      relationship=relationshipFactor,
      descriptor=t(chol(descriptorCovariance)),
      relationshipCount=relationshipCount,
      descriptorCount=descriptorCount)
  }
  covarianceFactors
}

.laplace_draw_random_effects <- function(covarianceFactors){
  unlist(lapply(covarianceFactors, function(factor){
    independent <- matrix(stats::rnorm(
      factor$relationshipCount * factor$descriptorCount),
      nrow=factor$relationshipCount,
      ncol=factor$descriptorCount)
    relationshipDraw <- solve(factor$relationship, independent,
      system="Lt")
    relationshipDraw <- solve(factor$relationship, relationshipDraw,
      system="Pt")
    as.numeric(relationshipDraw %*% t(factor$descriptor))
  }), use.names=FALSE)
}

.laplace_simulate_response <- function(object, nsim,
                                       includeRandomEffects=FALSE,
                                       maxMemoryMB=256){
  simulations <- matrix(NA_real_, nrow=object$nobs, ncol=nsim)
  familyIndex <- object$laplaceFamilyIndex
  randomFactors <- if(includeRandomEffects)
    .laplace_random_effect_sampling_factors(object, maxMemoryMB) else NULL
  for(simulationIndex in seq_len(nsim)){
    meanValue <- object$fitted.values
    if(includeRandomEffects){
      randomDraw <- .laplace_draw_random_effects(randomFactors)
      linearPredictor <- object$offset +
        as.numeric(object$fixedDesign %*% object$b) +
        as.numeric(object$randomDesign %*% randomDraw)
      meanValue <- .pql_apply(object$laplaceFamilies, familyIndex,
        "linkinv", linearPredictor)
    }
    for(familyPosition in seq_along(object$laplaceFamilies)){
      rows <- which(familyIndex == familyPosition)
      if(!length(rows)) next
      specification <- object$familySpecifications[[familyPosition]]
      currentMean <- meanValue[rows]
      trialCount <- object$priorWeights[rows]
      responseFamily <- specification$family
      simulations[rows, simulationIndex] <- switch(responseFamily,
        binomial=stats::rbinom(length(rows), size=trialCount,
          prob=currentMean) / trialCount,
        poisson=stats::rpois(length(rows), lambda=currentMean),
        ziPoisson={
          countMean <- currentMean / (1 -
            .laplace_zero_probability(specification, rows))
          inflated <- stats::rbinom(length(rows), size=1,
            prob=.laplace_zero_probability(specification, rows)) == 1
          draw <- stats::rpois(length(rows), lambda=countMean)
          draw[inflated] <- 0
          draw
        },
        nbinom2=stats::rnbinom(length(rows),
          size=exp(specification$logParameter), mu=currentMean),
        gaussian=stats::rnorm(length(rows), mean=currentMean,
          sd=sqrt(exp(specification$logParameter))),
        Gamma=stats::rgamma(length(rows),
          shape=exp(-specification$logParameter),
          scale=currentMean * exp(specification$logParameter)),
        stop("Unsupported family in native simulation: ", responseFamily, ".",
          call.=FALSE))
    }
  }
  simulations
}

.laplace_quantile_residuals <- function(object){
  response <- object$y
  meanValue <- object$fitted.values
  familyIndex <- object$laplaceFamilyIndex
  trials <- object$priorWeights
  randomizedValues <- stats::runif(length(response))
  probability <- numeric(length(response))
  for(familyPosition in seq_along(object$laplaceFamilies)){
    rows <- which(familyIndex == familyPosition)
    if(!length(rows)) next
    specification <- object$familySpecifications[[familyPosition]]
    currentMean <- meanValue[rows]
    currentResponse <- response[rows]
    currentTrials <- trials[rows]
    probability[rows] <- switch(specification$family,
      binomial={
        successes <- round(currentResponse * currentTrials)
        lower <- stats::pbinom(successes - 1, size=currentTrials,
          prob=currentMean)
        upper <- stats::pbinom(successes, size=currentTrials,
          prob=currentMean)
        lower + randomizedValues[rows] * (upper - lower)
      },
      poisson={
        lower <- stats::ppois(currentResponse - 1, lambda=currentMean)
        upper <- stats::ppois(currentResponse, lambda=currentMean)
        lower + randomizedValues[rows] * (upper - lower)
      },
      ziPoisson={
        probabilityZero <- .laplace_zero_probability(specification, rows)
        countMean <- currentMean / (1 - probabilityZero)
        lower <- probabilityZero + (1 - probabilityZero) *
          stats::ppois(currentResponse - 1, lambda=countMean)
        lower[currentResponse == 0] <- 0
        upper <- probabilityZero + (1 - probabilityZero) *
          stats::ppois(currentResponse, lambda=countMean)
        lower + randomizedValues[rows] * (upper - lower)
      },
      nbinom2={
        shape <- exp(specification$logParameter)
        lower <- stats::pnbinom(currentResponse - 1, size=shape, mu=currentMean)
        upper <- stats::pnbinom(currentResponse, size=shape, mu=currentMean)
        lower + randomizedValues[rows] * (upper - lower)
      },
      gaussian=stats::pnorm(currentResponse, mean=currentMean,
        sd=sqrt(exp(specification$logParameter))),
      Gamma=stats::pgamma(currentResponse,
        shape=exp(-specification$logParameter),
        scale=currentMean * exp(specification$logParameter)),
      stop("Unsupported family in native quantile residuals: ",
        specification$family, ".", call.=FALSE))
  }
  probability <- pmin(1 - .Machine$double.eps,
    pmax(.Machine$double.eps, probability))
  stats::qnorm(probability)
}

.laplace_population_variance <- function(randomDesign, structures,
                                         relationships, blockWidths){
  if(length(structures) != length(relationships) ||
     length(structures) != length(blockWidths) ||
     sum(blockWidths) != ncol(randomDesign)){
    stop("Population-marginal prediction could not align random-effect covariance blocks.",
         call.=FALSE)
  }
  variance <- numeric(nrow(randomDesign))
  columnStart <- 1L
  for(termIndex in seq_along(structures)){
    structure <- structures[[termIndex]]
    blockWidth <- blockWidths[termIndex]
    blockColumns <- columnStart + seq_len(blockWidth) - 1L
    block <- randomDesign[, blockColumns, drop=FALSE]
    relationshipPrecision <- relationships[[termIndex]]
    if(structure$dim != 1L ||
       nrow(relationshipPrecision) != blockWidth ||
       ncol(relationshipPrecision) != blockWidth ||
       !Matrix::isDiagonal(relationshipPrecision) ||
       any(abs(Matrix::diag(relationshipPrecision) - 1) > 1e-10)){
      stop("Population-marginal response means currently require independent random-intercept terms with identity relationship precision; use re.form=NA for other covariance structures.",
           call.=FALSE)
    }
    covariance <- evaluate_covstruct_cpp(structure, structure$par)
    if(!all(dim(covariance) == c(1L, 1L)) ||
       !is.finite(covariance[1L, 1L]) || covariance[1L, 1L] < 0){
      stop("Population-marginal prediction requires a finite non-negative random-intercept variance.",
           call.=FALSE)
    }
    variance <- variance + as.numeric(Matrix::rowSums(block^2)) *
      covariance[1L, 1L]
    columnStart <- columnStart + blockWidth
  }
  variance
}

.laplace_population_mean <- function(object, linearPredictor,
                                     familyIndex, randomVariance){
  if(length(linearPredictor) != length(randomVariance)){
    stop("Population-marginal predictor and random-effect variance have incompatible lengths.",
         call.=FALSE)
  }
  quadratureOrder <- 24L
  jacobiMatrix <- matrix(0, nrow=quadratureOrder, ncol=quadratureOrder)
  offDiagonal <- sqrt(seq_len(quadratureOrder - 1L))
  jacobiMatrix[cbind(seq_len(quadratureOrder - 1L),
    2:quadratureOrder)] <- offDiagonal
  jacobiMatrix[cbind(2:quadratureOrder,
    seq_len(quadratureOrder - 1L))] <- offDiagonal
  quadrature <- eigen(jacobiMatrix, symmetric=TRUE)
  quadratureWeights <- quadrature$vectors[1L, ]^2
  quadratureNodes <- quadrature$values
  result <- numeric(length(linearPredictor))
  for(familyPosition in seq_along(object$laplaceFamilies)){
    rows <- which(familyIndex == familyPosition)
    if(!length(rows)) next
    familyName <- object$familySpecifications[[familyPosition]]$family
    linkName <- object$familySpecifications[[familyPosition]]$link
    if(linkName == "log" && familyName %in% c("poisson", "ziPoisson", "nbinom2", "Gamma")){
      result[rows] <- exp(linearPredictor[rows] + randomVariance[rows] / 2)
      if(familyName == "ziPoisson"){
        result[rows] <- result[rows] *
          (1 - .laplace_zero_probability(object$familySpecifications[[familyPosition]], rows))
      }
    }else if(familyName == "gaussian"){
      result[rows] <- linearPredictor[rows]
    }else{
      for(rowPosition in rows){
        shiftedPredictors <- linearPredictor[rowPosition] +
          sqrt(randomVariance[rowPosition]) * quadratureNodes
        result[rowPosition] <- sum(quadratureWeights *
          object$laplaceFamilies[[familyPosition]]$linkinv(shiftedPredictors))
      }
    }
  }
  result
}

simulate.mmes.glmm.laplace <- function(object, nsim=1, seed=NULL,
                                       re.form=c("conditional", "population"),
                                       maxMemoryMB=NULL, ...){
  if(length(list(...))){
    stop("Native Laplace simulate() does not accept additional arguments.",
         call.=FALSE)
  }
  nsim <- .mme_positive(nsim, "nsim")
  if(length(nsim) != 1L || nsim != floor(nsim)){
    stop("nsim must be one positive integer.", call.=FALSE)
  }
  re.form <- match.arg(re.form)
  if(is.null(maxMemoryMB)){
    settings <- object$laplacePostPEVSettings
    maxMemoryMB <- if(is.null(settings$denseMemoryMB)) 256 else
      settings$denseMemoryMB
  }
  maxMemoryMB <- .mme_positive(maxMemoryMB, "maxMemoryMB")
  if(!is.null(seed)){
    if(length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
       seed != floor(seed)){
      stop("seed must be NULL or one integer.", call.=FALSE)
    }
    hadRandomState <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if(hadRandomState) oldRandomState <- get(".Random.seed", envir=.GlobalEnv)
    on.exit(if(hadRandomState){
      assign(".Random.seed", oldRandomState, envir=.GlobalEnv)
    }else if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)){
      rm(".Random.seed", envir=.GlobalEnv)
    }, add=TRUE)
    set.seed(seed)
  }
  outputMB <- as.double(object$nobs) * nsim * 8 / 1048576
  if(!is.finite(outputMB) || outputMB > maxMemoryMB){
    stop("Simulation output needs approximately ",
      format(round(outputMB, 1), scientific=FALSE),
      " MB, exceeding maxMemoryMB=", maxMemoryMB, ".", call.=FALSE)
  }
  simulationMatrix <- .laplace_simulate_response(object, nsim,
    includeRandomEffects=re.form == "population", maxMemoryMB=maxMemoryMB)
  result <- as.data.frame(simulationMatrix, check.names=FALSE)
  names(result) <- paste0("sim_", seq_len(nsim))
  rownames(result) <- rownames(object$fixedDesign)
  result
}
