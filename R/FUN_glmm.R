# Column of dataWork referenced via glm()'s non-standard evaluation of 'weights'.
utils::globalVariables(".sommer_pql_weights")

binm <- function(successes, failures=NULL, trials=NULL){
  if(is.null(failures) == is.null(trials)){
    stop("Provide exactly one of 'failures' or 'trials'.", call.=FALSE)
  }
  if(!is.numeric(successes)) stop("successes must be numeric.", call.=FALSE)
  total <- if(is.null(trials)) successes + failures else trials
  if(!is.numeric(total) || length(total) != length(successes)){
    stop("successes and failures/trials must be numeric vectors of the same length.",
         call.=FALSE)
  }
  ok <- !is.na(successes) & !is.na(total)
  if(any(!is.finite(successes[ok])) || any(!is.finite(total[ok])) ||
     any(successes[ok] != floor(successes[ok])) || any(total[ok] != floor(total[ok])) ||
     any(successes[ok] < 0) || any(total[ok] <= 0) || any(successes[ok] > total[ok])){
    stop("Counts must be finite integers satisfying 0 <= successes <= trials and trials > 0.", call.=FALSE)
  }
  structure(as.numeric(successes / total), trials=as.numeric(total),
            class=c("binm", "numeric"))
}

"[.binm" <- function(x, i){
  structure(unclass(x)[i], trials=attr(x, "trials")[i], class=class(x))
}

familym <- function(by, ..., trials=NULL){
  families <- list(...)
  if(!is.character(by) || length(by) != 1L || is.na(by) || !nzchar(by)){
    stop("by must be the name of the column that identifies the traits.", call.=FALSE)
  }
  if(length(families) < 1L || is.null(names(families)) || anyNA(names(families)) ||
     any(!nzchar(names(families))) || anyDuplicated(names(families))){
    stop("Families must be named by the levels of the trait column, e.g. Yield = gaussian().",
         call.=FALSE)
  }
  if(!all(vapply(families, inherits, logical(1), what="family"))){
    stop("Every element of familym() must be a family object.", call.=FALSE)
  }
  if(!is.null(trials) && (!is.character(trials) || length(trials) != 1L ||
     is.na(trials) || !nzchar(trials))){
    stop("familym() trials must be NULL or a column name.", call.=FALSE)
  }
  structure(list(by=by, families=families, trials=trials), class="sommer_familym")
}

.pql_fixed_dispersion <- function(fam){
  fam$family %in% c("binomial", "poisson") || grepl("^Negative Binomial", fam$family)
}

# Fix the residual variance of fixed-dispersion traits at one: the residual
# scale (variance of the first trait) and the variance ratios of the others.
.pql_fix_levels <- function(cs, levels, by){
  hit <- which(vapply(cs$factors, function(f) !is.null(f$ratio_index) &&
                        all(levels %in% f$levels), logical(1)))
  if(!length(hit)){
    stop("Mixed-family models need a residual structure with trait-specific variances over '", by,
         "', e.g. rcov = ~vsm(corgm(", by, ", variance=\"heterogeneous\"), ism(record)) or ",
         "~vsm(dsm(", by, "), ism(units)).", call.=FALSE)
  }
  f <- cs$factors[[hit[1L]]]
  if(!(f$levels[1L] %in% levels)){
    stop("The first level of '", by, "' must be a binomial/Poisson/negative-binomial trait ",
         "(its residual variance is the fixed reference); reorder the levels of '", by, "'.",
         call.=FALSE)
  }
  cs$par[1L] <- 0
  cs$free[1L] <- FALSE
  for(l in setdiff(match(levels, f$levels), 1L)){
    k <- f$par_start - 1L + f$ratio_index[l - 1L]
    cs$par[k] <- 0
    cs$free[k] <- FALSE
  }
  cs
}

# Apply a family component row-wise when rows belong to different families.
.pql_apply <- function(families, index, what, x, ...){
  out <- rep(NA_real_, length(x))
  extra <- list(...)
  for(g in seq_along(families)){
    rows <- which(index == g)
    if(!length(rows)) next
    args <- c(list(x[rows]), lapply(extra, function(e) e[rows]))
    out[rows] <- do.call(families[[g]][[what]], args)
  }
  out
}

.mmes_pql <- function(fixed, random, rcov, data, W, family, pqlControl,
                      mmesArgs){
  defaults <- .mme_pql_defaults
  if(!is.list(pqlControl) || is.null(names(pqlControl)) && length(pqlControl)){
    stop("pqlControl must be a named list.", call.=FALSE)
  }
  control <- utils::modifyList(defaults, pqlControl)
  if(length(control$warmStart) != 1L || !is.logical(control$warmStart) || is.na(control$warmStart)){
    stop("pqlControl$warmStart must be TRUE or FALSE.", call.=FALSE)
  }
  if(length(control$maxit) != 1L || !is.finite(control$maxit) || control$maxit < 1L ||
     control$maxit != as.integer(control$maxit)){
    stop("pqlControl$maxit must be a positive integer.", call.=FALSE)
  }
  if(length(control$tol) != 1L || !is.finite(control$tol) || control$tol <= 0){
    stop("pqlControl$tol must be a positive finite number.", call.=FALSE)
  }
  mixed <- inherits(family, "sommer_familym")
  if(mixed && !is.null(family$trials)){
    stop("familym(trials=...) currently requires native Laplace fitting; PQL trials use binm().", call.=FALSE)
  }
  families <- if(mixed) family$families else list(family)
  if(any(vapply(families, function(fam) identical(fam$family, "ziPoisson"),
      logical(1)))){
    stop("zipoisson() is currently supported only by native Laplace fitting.",
         call.=FALSE)
  }
  isNegBin <- !mixed && grepl("^Negative Binomial", family$family)
  estimateTheta <- isTRUE(control$estimateTheta)
  if(estimateTheta && !isNegBin){
    stop("pqlControl$estimateTheta requires family=MASS::negative.binomial(theta) (single family).",
         call.=FALSE)
  }
  theta <- if(isNegBin){
    if(!is.null(family$theta)) family$theta else
      get(".Theta", envir=environment(family$variance))
  }else NULL

  modelFrame <- stats::model.frame(fixed, data=data, na.action=stats::na.pass)
  response <- stats::model.response(modelFrame)
  if(!is.null(dim(response)) || !is.numeric(response)){
    stop("Non-Gaussian mmes() fits currently require one numeric response.", call.=FALSE)
  }
  n <- length(response)
  if(mixed){
    if(is.null(data) || !(family$by %in% names(data))){
      stop("The trait column '", family$by, "' of familym() must be in data.", call.=FALSE)
    }
    byValues <- as.character(data[[family$by]])
    famIndex <- match(byValues, names(families))
    unknown <- unique(byValues[is.na(famIndex) & !is.na(byValues)])
    if(length(unknown)){
      stop("Traits without a family in familym(): ", paste(unknown, collapse=", "), call.=FALSE)
    }
  }else{
    famIndex <- rep(1L, n)
  }
  priorWeights <- attr(response, "trials")
  response <- as.numeric(response)
  if(is.null(priorWeights)){
    priorWeights <- rep(1, n)
    isBinomialRow <- famIndex %in% which(vapply(families, function(f)
      identical(f$family, "binomial"), logical(1)))
    if(any(isBinomialRow & response %% 1 != 0 & !is.na(response))){
      warning("A binomial response with non-0/1 values but no number of trials was supplied; ",
              "use binm(successes, failures) so that proportions are weighted by their trials.",
              call.=FALSE)
    }
  }else if(length(priorWeights) != n){
    stop("The number of trials must have one value per observation.", call.=FALSE)
  }
  priorWeights[is.na(priorWeights)] <- 1
  offset <- stats::model.offset(modelFrame)
  if(is.null(offset)) offset <- rep(0, n)
  if(length(offset) != n){
    stop("The model offset must have one value per observation.", call.=FALSE)
  }
  if(is.null(data)){
    dataWork <- as.data.frame(modelFrame)
  }else{
    dataWork <- as.data.frame(data)
  }
  if(nrow(dataWork) != n){
    stop("The fixed formula and 'data' do not describe the same number of observations.",
         call.=FALSE)
  }

  fixedWork <- fixed
  fixedWork[[2L]] <- as.name(".sommer_pql_response")
  environment(fixedWork) <- environment(fixed)
  dataWork$.sommer_pql_response <- response
  dataWork$.sommer_pql_weights <- priorWeights
  if(mixed){
    # Standard IRLS starting means, family by family.
    mustart <- response
    for(g in seq_along(families)){
      rows <- which(famIndex == g)
      fg <- families[[g]]$family
      if(fg == "binomial"){
        mustart[rows] <- (priorWeights[rows]*response[rows] + 0.5) / (priorWeights[rows] + 1)
      }else if(fg == "poisson" || grepl("^Negative Binomial", fg)){
        mustart[rows] <- response[rows] + 0.1
      }
    }
    eta <- .pql_apply(families, famIndex, "linkfun", mustart)
  }else{
    initial <- stats::glm(fixedWork, data=dataWork, family=family,
                          weights=.sommer_pql_weights,
                          na.action=stats::na.exclude)
    eta <- as.numeric(stats::predict(initial, newdata=dataWork, type="link"))
  }

  if(!is.null(W)){
    if(nrow(W) != n || ncol(W) != n){
      stop("For non-Gaussian mmes() fits, W must have one row and column per original observation.",
           call.=FALSE)
    }
  }

  fixedByFamily <- vapply(families, .pql_fixed_dispersion, logical(1))
  fixedDispersion <- if(all(fixedByFamily)) TRUE else if(any(fixedByFamily) && mixed)
    list(by=family$by, levels=names(families)[fixedByFamily]) else FALSE
  monitor <- vector("list", control$maxit)
  fit <- NULL
  covarianceStart <- NULL
  baseFactor <- NULL
  previousDeviance <- Inf
  pqlConverged <- FALSE
  for(iteration in seq_len(control$maxit)){
    mu <- .pql_apply(families, famIndex, "linkinv", eta)
    derivative <- .pql_apply(families, famIndex, "mu.eta", eta)
    variance <- .pql_apply(families, famIndex, "variance", mu)
    # Links such as Gamma's inverse have negative mu.eta; only its square enters the weight.
    valid <- is.finite(response) & is.finite(offset) & is.finite(eta) & is.finite(mu) &
      is.finite(derivative) & is.finite(variance) & derivative != 0 & variance > 0 &
      is.finite(priorWeights) & priorWeights > 0
    workingResponse <- rep(NA_real_, n)
    workingResponse[valid] <- eta[valid] +
      (response[valid] - mu[valid]) / derivative[valid] - offset[valid]
    workingWeight <- rep(NA_real_, n)
    workingWeight[valid] <- priorWeights[valid] * derivative[valid]^2 / variance[valid]
    if(any(!is.finite(workingWeight[valid])) || any(workingWeight[valid] <= 0)){
      stop("IRLS produced invalid working weights.", call.=FALSE)
    }
    dataWork$.sommer_pql_response <- workingResponse

    innerArgs <- c(list(fixed=fixedWork, data=dataWork,
                        family=stats::gaussian(), pqlControl=list(),
                        .internal=list(pqlInner=TRUE, pqlFixedDispersion=fixedDispersion,
                          pqlWorkingPrecision=workingWeight,
                          pqlBaseW=W, pqlBaseFactor=baseFactor,
                          pqlStart=covarianceStart)),
                   mmesArgs)
    if(!is.null(random)) innerArgs$random <- random
    if(!is.null(rcov)) innerArgs$rcov <- rcov
    fit <- do.call(.mmes_fit, innerArgs)
    if(control$warmStart && !is.null(fit$covParWorking)){
      covarianceStart <- list(parameters=fit$covParWorking, covStruct=fit$covStruct,
                              included=fit$obsInfo$included, ldltCache=fit$.ldltCache,
                              cholmodCache=fit$.cholmodCache)
    }

    included <- fit$obsInfo$included
    if(is.null(baseFactor) && !is.null(W)){
      retainedW <- W[included, included, drop=FALSE]
      retainedW <- as(as(as(retainedW, "dMatrix"), "generalMatrix"),
                      "CsparseMatrix")
      baseFactor <- tryCatch(Matrix::chol(retainedW), error=function(e) NULL)
      if(is.null(baseFactor)){
        stop("PQL base W must be positive definite after observation filtering.",
             call.=FALSE)
      }
    }
    eta[included] <- offset[included] + as.numeric(fit$W %*% fit$bu)
    muIncluded <- .pql_apply(families, famIndex[included], "linkinv", eta[included])
    thetaChange <- 0
    if(estimateTheta){
      newTheta <- as.numeric(suppressWarnings(MASS::theta.ml(
        response[included], muIncluded, n=sum(priorWeights[included]),
        weights=priorWeights[included])))
      if(!is.finite(newTheta) || newTheta <= 0){
        stop("Estimation of the negative binomial theta failed.", call.=FALSE)
      }
      thetaChange <- abs(newTheta - theta) / theta
      theta <- newTheta
      family <- MASS::negative.binomial(theta, link=family$link)
      families <- list(family)
    }
    deviance <- sum(.pql_apply(families, famIndex[included], "dev.resids",
                               response[included], muIncluded, priorWeights[included]))
    monitor[[iteration]] <- c(iteration=iteration, deviance=deviance,
                              theta=if(isNegBin) theta else NA_real_,
                              warmStarted=as.integer(fit$pqlWarmStarted),
                              innerIterations=ncol(fit$monitor),
                              symbolicAnalyses=if(is.null(fit$engineDiagnostics$CsymbolicAnalyses))
                                NA_real_ else fit$engineDiagnostics$CsymbolicAnalyses)
    if(is.finite(previousDeviance) && thetaChange <= control$tol &&
       abs(deviance - previousDeviance) <= control$tol * (1 + abs(previousDeviance))){
      monitor <- monitor[seq_len(iteration)]
      pqlConverged <- TRUE
      break
    }
    previousDeviance <- deviance
  }

  monitor <- do.call(rbind, monitor)
  fit$workingResponse <- Matrix::Matrix(dataWork$.sommer_pql_response[fit$obsInfo$included],
                                        sparse=TRUE)
  fit$workingWeights <- workingWeight[fit$obsInfo$included]
  included <- fit$obsInfo$included
  fit$y <- Matrix::Matrix(response[included], sparse=TRUE)
  fit$priorWeights <- priorWeights[included]
  fit$family <- family
  fit$pqlFamilies <- families
  fit$pqlFamilyIndex <- famIndex[included]
  fit$offset <- offset[fit$obsInfo$included]
  fit$linear.predictors <- eta[fit$obsInfo$included]
  fit$linear.predictorsNoOffset <- fit$linear.predictors - fit$offset
  fit$fitted.values <- .pql_apply(families, fit$pqlFamilyIndex, "linkinv", fit$linear.predictors)
  fit$deviance <- unname(monitor[nrow(monitor), "deviance"])
  if(isNegBin){
    fit$theta.nb <- theta
    fit$theta.nb.se <- if(estimateTheta) attr(MASS::theta.ml(
      response[included], fit$fitted.values, n=sum(fit$priorWeights),
      weights=fit$priorWeights), "SE") else NA_real_
  }
  pearson <- fit$priorWeights * (response[included] - fit$fitted.values)^2 /
    .pql_apply(families, fit$pqlFamilyIndex, "variance", fit$fitted.values)
  nIncluded <- sum(included)
  fixedRank <- .fixed_rank(fit)
  fit$dispersion <- vapply(seq_along(families), function(g){
    if(fixedByFamily[g]) return(1)
    rows <- fit$pqlFamilyIndex == g
    sum(pearson[rows]) / max(1, sum(rows) - fixedRank * sum(rows) / nIncluded)
  }, numeric(1))
  if(mixed){
    # The trait-specific residual variances are estimated by REML in the
    # working model; report them when R is indexed by the trait levels.
    Rk <- tryCatch(covmatrix_mmes(fit, length(fit$covStruct), se=FALSE)$covariance,
                   error=function(e) NULL)
    if(!is.null(Rk) && all(names(families) %in% rownames(Rk))){
      est <- !fixedByFamily
      fit$dispersion[est] <- diag(Rk)[names(families)[est]]
    }
  }
  names(fit$dispersion) <- if(mixed) names(families) else NULL
  if(!mixed) fit$dispersion <- unname(fit$dispersion)
  # The working-model likelihood is not a likelihood of the GLMM.
  fit$llik[] <- NA_real_
  fit$AIC <- NA_real_
  fit$BIC <- NA_real_
  fit$pqlMonitor <- monitor
  fit$pqlConverged <- pqlConverged
  fit$pqlControl <- control
  fit$.ldltCache <- NULL
  fit$.cholmodCache <- NULL
  attr(fit$covStruct, "ldltCache") <- NULL
  attr(fit$covStruct, "cholmodCache") <- NULL
  class(fit) <- c("mmes.glmm", class(fit))
  fit
}

.mmes_laplace_curvature <- function(hessian){
  if(!nrow(hessian)){
    return(list(valid=TRUE, reason="no estimated parameters",
      eigenvalues=numeric(), covariance=matrix(numeric(), 0L, 0L)))
  }
  if(any(!is.finite(hessian))){
    return(list(valid=FALSE, reason="non-finite numerical Hessian"))
  }
  hessian <- (hessian + t(hessian)) / 2
  eigenvalues <- eigen(hessian, symmetric=TRUE, only.values=TRUE)$values
  if(!length(eigenvalues) ||
     min(eigenvalues) <= sqrt(.Machine$double.eps) * max(abs(eigenvalues))){
    return(list(valid=FALSE, reason="singular or indefinite numerical Hessian",
      eigenvalues=eigenvalues))
  }
  list(valid=TRUE, reason="positive-definite", eigenvalues=eigenvalues,
    covariance=chol2inv(chol(hessian)))
}

.mmes_laplace_boundary <- function(structures, reference){
  boundaryParameters <- character()
  threshold <- sqrt(.Machine$double.eps)
  for(termIndex in seq_along(structures)){
    descriptor <- structures[[termIndex]]
    covariance <- evaluate_covstruct_cpp(descriptor, descriptor$par)
    covarianceEigenvalues <- eigen(covariance, symmetric=TRUE,
      only.values=TRUE)$values
    if(min(covarianceEigenvalues) <= threshold * max(covarianceEigenvalues)){
      boundaryParameters <- c(boundaryParameters,
        paste0(termIndex, ":near-singular covariance"))
    }
    for(parameterIndex in which(descriptor$free)){
      coordinate <- descriptor$par[parameterIndex]
      kind <- .vc_kind(descriptor, parameterIndex)
      nearBoundary <- switch(kind$kind,
        scale=coordinate - reference[[termIndex]]$par[parameterIndex] <
          log(threshold),
        exp=coordinate - reference[[termIndex]]$par[parameterIndex] <
          log(threshold),
        tanh=1 - abs(tanh(coordinate)) < threshold,
        bounded_logit=min(stats::plogis(coordinate),
          stats::plogis(-coordinate)) < threshold,
        FALSE)
      if(nearBoundary){
        boundaryParameters <- c(boundaryParameters,
          paste0(termIndex, ":", descriptor$par_names[parameterIndex]))
      }
    }
  }
  boundaryParameters
}

.mmes_laplace_uncertainty <- function(object){
  if(!inherits(object, "mmes.glmm.laplace")){
    stop("Post-fit Laplace uncertainty requires a native Laplace fit.",
         call.=FALSE)
  }
  fixedCount <- length(object$b)
  betaSummary <- cbind(Estimate=object$b, `Std. Error`=NA_real_,
    `z value`=NA_real_)
  dimnames(betaSummary) <- list(names(object$b),
    c("Estimate", "Std. Error", "z value"))
  covarianceTable <- object$covParNative
  covarianceTable$StdError <- NA_real_
  covarianceTable$Zratio <- NA_real_
  covariance <- matrix(NA_real_, nrow=nrow(covarianceTable),
    ncol=nrow(covarianceTable),
    dimnames=list(paste(covarianceTable$term, covarianceTable$parameter, sep=":"),
      paste(covarianceTable$term, covarianceTable$parameter, sep=":")))
  familyParameters <- object$familyParameters
  zeroCoefficients <- object$zeroCoefficients
  if(!is.null(familyParameters)) familyParameters$StdError[] <- NA_real_
  unavailable <- function(reason){
    warning("Native Laplace standard errors are unavailable: ", reason, ".",
      call.=FALSE)
    list(valid=FALSE, reason=reason, fixedCovariance=NULL,
      betaSummary=betaSummary, covParNative=covarianceTable,
      covParCovariance=covariance, familyParameters=familyParameters,
      zeroCoefficients=zeroCoefficients)
  }
  if(!isTRUE(object$convergence)){
    return(unavailable("outer optimization did not converge"))
  }
  if(length(object$laplaceDiagnostics$boundaryParameters)){
    return(unavailable(paste0("near-boundary parameter(s): ",
      paste(object$laplaceDiagnostics$boundaryParameters, collapse=", "))))
  }
  layout <- object$laplaceParameterLayout
  if(is.null(layout) || is.null(object$laplaceNativeResponse) ||
     length(object$optimizer$par) != fixedCount + layout$covarianceCount +
       sum(vapply(object$familySpecifications, .laplace_family_count, integer(1)))){
    stop("Native Laplace uncertainty inputs are incomplete; refit the model before requesting standard errors.",
      call.=FALSE)
  }
  settings <- object$laplacePostPEVSettings
  denseMemoryMB <- if(is.null(settings$denseMemoryMB)) 1000 else
    settings$denseMemoryMB
  modeMaxit <- if(is.null(settings$modeMaxit)) 100L else settings$modeMaxit
  modeTol <- if(is.null(settings$modeTol)) 1e-8 else settings$modeTol
  context <- sommer_laplace_objective_create(object$fixedDesign,
    object$randomDesign, object$laplaceNativeResponse, object$offset, "",
    object$covStruct, object$randomRelationships, denseMemoryMB,
    familySpecs=object$familySpecifications,
    familyIndex=object$laplaceFamilyIndex, trials=object$priorWeights)
  objective <- function(parameters){
    evaluation <- sommer_laplace_objective_eval(context, parameters,
      maxModeIterations=modeMaxit, tolerance=modeTol)
    if(!isTRUE(evaluation$valid)) return(Inf)
    -evaluation$logLik
  }
  hessianFailure <- NULL
  hessian <- tryCatch(stats::optimHess(object$optimizer$par, objective),
    error=function(condition){
      if(!grepl("^non-finite finite-difference value", conditionMessage(condition))){
        stop(condition)
      }
      hessianFailure <<- conditionMessage(condition)
      NULL
    })
  curvature <- if(is.null(hessian)){
    list(valid=FALSE, reason=hessianFailure)
  }else .mmes_laplace_curvature(hessian)
  if(!curvature$valid) return(unavailable(curvature$reason))

  fixedCovariance <- curvature$covariance[
    seq_len(fixedCount), seq_len(fixedCount), drop=FALSE]
  dimnames(fixedCovariance) <- list(names(object$b), names(object$b))
  fixedSE <- sqrt(diag(fixedCovariance))
  betaSummary[, "Std. Error"] <- fixedSE
  betaSummary[, "z value"] <- object$b / fixedSE

  varianceTable <- .vc_param_table(object$covStruct, object$randomMetadata)
  freeCoordinatePositions <- integer(nrow(covarianceTable))
  for(termIndex in seq_along(layout$freePositions)){
    termRows <- which(varianceTable$structure == termIndex &
      varianceTable$free)
    if(length(termRows)){
      freeCoordinatePositions[termRows] <- fixedCount +
        layout$parameterStarts[termIndex] + seq_along(termRows)
    }
  }
  estimatedRows <- which(freeCoordinatePositions > 0L)
  if(length(estimatedRows)){
    coordinateJacobian <- vapply(estimatedRows, function(tableRow){
      structureIndex <- varianceTable$structure[tableRow]
      parameterIndex <- varianceTable$position[tableRow]
      descriptor <- object$covStruct[[structureIndex]]
      kind <- .vc_kind(descriptor, parameterIndex)
      .mmes_laplace_vc_derivative(covarianceTable$estimate[tableRow],
        kind$kind, kind$lower, kind$upper)
    }, numeric(1))
    coordinateCovariance <- curvature$covariance[
      freeCoordinatePositions[estimatedRows],
      freeCoordinatePositions[estimatedRows], drop=FALSE]
    naturalCovariance <- outer(coordinateJacobian, coordinateJacobian) *
      coordinateCovariance
    naturalVariance <- diag(naturalCovariance)
    if(any(!is.finite(naturalVariance)) || any(naturalVariance < 0)){
      stop("Natural-scale covariance-parameter uncertainty has an invalid variance; no standard error was returned.",
        call.=FALSE)
    }
    covariance[estimatedRows, estimatedRows] <- naturalCovariance
    covarianceTable$StdError[estimatedRows] <-
      sqrt(naturalVariance)
    positiveSE <- estimatedRows[is.finite(covarianceTable$StdError[estimatedRows]) &
      covarianceTable$StdError[estimatedRows] > 0]
    covarianceTable$Zratio[positiveSE] <- covarianceTable$estimate[positiveSE] /
      covarianceTable$StdError[positiveSE]
  }

  if(!is.null(familyParameters) && nrow(familyParameters)){
    for(familyRow in seq_len(nrow(familyParameters))){
      if(!familyParameters$estimated[familyRow]) next
      familyPosition <- match(familyParameters$trait[familyRow],
        names(object$familySpecifications))
      counts <- vapply(object$familySpecifications, .laplace_family_count, integer(1))
      coordinate <- fixedCount + layout$covarianceCount +
        sum(head(counts, familyPosition - 1L)) + 1L
      parameterDerivative <- .laplace_family_parameter_derivative(
        object$familySpecifications[[familyPosition]])
      familyParameters$StdError[familyRow] <- parameterDerivative *
        sqrt(curvature$covariance[coordinate, coordinate])
    }
  }
  if(!is.null(zeroCoefficients)){
    counts <- vapply(object$familySpecifications, .laplace_family_count, integer(1))
    for(position in seq_along(object$familySpecifications)){
      specification <- object$familySpecifications[[position]]
      if(is.null(specification$zeroDesign)) next
      coordinates <- fixedCount + layout$covarianceCount +
        sum(head(counts, position - 1L)) + seq_len(counts[position])
      rows <- zeroCoefficients$trait == names(object$familySpecifications)[position]
      zeroCoefficients$StdError[rows] <- sqrt(diag(curvature$covariance)[coordinates])
    }
  }
  list(valid=TRUE, reason="positive-definite", fixedCovariance=fixedCovariance,
    betaSummary=betaSummary, covParNative=covarianceTable,
    covParCovariance=covariance, familyParameters=familyParameters,
    zeroCoefficients=zeroCoefficients,
    hessian=hessian, parameterCovariance=curvature$covariance)
}

.mmes_laplace_vc_derivative <- function(value, kind, lower, upper){
  switch(kind,
    scale=value,
    exp=value,
    identity=1,
    tanh=1 - value^2,
    bounded_logit={
      proportion <- (value - lower) / (upper - lower)
      (upper - lower) * proportion * (1 - proportion)
    },
    1
  )
}

.mmes_prepare_model <- function(fixed, random, rcov, data, family,
                                naMethodY, mmesArgs){
  preparationArguments <- list(fixed=fixed, family=family,
    naMethodY=naMethodY, naMethodX=mmesArgs$naMethodX,
    naMethodRandom=mmesArgs$naMethodRandom, naMethodR=mmesArgs$naMethodR,
    contrasts=mmesArgs$contrasts, engine=mmesArgs$engine,
    dateWarning=FALSE, verbose=FALSE, returnParam=TRUE, .prepareOnly=TRUE)
  if(!is.null(random)) preparationArguments$random <- random
  if(!is.null(rcov)) preparationArguments$rcov <- rcov
  if(!is.null(data)) preparationArguments$data <- data
  do.call(.mmes_fit, preparationArguments)
}

.mmes_laplace <- function(fixed, random, rcov, data, W, weights, family,
                          vcc, REML, maxit, modeMaxit, modeTol,
                          familyParameters=list(),
                          naMethodY, henderson,
                          solveOnly, computeCi, mmesArgs){
  if(isTRUE(REML)){
    stop("Laplace GLMMs use maximum likelihood; set task$objective='ml'.",
         call.=FALSE)
  }
  if(!isTRUE(henderson)){
    stop("The native Laplace module currently uses the Henderson formulation.",
         call.=FALSE)
  }
  if(isTRUE(solveOnly)){
    stop("Native Laplace fitting does not support solve mode.", call.=FALSE)
  }
  inverseMode <- as.integer(computeCi)
  if(inverseMode != 0L){
    stop("Native Laplace fits defer PEV calculations until after fitting; use postPEV() on the fitted model.",
         call.=FALSE)
  }
  panelSize <- as.integer(mmesArgs$engine$experimental$panelCi2)
  denseMemoryMB <- as.double(mmesArgs$engine$denseMemoryMB)
  familyModel <- .laplace_family_specs(family, familyParameters)
  families <- familyModel$families
  specifications <- familyModel$specifications
  if(!is.null(W) || !is.null(weights) || !is.null(vcc)){
    stop("Native Laplace fitting does not yet support W, weights formulas, or vcc.",
         call.=FALSE)
  }
  if(!is.null(rcov) &&
     !identical(attr(stats::terms(rcov), "term.labels"), "units")){
    stop("Native Laplace fitting uses the family likelihood for independent observation errors; structured rcov terms are not yet supported.",
         call.=FALSE)
  }
  if(is.null(random)){
    stop("Native Laplace fitting requires at least one random-effect term.",
         call.=FALSE)
  }
  if(!identical(naMethodY, "exclude") && !identical(naMethodY, "omit")){
    stop("Native Laplace fitting currently supports naMethodY='exclude' or 'omit'.",
         call.=FALSE)
  }
  if(identical(mmesArgs$engine$updateMethod, "conditional-fa") ||
     identical(mmesArgs$engine$conditionalEquations, TRUE) ||
     !identical(mmesArgs$engine$factorScoreAugmentation, "none") ||
     !identical(mmesArgs$engine$compression, "off") ||
     isTRUE(mmesArgs$engine$gpu) ||
     !identical(mmesArgs$engine$solver, "auto") ||
     !identical(mmesArgs$engine$acceleration, "none") ||
     !identical(mmesArgs$engine$updateMethod, "ai") ||
     !identical(mmesArgs$engine$fixedRank, "error") ||
     !identical(mmesArgs$engine$covarianceInformation, "full") ||
     !isTRUE(mmesArgs$engine$retainCoefficientMatrix) ||
     any(vapply(mmesArgs$engine$experimental, function(mode)
       identical(mode, "force"), logical(1)))){
    stop("Native Laplace fitting does not use non-default Gaussian solver, acceleration, update, fixed-rank, conditional FA, factor-score, compression, GPU, or forced experimental controls.",
         call.=FALSE)
  }

  prepared <- .mmes_prepare_model(fixed, random, rcov, data, family,
    naMethodY, mmesArgs)
  randomCount <- length(prepared$rtermss)
  if(!randomCount){
    stop("Native Laplace fitting requires at least one random-effect term.",
         call.=FALSE)
  }
  fixedDesign <- prepared$X
  randomDesign <- Matrix::Matrix(do.call(cbind, prepared$Z), sparse=TRUE)
  responseFrame <- prepared$modelFrame
  response <- stats::model.response(responseFrame)
  included <- prepared$obsInfo$included
  familyIndex <- .laplace_family_index(family, data, nrow(responseFrame))[included]
  if(anyNA(familyIndex)){
    stop("Retained observations must have a declared familym trait.", call.=FALSE)
  }
  if(any(tabulate(familyIndex, nbins=length(families)) == 0L)){
    stop("Each declared family must have retained observations.", call.=FALSE)
  }
  specifications <- .laplace_zero_designs(families, specifications, data,
    included, familyIndex)
  retainedFrame <- droplevels(responseFrame[included, , drop=FALSE])
  fixedTerms <- stats::delete.response(attr(responseFrame, "terms"))
  fixedXlevels <- stats::.getXlevels(fixedTerms, retainedFrame)
  fittedContrasts <- attr(fixedDesign, "contrasts")
  predictionColumns <- colnames(stats::model.matrix(fixedTerms, retainedFrame,
    contrasts.arg=fittedContrasts))
  fixedAssignment <- attr(stats::model.matrix(fixedTerms, retainedFrame,
    contrasts.arg=fittedContrasts), "assign")
  declaredTrials <- if(familyModel$mixed && !is.null(family$trials))
    data[[family$trials]][included] else NULL
  responseData <- .laplace_response(response, included, specifications, familyIndex,
    declaredTrials=declaredTrials)
  response <- responseData$response
  trials <- responseData$trials
  offset <- stats::model.offset(responseFrame)
  offset <- if(is.null(offset)) numeric(length(response)) else offset[included]
  if(length(response) != nrow(fixedDesign) || length(offset) != length(response)){
    stop("Prepared response, offset, and fixed design have incompatible dimensions.",
         call.=FALSE)
  }
  if(any(!is.finite(response)) || any(!is.finite(offset)) ||
     any(!is.finite(fixedDesign)) || any(!is.finite(randomDesign))){
    stop("Native Laplace fitting requires finite retained responses, offsets, and designs; use exclude/omit or fail missing-data policies.",
      call.=FALSE)
  }
  if(ncol(fixedDesign) &&
     Matrix::rankMatrix(fixedDesign)[1L] != ncol(fixedDesign)){
    stop("Native Laplace fitting requires a full-rank fixed-effects design.",
         call.=FALSE)
  }

  randomStructures <- prepared$covStruct[seq_len(randomCount)]
  freePositions <- lapply(randomStructures, function(structure)
    which(structure$free))
  covarianceParameterCount <- sum(vapply(freePositions, length, integer(1)))
  parameterCount <- ncol(fixedDesign)
  parameterStarts <- cumsum(c(0L,
    head(vapply(freePositions, length, integer(1)), -1L)))
  initialGLMs <- lapply(seq_along(families), function(familyPosition){
    rows <- familyIndex == familyPosition
    initialFamily <- if(specifications[[familyPosition]]$family == "ziPoisson")
      stats::poisson() else families[[familyPosition]]
    stats::glm.fit(as.matrix(fixedDesign[rows, , drop=FALSE]), response[rows],
      weights=trials[rows], family=initialFamily, offset=offset[rows])
  })
  if(length(families) == 1L){
    initialFixed <- initialGLMs[[1L]]$coefficients
  }else{
    initialPredictor <- numeric(length(response))
    for(familyPosition in seq_along(families)){
      rows <- familyIndex == familyPosition
      initialPredictor[rows] <- initialGLMs[[familyPosition]]$linear.predictors
    }
    initialFixed <- stats::lm.fit(as.matrix(fixedDesign), initialPredictor - offset)$coefficients
  }
  if(any(!is.finite(initialFixed))){
    stop("Native Laplace initialization produced non-finite fixed effects; check the design and separation.",
      call.=FALSE)
  }
  initialCovariance <- unlist(lapply(seq_len(randomCount), function(termIndex)
    randomStructures[[termIndex]]$par[freePositions[[termIndex]]]),
    use.names=FALSE)
  initialParameters <- c(initialFixed, initialCovariance)
  for(familyPosition in seq_along(specifications)){
    specification <- specifications[[familyPosition]]
    if(specification$free && !specification$supplied &&
       specification$family %in% c("gaussian", "Gamma")){
      initialGLM <- initialGLMs[[familyPosition]]
      rows <- familyIndex == familyPosition
      initialDispersion <- sum(families[[familyPosition]]$dev.resids(
        response[rows], initialGLM$fitted.values, trials[rows])) /
        max(1L, sum(rows) - initialGLM$rank)
      if(!is.finite(initialDispersion) || initialDispersion <= 0){
        stop("Native Laplace dispersion initialization is non-positive; supply a positive starting dispersion and check the data.", call.=FALSE)
      }
      specifications[[familyPosition]]$logParameter <- log(initialDispersion)
    }
  }
  freeFamilies <- which(vapply(specifications, function(specification) specification$free, logical(1)))
  initialFamilyParameters <- unlist(lapply(specifications[freeFamilies],
    function(specification) if(is.null(specification$zeroDesign))
      specification$logParameter else specification$zeroCoefficients),
    use.names=FALSE)
  initialParameters <- c(initialParameters, initialFamilyParameters)
  nativeContext <- sommer_laplace_objective_create(fixedDesign, randomDesign,
    responseData$nativeResponse, offset, "", randomStructures, prepared$Ai,
    mmesArgs$engine$denseMemoryMB, familySpecs=specifications,
    familyIndex=familyIndex, trials=trials)

  evaluateParameters <- function(parameters, gradient=FALSE){
    if(any(!is.finite(parameters))){
      return(list(valid=FALSE, reason="non-finite parameter proposal"))
    }
    structures <- randomStructures
    for(termIndex in seq_len(randomCount)){
      freeIndex <- freePositions[[termIndex]]
      if(length(freeIndex)){
        valueIndex <- parameterCount + parameterStarts[termIndex] +
          seq_along(freeIndex)
        structures[[termIndex]]$par[freeIndex] <- parameters[valueIndex]
      }
    }
    evaluation <- sommer_laplace_objective_eval(nativeContext, parameters,
      maxModeIterations=modeMaxit, tolerance=modeTol, gradient=gradient)
    evaluation$covStruct <- structures
    evaluation
  }
  initialEvaluation <- evaluateParameters(initialParameters)
  if(!isTRUE(initialEvaluation$valid)){
    stop("Native Laplace starting parameters are invalid: ",
         initialEvaluation$reason, ".",
         call.=FALSE)
  }
  rejectedEvaluations <- integer()
  objective <- function(parameters){
    evaluation <- evaluateParameters(parameters)
    if(!isTRUE(evaluation$valid)){
      rejectionReason <- evaluation$reason
      if(!rejectionReason %in% names(rejectedEvaluations)){
        rejectedEvaluations[rejectionReason] <<- 0L
      }
      rejectedEvaluations[rejectionReason] <<-
        rejectedEvaluations[rejectionReason] + 1L
      return(Inf)
    }
    -evaluation$logLik
  }
  objectiveGradient <- function(parameters){
    evaluation <- evaluateParameters(parameters, gradient=TRUE)
    if(!isTRUE(evaluation$valid)){
      stop("Native Laplace gradient requested at invalid parameters: ",
        evaluation$reason, ".", call.=FALSE)
    }
    -as.numeric(evaluation$gradient)
  }
  optimization <- if(length(initialParameters)){
    stats::nlminb(start=initialParameters, objective=objective,
      gradient=objectiveGradient,
      control=list(iter.max=maxit, eval.max=max(1000, 10 * maxit)))
  }else list(par=numeric(), objective=-initialEvaluation$logLik,
    convergence=0L, message="no estimated parameters")
  finalEvaluation <- evaluateParameters(optimization$par, gradient=TRUE)
  if(!isTRUE(finalEvaluation$valid)){
    stop("Native Laplace optimization ended at invalid parameters: ",
         finalEvaluation$reason, ".",
         call.=FALSE)
  }
  converged <- identical(optimization$convergence, 0L)
  if(!converged){
    warning("Native Laplace optimizer did not converge: ",
      if(is.null(optimization$message)) paste("code", optimization$convergence)
      else optimization$message, call.=FALSE)
  }

  fixedEffects <- optimization$par[seq_len(parameterCount)]
  boundaryParameters <- .mmes_laplace_boundary(finalEvaluation$covStruct,
    randomStructures)
  fittedSpecifications <- specifications
  familyCounts <- vapply(specifications, .laplace_family_count, integer(1))
  familyStarts <- cumsum(c(0L, head(familyCounts, -1L)))
  for(freePosition in seq_along(freeFamilies)){
    familyPosition <- freeFamilies[freePosition]
    positions <- parameterCount + covarianceParameterCount +
      familyStarts[familyPosition] + seq_len(familyCounts[familyPosition])
    if(!is.null(specifications[[familyPosition]]$zeroDesign)){
      fittedSpecifications[[familyPosition]]$zeroCoefficients[] <- optimization$par[positions]
      probabilities <- .laplace_zero_probability(fittedSpecifications[[familyPosition]],
        which(familyIndex == familyPosition))
      if(any(probabilities < sqrt(.Machine$double.eps) |
             probabilities > 1 - sqrt(.Machine$double.eps))){
        boundaryParameters <- c(boundaryParameters,
          paste0(names(specifications)[familyPosition], ":zeroInflation"))
      }
      next
    }
    coordinate <- optimization$par[positions]
    fittedSpecifications[[familyPosition]]$logParameter <- coordinate
    startCoordinate <- specifications[[familyPosition]]$logParameter
    specification <- fittedSpecifications[[familyPosition]]
    atBoundary <- if(specification$family == "nbinom2"){
      coordinate - startCoordinate > -log(sqrt(.Machine$double.eps)) ||
        coordinate - startCoordinate < log(sqrt(.Machine$double.eps))
    }else if(specification$family == "ziPoisson"){
      probability <- stats::plogis(coordinate)
      probability < sqrt(.Machine$double.eps) ||
        probability > 1 - sqrt(.Machine$double.eps)
    }else coordinate - startCoordinate < log(sqrt(.Machine$double.eps))
    if(atBoundary){
      boundaryParameters <- c(boundaryParameters,
        paste0(names(specifications)[familyPosition], ":",
          specification$parameter))
    }
  }
  fittedFamilies <- .laplace_update_families(families, fittedSpecifications, familyIndex)
  names(fittedFamilies) <- names(families)
  dispersion <- vapply(fittedSpecifications, function(specification)
    if(specification$family %in% c("gaussian", "Gamma")) exp(specification$logParameter)
    else 1, numeric(1))
  fixedStandardErrors <- rep(NA_real_, parameterCount)
  binaryRows <- familyIndex %in% which(vapply(specifications, function(specification)
    specification$family == "binomial", logical(1)))
  extremeMeans <- any(binaryRows & (finalEvaluation$fitted < .Machine$double.eps |
      finalEvaluation$fitted > 1 - .Machine$double.eps))
  if(extremeMeans){
    warning("Native Laplace fitted probabilities are numerically 0 or 1; check separation and identifiability.",
      call.=FALSE)
  }
  names(fixedEffects) <- colnames(fixedDesign)
  betaSummary <- cbind(
    Estimate=fixedEffects,
    `Std. Error`=fixedStandardErrors,
    `z value`=fixedEffects / fixedStandardErrors
  )
  colnames(betaSummary) <- c("Estimate", "Std. Error", "z value")
  rownames(betaSummary) <- colnames(fixedDesign)
  fixedCovariance <- NULL
  randomVariables <- unique(all.vars(random))
  randomFactorLevels <- list()
  if(!is.null(data)){
    randomFactorLevels <- lapply(randomVariables, function(variable){
      value <- data[[variable]]
      if(is.factor(value)) levels(value) else if(is.character(value))
        unique(value[!is.na(value)]) else NULL
    })
    names(randomFactorLevels) <- randomVariables
    randomFactorLevels <- Filter(Negate(is.null), randomFactorLevels)
  }
  fittedStructures <- finalEvaluation$covStruct
  varianceTable <- .vc_param_table(fittedStructures, prepared$rtermss)
  covarianceTable <- data.frame(
    term=varianceTable$term, parameter=varianceTable$parameter,
    estimate=varianceTable$start, free=varianceTable$free,
    StdError=NA_real_, Zratio=NA_real_,
    status=rep("Laplace estimate", nrow(varianceTable)),
    stringsAsFactors=FALSE
  )
  covarianceCoordinateCovariance <- matrix(NA_real_,
    nrow=nrow(covarianceTable), ncol=nrow(covarianceTable),
    dimnames=list(paste(covarianceTable$term, covarianceTable$parameter, sep=":"),
      paste(covarianceTable$term, covarianceTable$parameter, sep=":")))
  covarianceParameters <- lapply(fittedStructures, function(structure){
    kinds <- lapply(seq_along(structure$par), function(position)
      .vc_kind(structure, position))
    vapply(seq_along(structure$par), function(position)
      .vc_to_natural(structure$par[position], kinds[[position]]$kind,
        kinds[[position]]$lower, kinds[[position]]$upper), numeric(1))
  })
  names(covarianceParameters) <- prepared$rtermss
  fittedValues <- as.numeric(finalEvaluation$fitted)
  linearPredictors <- as.numeric(finalEvaluation$linearPredictors)
  logLikelihood <- as.numeric(finalEvaluation$logLik)
  parameterTotal <- parameterCount + covarianceParameterCount + sum(familyCounts)
  familyParameterTable <- do.call(rbind, lapply(seq_along(fittedSpecifications),
    function(familyPosition){
      specification <- fittedSpecifications[[familyPosition]]
      if(is.null(specification$parameter) || !is.null(specification$zeroDesign)) return(NULL)
      estimate <- .laplace_family_parameter_value(specification)
      data.frame(trait=names(families)[familyPosition], parameter=specification$parameter,
        estimate=estimate, StdError=NA_real_,
        estimated=specification$free)
    }))
  fittedObject <- list(
    call=match.call(),
    family=if(familyModel$mixed) structure(list(by=family$by,
      families=fittedFamilies, trials=family$trials), class="sommer_familym") else fittedFamilies[[1L]],
    b=fixedEffects,
    betaSummary=betaSummary, fitted.values=fittedValues,
    fixedCovariance=fixedCovariance,
    linear.predictors=linearPredictors, y=response, offset=offset,
    AIC=-2 * logLikelihood + 2 * parameterTotal,
    BIC=-2 * logLikelihood + log(length(response)) * parameterTotal,
    llik=logLikelihood,
    deviance=sum(.laplace_deviance(response, fittedValues, fittedFamilies,
      familyIndex, trials, dispersion)),
    laplaceFamilies=fittedFamilies, laplaceFamilyIndex=familyIndex,
    familySpecifications=fittedSpecifications,
    familyParameters=familyParameterTable,
    zeroCoefficients=.laplace_zero_table(fittedSpecifications), priorWeights=trials,
    dispersion=if(familyModel$mixed) dispersion else unname(dispersion),
    theta.nb=if(!familyModel$mixed && fittedSpecifications[[1L]]$family == "nbinom2")
      exp(fittedSpecifications[[1L]]$logParameter) else NULL,
    theta.nb.se=if(!familyModel$mixed && fittedSpecifications[[1L]]$family == "nbinom2")
      familyParameterTable$StdError else NULL,
    convergence=converged, covParNative=covarianceTable,
    covParCovariance=covarianceCoordinateCovariance,
    engine="native Laplace", REML=FALSE, nobs=length(response),
    npar=parameterTotal,
    fixedFormula=fixed, fixedTerms=fixedTerms, fixedXlevels=fixedXlevels,
    randomFormula=random, predictionColumns=predictionColumns,
    fixedAssignment=fixedAssignment,
    randomBlockColumns=lapply(prepared$Z, colnames),
    randomBlockWidths=vapply(seq_len(randomCount), function(termIndex)
      sum(vapply(prepared$Z[prepared$Zind == termIndex], ncol, integer(1))),
      integer(1)),
    randomDesignBlockWidths=vapply(prepared$Z, ncol, integer(1)),
    randomFactorLevels=randomFactorLevels,
    predictionSettings=list(rcov=rcov, family=family,
      naMethodY=naMethodY, mmesArgs=mmesArgs),
    fixedDesign=fixedDesign, randomDesign=randomDesign,
    randomMode=as.numeric(finalEvaluation$mode),
    laplaceNativeResponse=responseData$nativeResponse,
    laplaceModeHessian=finalEvaluation$modeHessian,
    laplaceObservationCurvature=as.numeric(finalEvaluation$observationCurvature),
    randomMetadata=prepared$rtermss, covStruct=fittedStructures,
    randomRelationships=prepared$Ai,
    covPar=covarianceParameters, u=matrix(as.numeric(finalEvaluation$mode), ncol=1L),
    contrasts=fittedContrasts, optimizer=optimization,
    laplaceParameterLayout=list(fixedCount=parameterCount,
      covarianceCount=covarianceParameterCount, freePositions=freePositions,
      parameterStarts=parameterStarts, freeFamilies=freeFamilies),
    obsInfo=prepared$obsInfo,
    laplaceDiagnostics=list(outerConverged=converged, modeConverged=TRUE,
      modeScore=finalEvaluation$modeScore,
      curvatureValid=NA, curvatureReason="not calculated; use vcov() or summary()",
      uncertaintyValid=FALSE, uncertaintyComputed=FALSE,
      boundary=length(boundaryParameters) > 0L,
      boundaryParameters=boundaryParameters,
      initialGLMConverged=all(vapply(initialGLMs, function(initialGLM)
        isTRUE(initialGLM$converged), logical(1))),
      gradient=as.numeric(finalEvaluation$gradient),
      nativeWork=finalEvaluation$work,
      extremeProbabilities=extremeMeans, rejectedEvaluations=rejectedEvaluations),
    modeIterations=finalEvaluation$iterations,
    logDetHessian=finalEvaluation$logDetHessian,
    executionPlan=list(glmm=list(selected="Native Laplace"))
  )
  class(fittedObject) <- c("mmes.glmm.laplace", "mmes.glmm", "mmes")
  fittedObject$Ci <- NULL
  fittedObject$VarU <- NULL
  fittedObject$uPevList <- NULL
  fittedObject$uPevMode <- 0L
  fittedObject$laplacePostPEVSettings <- list(inverse=inverseMode,
    panelSize=panelSize, denseMemoryMB=denseMemoryMB,
    modeMaxit=modeMaxit, modeTol=modeTol)
  fittedObject
}