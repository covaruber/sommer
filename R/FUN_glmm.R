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
  if(any(successes[ok] < 0) || any(total[ok] <= 0) || any(successes[ok] > total[ok])){
    stop("Counts must satisfy 0 <= successes <= trials and trials > 0.", call.=FALSE)
  }
  structure(as.numeric(successes / total), trials=as.numeric(total),
            class=c("binm", "numeric"))
}

"[.binm" <- function(x, i){
  structure(unclass(x)[i], trials=attr(x, "trials")[i], class=class(x))
}

familym <- function(by, ...){
  families <- list(...)
  if(!is.character(by) || length(by) != 1L){
    stop("by must be the name of the column that identifies the traits.", call.=FALSE)
  }
  if(length(families) < 1L || is.null(names(families)) || any(!nzchar(names(families)))){
    stop("Families must be named by the levels of the trait column, e.g. Yield = gaussian().",
         call.=FALSE)
  }
  if(!all(vapply(families, inherits, logical(1), what="family"))){
    stop("Every element of familym() must be a family object.", call.=FALSE)
  }
  structure(list(by=by, families=families), class="sommer_familym")
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
  defaults <- list(maxit=20L, tol=1e-5, estimateTheta=FALSE, warmStart=TRUE)
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
  families <- if(mixed) family$families else list(family)
  isNegBin <- !mixed && grepl("^Negative Binomial", family$family)
  estimateTheta <- isTRUE(control$estimateTheta)
  if(estimateTheta && !isNegBin){
    stop("pqlControl$estimateTheta requires family=MASS::negative.binomial(theta) (single family).",
         call.=FALSE)
  }
  theta <- if(isNegBin) get(".Theta", envir=environment(family$variance)) else NULL

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
                        .pqlInner=TRUE, .pqlFixedDispersion=fixedDispersion,
                .pqlWorkingPrecision=workingWeight,
                        .pqlBaseW=W, .pqlBaseFactor=baseFactor,
                        .pqlStart=covarianceStart),
                   mmesArgs)
    if(!is.null(random)) innerArgs$random <- random
    if(!is.null(rcov)) innerArgs$rcov <- rcov
    fit <- do.call(mmes, innerArgs)
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
  fit$pqlConverged <- nrow(monitor) < control$maxit
  fit$pqlControl <- control
  fit$.ldltCache <- NULL
  fit$.cholmodCache <- NULL
  attr(fit$covStruct, "ldltCache") <- NULL
  attr(fit$covStruct, "cholmodCache") <- NULL
  class(fit) <- c("mmes.glmm", class(fit))
  fit
}