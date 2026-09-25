.mmes_pql <- function(fixed, random, rcov, data, W, family, pqlControl,
                      mmesArgs){
  defaults <- list(maxit=20L, tol=1e-5)
  if(!is.list(pqlControl) || is.null(names(pqlControl)) && length(pqlControl)){
    stop("pqlControl must be a named list.", call.=FALSE)
  }
  control <- utils::modifyList(defaults, pqlControl)
  if(length(control$maxit) != 1L || !is.finite(control$maxit) || control$maxit < 1L ||
     control$maxit != as.integer(control$maxit)){
    stop("pqlControl$maxit must be a positive integer.", call.=FALSE)
  }
  if(length(control$tol) != 1L || !is.finite(control$tol) || control$tol <= 0){
    stop("pqlControl$tol must be a positive finite number.", call.=FALSE)
  }

  modelFrame <- stats::model.frame(fixed, data=data, na.action=stats::na.pass)
  response <- stats::model.response(modelFrame)
  if(!is.null(dim(response)) || !is.numeric(response)){
    stop("Non-Gaussian mmes() fits currently require one numeric response.", call.=FALSE)
  }
  n <- length(response)
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
  initial <- stats::glm(fixedWork, data=dataWork, family=family,
                        na.action=stats::na.exclude)
  eta <- as.numeric(stats::predict(initial, newdata=dataWork, type="link"))

  if(!is.null(W)){
    if(nrow(W) != n || ncol(W) != n){
      stop("For non-Gaussian mmes() fits, W must have one row and column per original observation.",
           call.=FALSE)
    }
  }

  fixedDispersion <- family$family %in% c("binomial", "poisson")
  monitor <- vector("list", control$maxit)
  fit <- NULL
  baseFactor <- NULL
  previousDeviance <- Inf
  for(iteration in seq_len(control$maxit)){
    mu <- family$linkinv(eta)
    derivative <- family$mu.eta(eta)
    variance <- family$variance(mu)
    valid <- is.finite(response) & is.finite(offset) & is.finite(eta) & is.finite(mu) &
      is.finite(derivative) & is.finite(variance) & derivative > 0 & variance > 0
    workingResponse <- rep(NA_real_, n)
    workingResponse[valid] <- eta[valid] +
      (response[valid] - mu[valid]) / derivative[valid] - offset[valid]
    workingWeight <- rep(NA_real_, n)
    workingWeight[valid] <- derivative[valid]^2 / variance[valid]
    if(any(!is.finite(workingWeight[valid])) || any(workingWeight[valid] <= 0)){
      stop("IRLS produced invalid working weights.", call.=FALSE)
    }
    dataWork$.sommer_pql_response <- workingResponse

    innerArgs <- c(list(fixed=fixedWork, data=dataWork,
                        family=stats::gaussian(), pqlControl=list(),
                        .pqlInner=TRUE, .pqlFixedDispersion=fixedDispersion,
                .pqlWorkingPrecision=workingWeight,
                        .pqlBaseW=W, .pqlBaseFactor=baseFactor),
                   mmesArgs)
    if(!is.null(random)) innerArgs$random <- random
    if(!is.null(rcov)) innerArgs$rcov <- rcov
    fit <- do.call(mmes, innerArgs)

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
    muIncluded <- family$linkinv(eta[included])
    deviance <- sum(family$dev.resids(response[included], muIncluded,
                      rep(1, sum(included))))
    monitor[[iteration]] <- c(iteration=iteration, deviance=deviance)
    if(is.finite(previousDeviance) &&
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
  fit$y <- Matrix::Matrix(response[fit$obsInfo$included], sparse=TRUE)
  fit$family <- family
  fit$offset <- offset[fit$obsInfo$included]
  fit$linear.predictors <- eta[fit$obsInfo$included]
  fit$linear.predictorsNoOffset <- fit$linear.predictors - fit$offset
  fit$fitted.values <- family$linkinv(fit$linear.predictors)
  fit$deviance <- monitor[nrow(monitor), "deviance"]
  fit$pqlMonitor <- monitor
  fit$pqlConverged <- nrow(monitor) < control$maxit
  fit$pqlControl <- control
  class(fit) <- c("mmes.glmm", class(fit))
  fit
}