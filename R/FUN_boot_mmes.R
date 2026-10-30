.boot_mmes_parameters <- function(object, parameters=object$optimizer$par){
  layout <- object$laplaceParameterLayout
  result <- stats::setNames(parameters[seq_len(layout$fixedCount)],
    if(layout$fixedCount) paste0("count:", names(object$b)) else character())
  for(term in seq_along(layout$freePositions)){
    positions <- layout$freePositions[[term]]
    coordinates <- layout$fixedCount + layout$parameterStarts[term] +
      seq_along(positions)
    if(length(positions)){
      values <- parameters[coordinates]
      names(values) <- paste0("covariance:", term, ":",
        object$covStruct[[term]]$par_names[positions])
      result <- c(result, values)
    }
  }
  specifications <- .laplace_proposal_specs(object$familySpecifications,
    parameters, layout)
  for(position in seq_along(specifications)){
    specification <- specifications[[position]]
    if(!specification$free) next
    if(specification$family == "ordinal"){
      values <- .laplace_ordinal_coordinates(specification$thresholds)
      names(values) <- paste0("family:", names(specifications)[position], ":",
        c("threshold1", if(length(values) > 1L) paste0("logGap", 2:length(values))))
      result <- c(result, values)
      next
    }
    values <- if(is.null(specification$zeroDesign)) specification$logParameter else
      specification$zeroCoefficients
    names(values) <- paste0("family:", names(specifications)[position], ":",
      if(is.null(specification$zeroDesign)) specification$parameter else names(values))
    result <- c(result, values)
  }
  result
}

.boot_mmes_refit <- function(object, response, maxit){
  if(.laplace_is_nominal(object)){
    prepared <- object$nominalPreparation
    prepared$nominal$response <- response
    prepared$covStruct <- object$covStruct
    prepared$initialParameters <- object$optimizer$par
    objective <- .laplace_nominal_objective(prepared)
    rejected <- 0L
    value <- function(parameters){
      evaluation <- objective$evaluate(parameters, gradient=FALSE)
      if(!isTRUE(evaluation$valid)){
        rejected <<- rejected + 1L
        return(Inf)
      }
      -evaluation$logLik
    }
    gradient <- function(parameters){
      evaluation <- objective$evaluate(parameters, gradient=TRUE)
      if(!isTRUE(evaluation$valid)){
        stop("Nominal bootstrap gradient evaluation failed: ",
          evaluation$reason, call.=FALSE)
      }
      -as.numeric(evaluation$gradient)
    }
    fit <- stats::nlminb(object$optimizer$par, value, gradient,
      control=list(iter.max=maxit,
        eval.max=max(1000, 10 * as.double(maxit))))
    evaluation <- objective$evaluate(fit$par, gradient=TRUE)
    if(!isTRUE(evaluation$valid)){
      stop("Nominal bootstrap final evaluation failed: ",
        evaluation$reason, call.=FALSE)
    }
    structures <- object$covStruct
    layout <- object$laplaceParameterLayout
    for(term in seq_along(layout$freePositions)){
      positions <- layout$freePositions[[term]]
      if(length(positions)){
        structures[[term]]$par[positions] <- fit$par[
          layout$fixedCount + layout$parameterStarts[term] +
            seq_along(positions)]
      }
    }
    boundary <- .mmes_laplace_boundary(structures, object$covStruct)
    return(list(parameters=.boot_mmes_parameters(object, fit$par),
      logLik=as.numeric(evaluation$logLik),
      converged=fit$convergence == 0L, message=fit$message,
      boundary=paste(boundary, collapse="; "),
      scoreMax=max(abs(evaluation$gradient)), rejected=rejected))
  }
  nativeResponse <- response
  for(position in seq_along(object$familySpecifications)){
    if(object$familySpecifications[[position]]$family != "binomial") next
    rows <- object$laplaceFamilyIndex == position
    nativeResponse[rows] <- round(response[rows] * object$priorWeights[rows])
  }
  settings <- object$laplacePostPEVSettings
  context <- sommer_laplace_objective_create(object$fixedDesign, object$randomDesign,
    nativeResponse, object$offset, "", object$covStruct, object$randomRelationships,
    settings$denseMemoryMB, familySpecs=object$familySpecifications,
    familyIndex=object$laplaceFamilyIndex, trials=object$priorWeights)
  evaluate <- function(parameters, gradient=FALSE){
    sommer_laplace_objective_eval(context, parameters,
      maxModeIterations=settings$modeMaxit, tolerance=settings$modeTol,
      gradient=gradient)
  }
  rejected <- 0L
  objective <- function(parameters){
    value <- evaluate(parameters)
    if(!isTRUE(value$valid)){
      rejected <<- rejected + 1L
      return(Inf)
    }
    -value$logLik
  }
  gradient <- function(parameters){
    value <- evaluate(parameters, TRUE)
    if(!isTRUE(value$valid)){
      stop("Bootstrap gradient evaluation failed: ", value$reason, call.=FALSE)
    }
    -as.numeric(value$gradient)
  }
  fit <- stats::nlminb(object$optimizer$par, objective, gradient,
    control=list(iter.max=maxit, eval.max=max(1000, 10 * as.double(maxit))))
  evaluation <- evaluate(fit$par, TRUE)
  if(!isTRUE(evaluation$valid)){
    stop("Bootstrap final evaluation failed: ", evaluation$reason, call.=FALSE)
  }
  layout <- object$laplaceParameterLayout
  structures <- object$covStruct
  for(term in seq_along(layout$freePositions)){
    positions <- layout$freePositions[[term]]
    structures[[term]]$par[positions] <- fit$par[layout$fixedCount +
      layout$parameterStarts[term] + seq_along(positions)]
  }
  boundary <- .mmes_laplace_boundary(structures, object$covStruct)
  specifications <- .laplace_proposal_specs(object$familySpecifications,
    fit$par, layout)
  for(position in seq_along(specifications)){
    specification <- specifications[[position]]
    if(specification$family == "ordinal" && specification$free &&
       any(diff(specification$thresholds) < sqrt(.Machine$double.eps))){
      boundary <- c(boundary, "response:threshold gap")
    }
    if(specification$family == "ziPoisson"){
      probability <- .laplace_zero_probability(specification,
        which(object$laplaceFamilyIndex == position))
      if(any(probability < sqrt(.Machine$double.eps) |
         probability > 1-sqrt(.Machine$double.eps))){
        boundary <- c(boundary, paste0("zero probability: ", names(specifications)[position]))
      }
    }
  }
  list(parameters=.boot_mmes_parameters(object, fit$par),
    logLik=as.numeric(evaluation$logLik), converged=fit$convergence == 0L,
    message=fit$message, boundary=paste(boundary, collapse="; "),
    scoreMax=max(abs(evaluation$gradient)), rejected=rejected)
}

boot_mmes <- function(object, object2=NULL, nsim=199L, seed=NULL, maxit=1000L){
  models <- if(is.null(object2)) list(object) else list(object, object2)
  if(any(!vapply(models, inherits, logical(1), what="mmes.glmm.laplace"))){
    stop("boot_mmes() requires native Laplace fits.", call.=FALSE)
  }
  if(any(!vapply(models, function(model) isTRUE(model$convergence), logical(1)))){
    stop("boot_mmes() requires converged generating and comparison fits.", call.=FALSE)
  }
  nsim <- .mme_positive(nsim, "nsim", TRUE)
  maxit <- .mme_positive(maxit, "maxit", TRUE)
  if(!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
     is.na(seed) || !is.finite(seed) || seed < 0 ||
     seed > .Machine$integer.max || seed != floor(seed))){
    stop("seed must be NULL or one nonnegative integer.", call.=FALSE)
  }
  comparison <- if(is.null(object2)) NULL else
    .anova_laplace_compare(object, object2, allowBoundary=TRUE)
  if(isTRUE(object$laplaceDiagnostics$boundary)){
    warning("The bootstrap generating fit has boundary parameters; results are conditional on this fitted boundary model and require calibration.",
      call.=FALSE)
  }
  hadSeed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if(!is.null(seed)){
    if(hadSeed) oldSeed <- get(".Random.seed", .GlobalEnv)
    on.exit(if(hadSeed) assign(".Random.seed", oldSeed, .GlobalEnv) else
      if(exists(".Random.seed", .GlobalEnv, inherits=FALSE))
        rm(".Random.seed", envir=.GlobalEnv), add=TRUE)
    set.seed(seed)
  }
  estimates <- .boot_mmes_parameters(object)
  bytes <- 8 * as.double(nsim) * (length(estimates) + 12)
  if(bytes > object$laplacePostPEVSettings$denseMemoryMB * 1048576){
    stop("Bootstrap output exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  draws <- matrix(NA_real_, nsim, length(estimates),
    dimnames=list(NULL, names(estimates)))
  diagnostics <- data.frame(replicate=seq_len(nsim), valid=FALSE,
    nullConverged=FALSE, alternativeConverged=NA, nullLogLik=NA_real_,
    alternativeLogLik=NA_real_, statistic=NA_real_, scoreMax=NA_real_,
    rejectedEvaluations=NA_integer_, boundary="", warnings="", error="")
  randomFactors <- .laplace_random_effect_sampling_factors(object,
    object$laplacePostPEVSettings$denseMemoryMB)
  for(index in seq_len(nsim)){
    warnings <- character()
    tryCatch(withCallingHandlers({
      response <- as.numeric(.laplace_simulate_response(object, nsim=1L,
        includeRandomEffects=TRUE,
        maxMemoryMB=object$laplacePostPEVSettings$denseMemoryMB,
        randomFactors=randomFactors))
      if(.laplace_is_nominal(object)){
        response <- match(object$nominal$levels[response],
          object$nominal$nativeLevels)
        if(anyNA(response)){
          stop("Nominal bootstrap simulation returned an undeclared category.",
            call.=FALSE)
        }
      }
      null <- .boot_mmes_refit(object, response, maxit)
      diagnostics$nullConverged[index] <- null$converged
      diagnostics$nullLogLik[index] <- null$logLik
      diagnostics$scoreMax[index] <- null$scoreMax
      diagnostics$rejectedEvaluations[index] <- null$rejected
      diagnostics$boundary[index] <- null$boundary
      draws[index, ] <- null$parameters
      if(!null$converged) stop("Generating-model refit did not converge: ", null$message)
      if(!is.null(object2)){
        alternative <- .boot_mmes_refit(object2, response, maxit)
        diagnostics$alternativeConverged[index] <- alternative$converged
        diagnostics$alternativeLogLik[index] <- alternative$logLik
        diagnostics$boundary[index] <- paste(null$boundary, alternative$boundary, sep=" | ")
        diagnostics$scoreMax[index] <- max(null$scoreMax, alternative$scoreMax)
        diagnostics$rejectedEvaluations[index] <- null$rejected + alternative$rejected
        if(!alternative$converged){
          stop("Alternative-model refit did not converge: ", alternative$message)
        }
        statistic <- 2*(alternative$logLik-null$logLik)
        if(statistic < -1e-7){
          stop("Alternative-model refit has lower likelihood than the null refit.")
        }
        diagnostics$statistic[index] <- max(0, statistic)
      }
      diagnostics$valid[index] <- TRUE
    }, warning=function(condition){
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }), error=function(condition){
      diagnostics$error[index] <<- conditionMessage(condition)
    })
    diagnostics$warnings[index] <- paste(unique(warnings), collapse=" | ")
  }
  valid <- diagnostics$valid
  failures <- sum(!valid)
  if(failures){
    warning(failures, " of ", nsim,
      " bootstrap replicates failed; inspect $diagnostics. Intervals use successful refits only; no unqualified p-value is returned.",
      call.=FALSE)
  }
  if(any(nzchar(diagnostics$warnings))){
    warning("Bootstrap refits produced warnings; inspect $diagnostics$warnings.", call.=FALSE)
  }
  intervals <- if(sum(valid) >= 2L)
    t(apply(draws[valid, , drop=FALSE], 2L, stats::quantile, probs=c(0.025, 0.975))) else
      matrix(NA_real_, length(estimates), 2L, dimnames=list(names(estimates), c("2.5%", "97.5%")))
  if(sum(valid) < 2L && !failures){
    warning("Fewer than two valid bootstrap refits; percentile intervals are unavailable.", call.=FALSE)
  }
  observed <- if(is.null(comparison)) NULL else comparison$Chisq[2L]
  exceedances <- if(is.null(observed)) NULL else
    sum(diagnostics$statistic[valid] >= observed)
  bounds <- if(is.null(observed)) NULL else
    c(lower=(1+exceedances)/(nsim+1), upper=(1+exceedances+failures)/(nsim+1))
  list(call=match.call(), nsim=nsim, successful=sum(valid), seed=seed,
    estimates=estimates, draws=draws, percentile=intervals,
    diagnostics=diagnostics, statistic=observed,
    p.value=if(is.null(observed) || failures) NA_real_ else bounds[[1L]],
    p.value.bounds=bounds,
    monteCarloSE=if(is.null(observed) || failures) NA_real_ else
      sqrt(bounds[[1L]]*(1-bounds[[1L]])/nsim),
    parameterScale="count coefficients and native covariance/family optimizer coordinates",
    method="Post-fit population parametric bootstrap; native Laplace refits")
}
