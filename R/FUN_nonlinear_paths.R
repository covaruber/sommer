.mme_nonlinear_mode_defaults <- list(maxIters=100L, scoreTol=1e-8, stepTol=1e-8,
  warmStart="off", starts="deterministic")
.mme_nonlinear_prediction_defaults <- list(order=20L, maxOrder=160L, tol=1e-7)
.mme_nonlinear_inference_defaults <- list(step=1e-4, derivativeTol=5e-3,
  curvatureTol=1e-8, maxCondition=1e10)
.mme_nonlinear_defaults <- list(optimizer="nlminb", gradient="analytic",
  initialization="auto", start=NULL, mode=.mme_nonlinear_mode_defaults,
  population=.mme_nonlinear_prediction_defaults,
  inference=.mme_nonlinear_inference_defaults)

.nonlinear_controls <- function(value){
  value <- .mme_module(value, .mme_nonlinear_defaults, "task$nonlinear")
  value$optimizer <- .mme_choice(value$optimizer, "nlminb", "task$nonlinear$optimizer")
  value$gradient <- .mme_choice(value$gradient, "analytic", "task$nonlinear$gradient")
  value$initialization <- .mme_choice(value$initialization, c("auto", "supplied"),
    "task$nonlinear$initialization")
  value$mode <- .mme_module(value$mode, .mme_nonlinear_mode_defaults,
    "task$nonlinear$mode")
  for(name in c("maxIters", "scoreTol", "stepTol")){
    value$mode[[name]] <- .mme_positive(value$mode[[name]],
      paste0("task$nonlinear$mode$", name), integer=name == "maxIters")
  }
  value$mode$warmStart <- .mme_choice(value$mode$warmStart, "off",
    "task$nonlinear$mode$warmStart")
  value$mode$starts <- .mme_choice(value$mode$starts, "deterministic",
    "task$nonlinear$mode$starts")
  value$population <- .mme_module(value$population, .mme_nonlinear_prediction_defaults,
    "task$nonlinear$population")
  for(name in c("order", "maxOrder", "tol")){
    value$population[[name]] <- .mme_positive(value$population[[name]],
      paste0("task$nonlinear$population$", name), integer=name != "tol")
  }
  if(value$population$order < 2L || value$population$maxOrder < 2*value$population$order ||
     value$population$maxOrder > 512L){
    stop("Population quadrature requires order>=2 and 2*order<=maxOrder<=512.",
      call.=FALSE)
  }
  value$inference <- .mme_module(value$inference, .mme_nonlinear_inference_defaults,
    "task$nonlinear$inference")
  for(name in names(value$inference))
    value$inference[[name]] <- .mme_positive(value$inference[[name]],
      paste0("task$nonlinear$inference$",name))
  if(value$inference$step>=0.1 || value$inference$derivativeTol>=1 ||
     value$inference$curvatureTol>=1 || value$inference$maxCondition<=1)
    stop("Nonlinear inference requires step<0.1, derivativeTol<1, curvatureTol<1 and maxCondition>1.",
      call.=FALSE)
  if(!is.null(value$start) &&
     (!is.numeric(value$start) || is.complex(value$start) ||
      !is.null(dim(value$start)) || !length(value$start) ||
      any(!is.finite(value$start)) || is.null(names(value$start)) ||
      anyNA(names(value$start)) || any(!nzchar(names(value$start))) ||
      anyDuplicated(names(value$start)))){
    stop("task$nonlinear$start must be a finite numeric vector with unique non-empty names.",
      call.=FALSE)
  }
  value
}

.nonlinear_exact_gaussian <- function(family, control){
  inherits(family, "family") && identical(family$family, "gaussian") &&
    identical(family$link, "identity") && !length(control$task$glmm$familyParameters)
}

.nonlinear_path_plan <- function(model, control, family, random, rcov, W,
    weights, covPar, internal){
  .nonlinear_validate(model)
  .laplace_family_specs(family, control$task$glmm$familyParameters)
  gaussian <- .nonlinear_exact_gaussian(family, control)
  if(control$task$objective != "ml" || control$task$mode != "fit" ||
     control$formulation != "henderson"){
    stop("Nonlinear fitting requires task$objective='ml', task$mode='fit' and formulation='henderson'.",
      call.=FALSE)
  }
  if(is.null(random) && !length(model$random) && is.null(model$joint)){
    stop("Nonlinear mixed models require at least one additive Gaussian random term.",
      call.=FALSE)
  }
  if(gaussian && !length(model$random) && is.null(model$joint) &&
     !identical(control$task$nonlinear$mode, .mme_nonlinear_mode_defaults)){
    stop("Nonlinear mode controls require random curve parameters.", call.=FALSE)
  }
  if(gaussian && !length(model$random) && is.null(model$joint) &&
     !identical(control$task$nonlinear$population, .mme_nonlinear_prediction_defaults)){
    stop("Nonlinear population controls require random curve parameters.", call.=FALSE)
  }
  if(!is.null(rcov) &&
     (!inherits(rcov, "formula") || length(rcov) != 2L ||
      !identical(rcov[[2L]], as.name("units")))){
    stop("Nonlinear fitting currently supports only independent homogeneous residuals (rcov=~units).",
      call.=FALSE)
  }
  if(!is.null(W) || !is.null(weights) || !is.null(covPar) || !is.null(internal)){
    stop("Nonlinear fitting does not support W, weights, covPar or internal working models; use random covariance descriptors for covariance starts.",
      call.=FALSE)
  }
  if(control$output$inverse != 0L){
    stop("Nonlinear fits defer PEV and uncertainty until post-fit; output$inverse must be 0.",
      call.=FALSE)
  }
  defaults <- .mmes_control_internal()
  checked <- control
  checked$task$glmm$familyParameters <- defaults$task$glmm$familyParameters
  allowed <- c("task$mode", "task$objective", "task$nonlinear",
    "task$glmm$familyParameters",
    "convergence$maxIters", "convergence$likelihoodTol", "convergence$parameterTol",
    "execution$denseMemoryMB", "execution$blasThreads", "diagnostics$verbose",
    "diagnostics$profile")
  paths <- setdiff(.mme_control_paths(control), allowed)
  changed <- paths[vapply(paths, function(path)
    !identical(.mme_path_value(checked, path), .mme_path_value(defaults, path)),
    logical(1))]
  # Age-message suppression is shared preparation policy, not a fitting method.
  changed <- setdiff(changed, "diagnostics$dateWarning")
  forced <- .mme_mode_paths[vapply(.mme_mode_paths, function(path)
    identical(.mme_path_value(control, path), "force"), logical(1))]
  changed <- unique(c(changed, setdiff(forced,allowed),
    if(control$execution$gpu) "execution$gpu"))
  if(length(changed)){
    stop("Nonlinear fitting does not use these controls: ",
      paste(changed, collapse=", "), ".", call.=FALSE)
  }
  list(requested="traditional nonlinear ML", eligible=TRUE,
    selected=if(!gaussian) "Laplace family ML" else
      if(length(model$random) || !is.null(model$joint)) "Laplace Gaussian ML" else "exact Gaussian ML",
    curve=model$curve,
    random=if(length(model$random) || !is.null(model$joint)) "Gaussian curve-coordinate deviations" else
        "additive linear Gaussian", backend="random",
    optimizer=control$task$nonlinear$optimizer,
    gradient=control$task$nonlinear$gradient,
    initialization=if(is.null(control$task$nonlinear$start)) "auto" else "supplied",
    uncertainty="not computed during fitting; use post_nonlinear() or postPEV()")
}

.nonlinear_no_inference <- function(object, method){
  if(inherits(object,"mmes.specification"))
    stop(method," is not supported for formula specification fits; coefficient covariance is available through vcov().",
      call.=FALSE)
  if(inherits(object, "mmes.nonlinear")){
    stop(method, " is not yet supported for nonlinear models; PEV and uncertainty are post-fit operations.",
      call.=FALSE)
  }
  invisible(NULL)
}

.nonlinear_covariance_eligibility <- function(structures, terms){
  for(index in seq_along(structures)){
    for(factor in structures[[index]]$factors){
      if((isTRUE(factor$zero_specific) ||
          isTRUE(factor$representation$loadings$zero_specific)) &&
         identical(factor$derivative$backend,"native") &&
         identical(factor$derivative$op,"rr")){
        next
      }
      if(!factor$evaluator$backend %in% c("native", "fixed") ||
         !factor$derivative$backend %in% c("native", "none")){
        stop("Nonlinear fitting requires native analytic covariance descriptors; unsupported evaluator or derivative in ",
          terms[index], ".", call.=FALSE)
      }
    }
  }
  invisible(NULL)
}
