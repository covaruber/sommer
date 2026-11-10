#' Traditional nonlinear mixed-model curves
#'
#' Declare a native parametric curve to add to the ordinary fixed and random
#' predictors in \code{mmes()}. No layers, callbacks or learned basis are used.
#' @param input Name of one numeric column in \code{data}.
#' @param curve One of \code{"exponential"}, \code{"logistic"} or
#'   \code{"saturation"}.
#' @param random Parameter-named one-sided random formulas for independent
#'   Gaussian curve deviations. Existing \code{vsm()} syntax is unchanged.
#' @param joint Optional list with \code{parameters} and one \code{random}
#'   formula using a covariance constructor on the internal \code{.parameter}
#'   axis, for example \code{~vsm(usm(.parameter),ism(group))}.
#'   Gaussian deviations act on the log scale for positive parameters and on
#'   the identity scale for midpoint.
#' @details The curves are \eqn{A\exp(-kx)},
#'   \eqn{A/(1+\exp(-k(x-m)))} and \eqn{Ax/(K+x)}, respectively.
#'   Parameters \code{asymptote}, \code{rate} and \code{half} are positive;
#'   \code{midpoint} is unconstrained. Positive parameters are optimized on
#'   the log scale. Exponential and saturation inputs must be non-negative.
#'   Starting values, on the natural scale, belong to
#'   \code{mmes(start=list(parameters=...))}.
#' @return A serializable \code{sommer_nonlinear_curve} specification.
#' @export
nonlinear_curve <- function(input, curve=c("exponential", "logistic", "saturation"),
    random=list(), joint=NULL){
  if(!is.character(input) || length(input) != 1L || is.na(input) || !nzchar(input)){
    stop("input must name one numeric data column.", call.=FALSE)
  }
  curve <- match.arg(curve)
  parameters <- switch(curve,
    exponential=c("asymptote", "rate"),
    logistic=c("asymptote", "rate", "midpoint"),
    saturation=c("asymptote", "half"))
  if(!is.list(random) ||
     (length(random) && (is.null(names(random)) || anyNA(names(random)) ||
      anyDuplicated(names(random)) || any(!names(random) %in% parameters) ||
      any(!vapply(random,function(formula) inherits(formula,"formula") &&
        length(formula)==2L && length(.mmes_split_plus(formula[[2L]]))==1L,
        logical(1)))))){
    stop("random must contain uniquely curve-parameter-named, one-term random formulas.",
      call.=FALSE)
  }
  if(!is.null(joint)){
    if(!is.list(joint) || length(joint)!=2L || anyDuplicated(names(joint)) ||
       !setequal(names(joint),c("parameters","random")) ||
       !is.character(joint$parameters) || length(joint$parameters)<2L ||
       anyNA(joint$parameters) || anyDuplicated(joint$parameters) ||
       any(!joint$parameters %in% parameters) ||
       any(joint$parameters %in% names(random)) ||
       !inherits(joint$random,"formula") || length(joint$random)!=2L ||
       length(.mmes_split_plus(joint$random[[2L]]))!=1L ||
       !".parameter" %in% all.vars(joint$random))
      stop("joint requires unique curve parameters disjoint from random and one random formula containing .parameter.",
        call.=FALSE)
    expression <- joint$random[[2L]]
    if(!is.call(expression) || !identical(expression[[1L]],as.name("vsm")))
      stop("joint$random must use vsm() with .parameter as a covariance-shaping factor.",
        call.=FALSE)
    arguments <- as.list(expression)[-1L]
    if(!is.null(names(arguments))) arguments <- arguments[!names(arguments) %in%
      c("Gu","sigma2","fixedSigma2","rotation","isFixed","verbose")]
    axis <- which(vapply(arguments,function(value)
      ".parameter" %in% all.vars(value),logical(1)))
    if(length(axis)!=1L || axis>=length(arguments))
      stop("Joint .parameter must identify exactly one covariance-shaping axis.",
        call.=FALSE)
    joint <- joint[c("parameters","random")]
    joint$parameters <- unname(joint$parameters)
  }
  out <- list(version=if(length(random)>1L || !is.null(joint)) 3L else
    if(length(random)) 2L else 1L, input=input, curve=curve,
    parameters=parameters, positive=parameters != "midpoint")
  if(length(random)) out$random <- random
  if(!is.null(joint)) out$joint <- joint
  structure(out, class="sommer_nonlinear_curve")
}

.nonlinear_validate <- function(model){
  if(!inherits(model, "sommer_nonlinear_curve") || !is.list(model) ||
     length(model$version)!=1L || is.na(model$version) ||
     !model$version %in% 1:3 ||
     !is.character(model$curve) || length(model$curve) != 1L ||
     is.na(model$curve) ||
     !model$curve %in% c("exponential", "logistic", "saturation")){
    stop("nonlinear must be a specification created by nonlinear_curve().", call.=FALSE)
  }
  expected <- nonlinear_curve(model$input, model$curve,
    if(is.null(model$random)) list() else model$random, model$joint)
  if(!identical(unclass(model), unclass(expected))){
    stop("The nonlinear curve parameter layout has been modified.", call.=FALSE)
  }
  model
}

.nonlinear_input <- function(model, data, missingOK=FALSE){
  if(!is.data.frame(data) || !model$input %in% names(data)){
    stop("Nonlinear input '", model$input, "' must be a column in data.", call.=FALSE)
  }
  x <- data[[model$input]]
  if(!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
     any(is.infinite(x)) || (!missingOK && anyNA(x))){
    stop("Nonlinear input '", model$input, "' must contain finite numeric values.",
      call.=FALSE)
  }
  if(model$curve != "logistic" && any(x < 0, na.rm=TRUE)){
    stop(model$curve, " curve inputs must be non-negative.", call.=FALSE)
  }
  x
}

.nonlinear_coordinates <- function(model, start){
  if(!is.numeric(start) || is.complex(start) || !is.null(dim(start)) ||
     length(start) != length(model$parameters) ||
     is.null(names(start)) || anyNA(names(start)) || anyDuplicated(names(start)) ||
     !setequal(names(start), model$parameters) || any(!is.finite(start))){
    stop("task$nonlinear$start must be a finite named vector containing exactly: ",
      paste(model$parameters, collapse=", "), ".", call.=FALSE)
  }
  start <- start[model$parameters]
  if(any(start[model$positive] <= 0)){
    stop("Positive nonlinear parameters require strictly positive starting values.",
      call.=FALSE)
  }
  start[model$positive] <- log(start[model$positive])
  start
}

.nonlinear_natural <- function(model, coordinates){
  names(coordinates) <- model$parameters
  coordinates[model$positive] <- exp(coordinates[model$positive])
  coordinates
}

.nonlinear_start <- function(model, x, y, settings){
  if(!is.null(settings$start)) return(.nonlinear_coordinates(model, settings$start))
  if(settings$initialization == "supplied"){
    stop("initialization='supplied' requires task$nonlinear$start.", call.=FALSE)
  }
  span <- diff(range(x))
  if(span <= 0){
    stop("Nonlinear fitting requires variation in the input.", call.=FALSE)
  }
  amplitude <- max(y)
  if(!is.finite(amplitude) || amplitude <= 0){
    stop("Automatic curve initialization requires a positive response maximum; supply task$nonlinear$start.",
      call.=FALSE)
  }
  start <- switch(model$curve,
    exponential=c(asymptote=amplitude, rate=1/span),
    logistic=c(asymptote=1.1*amplitude, rate=4/span, midpoint=stats::median(x)),
    saturation=c(asymptote=1.2*amplitude, half=span/2))
  .nonlinear_coordinates(model, start)
}

.nonlinear_prediction_formula <- function(formula, data){
  parent <- environment(formula)
  captured <- new.env(parent=parent)
  for(variable in setdiff(all.vars(formula), c(names(data), "units"))){
    if(exists(variable, envir=parent, inherits=TRUE)){
      assign(variable, get(variable, envir=parent, inherits=TRUE), envir=captured)
    }
  }
  environment(formula) <- captured
  formula
}
