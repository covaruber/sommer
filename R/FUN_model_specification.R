# The specification engine deliberately uses numerical derivatives and dense,
# bounded workspaces. Covariance declarations still use the common vsm engine.
parameter <- function(fixed=~1, random=NULL, link="identity"){
  if(!inherits(fixed, "formula") || length(fixed)!=2L ||
     (!is.null(random) && (!inherits(random,"formula") || length(random)!=2L)))
    stop("parameter fixed and random must be one-sided formulas.",call.=FALSE)
  if(!is.character(link) || length(link)!=1L || is.na(link) ||
     !link %in% c("identity","log"))
    stop("parameter link must be 'identity' or 'log'.",call.=FALSE)
  structure(list(fixed=fixed,random=random,link=link),class="sommer_parameter")
}

graph_term <- function(name){
  if(!is.character(name) || length(name)!=1L || is.na(name) || !nzchar(name))
    stop("graph_term requires one non-empty predictor name.",call.=FALSE)
  stop("graph_term() is a formula marker, not a standalone evaluator.",call.=FALSE)
}

.spec_dots <- function(...){
  if(length(list(...))) stop("Unknown specification options are not supported.",call.=FALSE)
}

.spec_named <- function(x, label){
  if(length(x) && (is.null(names(x)) || anyNA(names(x)) ||
     any(!nzchar(names(x))) || anyDuplicated(names(x))))
    stop(label," must have unique non-empty names.",call.=FALSE)
}

.spec_split <- function(x, sign=1){
  if(is.call(x) && as.character(x[[1L]]) %in% c("+","-") && length(x)==3L)
    return(c(.spec_split(x[[2L]],sign),
      .spec_split(x[[3L]],if(as.character(x[[1L]])=="-") -sign else sign)))
  list(list(expression=x,sign=sign))
}

.spec_graph_names <- function(x){
  if(!is.call(x)) return(character())
  if(identical(x[[1L]],as.name("graph_term"))){
    if(length(x)!=2L || !is.character(x[[2L]]) || length(x[[2L]])!=1L ||
       is.na(x[[2L]]) || !nzchar(x[[2L]]))
      stop("Use graph_term('name') with a literal predictor name.",call.=FALSE)
    return(x[[2L]])
  }
  unique(unlist(lapply(as.list(x)[-1L],.spec_graph_names),use.names=FALSE))
}

.spec_validate_expression <- function(x, symbols, depth=0L){
  if(depth>32L) stop("Predictor expression exceeds nesting limit 32.",call.=FALSE)
  if(is.numeric(x) && length(x)==1L && is.finite(x)) return(invisible(NULL))
  if(is.symbol(x) && as.character(x) %in% symbols) return(invisible(NULL))
  if(!is.call(x) || !is.symbol(x[[1L]]))
    stop("Unsupported predictor expression.",call.=FALSE)
  op <- as.character(x[[1L]])
  if(op=="graph_term"){
    .spec_graph_names(x)
    return(invisible(NULL))
  }
  arity <- length(x)-1L
  allowed <- c("+","-","*","/","^","exp","log","tanh","sqrt","sin","cos","(")
  if(!op %in% allowed ||
     (op %in% c("+","-") && !arity %in% 1:2) ||
     (op %in% c("*","/","^") && arity!=2L) ||
     (!op %in% c("+","-","*","/","^") && arity!=1L))
    stop("Unsupported expression operation: ",op,".",call.=FALSE)
  for(child in as.list(x)[-1L]) .spec_validate_expression(child,symbols,depth+1L)
  invisible(NULL)
}

.spec_matrix <- function(formula,data,template=NULL){
  frame <- stats::model.frame(if(is.null(template)) formula else template$terms,
    data=data,na.action=stats::na.fail,
    xlev=if(!is.null(template)) template$xlevels else NULL)
  terms <- stats::terms(frame)
  x <- stats::model.matrix(terms,frame,
    contrasts.arg=if(!is.null(template)) template$contrasts else NULL)
  if(any(!is.finite(x))) stop("Non-finite fixed design.",call.=FALSE)
  if(!is.null(template) && !identical(colnames(x),colnames(template$x)))
    stop("Prediction fixed design differs from the training design.",call.=FALSE)
  list(x=x,terms=terms,xlevels=stats::.getXlevels(terms,frame),
    contrasts=attr(x,"contrasts"),
    offset=if(is.null(stats::model.offset(frame))) numeric(nrow(data)) else
      stats::model.offset(frame))
}

.spec_random <- function(formula,data,control){
  response <- ".sommer_spec_response"
  while(response %in% names(data)) response <- paste0(response,"_")
  training <- data
  training[[response]] <- seq_len(nrow(training))/nrow(training)
  fixed <- stats::as.formula(call("~",as.name(response),0),
    env=environment(formula))
  prepared <- .nonlinear_prepare(fixed,formula,training,control,
    "fail","fail","fail","fail",NULL)
  count <- length(prepared$rtermss)
  if(!count) stop("Random formula has no covariance terms.",call.=FALSE)
  mapping <- .laplace_covariance_map(prepared$covStruct[seq_len(count)],
    prepared$rtermss)
  dimensions <- vapply(mapping$structures,function(s) as.integer(s$dim),integer(1))
  ends <- cumsum(dimensions)
  blocks <- lapply(seq_len(count),function(i){
    index <- seq_len(dimensions[i])+(if(i==1L) 0L else ends[i-1L])
    block <- do.call(cbind,prepared$Z[index])
    colnames(block) <- make.unique(colnames(block))
    block
  })
  if(8*as.double(nrow(data))*sum(vapply(blocks,ncol,integer(1)))>
     control$execution$denseMemoryMB*1048576)
    stop("Random design exceeds execution$denseMemoryMB.",call.=FALSE)
  list(formula=formula,blocks=blocks,z=as.matrix(do.call(cbind,blocks)),
    relationships=prepared$Ai[seq_len(count)],mapping=mapping,
    variables=all.vars(formula)[all.vars(formula) %in% names(data)])
}

.spec_random_design <- function(record,newdata,training,control){
  variables <- record$variables
  if(any(!variables %in% names(newdata)))
    stop("Prediction is missing random-effect columns.",call.=FALSE)
  first <- newdata[,variables,drop=FALSE]
  last <- training[,variables,drop=FALSE]
  for(name in variables){
    if(is.factor(last[[name]]) || is.character(last[[name]])){
      levels <- if(is.factor(last[[name]])) levels(last[[name]]) else
        unique(as.character(last[[name]]))
      values <- as.character(first[[name]])
      if(anyNA(values) || any(!values %in% levels))
        stop("Unseen or missing random-effect level in '",name,"'.",call.=FALSE)
      first[[name]] <- factor(values,levels=levels)
      last[[name]] <- factor(as.character(last[[name]]),levels=levels)
    }
  }
  combined <- .spec_random(record$formula,rbind(first,last),control)
  if(length(combined$blocks)!=length(record$blocks))
    stop("Prediction random covariance basis changed.",call.=FALSE)
  blocks <- lapply(seq_along(record$blocks),function(i){
    expected <- colnames(record$blocks[[i]])
    actual <- combined$blocks[[i]]
    if(is.null(expected) || !setequal(expected,colnames(actual)))
      stop("Prediction random covariance basis changed.",call.=FALSE)
    actual[seq_len(nrow(newdata)),match(expected,colnames(actual)),drop=FALSE]
  })
  as.matrix(do.call(cbind,blocks))
}

.spec_covariance <- function(record,coordinates){
  structures <- record$mapping$structures
  map <- record$mapping$map
  working <- as.numeric(map$offset+map$T %*% coordinates)
  start <- 0L
  covariances <- lapply(seq_along(structures),function(i){
    structure <- structures[[i]]
    positions <- seq_along(structure$par)+start
    start <<- start+length(structure$par)
    structure$par <- working[positions]
    sigma <- evaluate_covstruct_cpp(structure,structure$par)
    relationship <- as.matrix(record$relationships[[i]])
    invisible(chol(relationship))
    covariance <- kronecker(sigma,chol2inv(chol(relationship)))
    if(nrow(covariance)!=ncol(record$blocks[[i]]))
      stop("Shared covariance descriptor and random design disagree.",call.=FALSE)
    covariance
  })
  as.matrix(Matrix::bdiag(covariances))
}

.spec_predictor <- function(spec,coefficients,data,design,random,mode){
  n <- nrow(data)
  value <- as.numeric(design$fixed$x %*% coefficients[spec$fixedPositions])+
    design$fixed$offset
  bindings <- new.env(parent=baseenv())
  for(name in spec$inputs){
    x <- data[[name]]
    if(!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
       length(x)!=n || any(!is.finite(x)))
      stop("Expression input '",name,"' must be finite numeric data.",call.=FALSE)
    bindings[[name]] <- x
  }
  values <- list()
  for(name in names(spec$parameters)){
    p <- spec$parameters[[name]]
    eta <- as.numeric(design$parameters[[name]]$x %*% coefficients[p$positions])
    index <- p$randomIndex
    if(length(index)) eta <- eta+as.numeric(random[[index]] %*% mode[[index]])
    values[[name]] <- if(p$link=="log") exp(eta) else eta
    bindings[[name]] <- values[[name]]
  }
  bindings$graph_term <- function(name){
    graph <- spec$predictors[[name]]
    parameters <- do.call(cbind,values[paste0(name,"::",graph$parameters)])
    colnames(parameters) <- graph$parameters
    as.numeric(.nonlinear_graph_values(graph,data,parameters,spec$control))
  }
  for(term in spec$expressions){
    current <- suppressWarnings(eval(term$expression,bindings))
    if(!is.numeric(current) || !length(current) %in% c(1L,n))
      stop("Expression must return one numeric value per observation.",call.=FALSE)
    value <- value+term$sign*current
  }
  if(length(spec$additiveIndex))
    value <- value+as.numeric(random[[spec$additiveIndex]] %*%
      mode[[spec$additiveIndex]])
  value
}

.spec_logdensity <- function(y,eta,family,variance){
  mu <- family$linkinv(eta)
  if(any(!is.finite(mu)) || any(!family$validmu(mu))) return(-Inf)
  switch(family$family,
    gaussian=sum(stats::dnorm(y,mu,sqrt(variance),log=TRUE)),
    Gamma=sum(stats::dgamma(y,shape=1/variance,scale=mu*variance,log=TRUE)),
    poisson=sum(stats::dpois(y,mu,log=TRUE)),
    binomial=sum(stats::dbinom(y,1,mu,log=TRUE)))
}

.spec_gradient <- function(at,fn){
  vapply(seq_along(at),function(i){
    step <- .Machine$double.eps^(1/3)*(1+abs(at[i]))
    upper <- lower <- at
    upper[i] <- upper[i]+step
    lower[i] <- lower[i]-step
    (fn(upper)-fn(lower))/(2*step)
  },numeric(1))
}

.spec_mode <- function(initial,conditional,gradientFunction,settings){
  inner <- stats::optim(initial,conditional,gr=gradientFunction,method="BFGS",
    control=list(maxit=settings$maxIters,reltol=settings$scoreTol))
  at <- inner$par
  # A likelihood change criterion is not a mode certificate. Refine against
  # the actual numerical score, without regularizing conditional curvature.
  remaining <- max(0L,settings$maxIters-inner$counts[["gradient"]])
  for(iteration in seq_len(remaining+1L)){
    gradient <- gradientFunction(at)
    if(any(!is.finite(gradient)))
      stop("Conditional random-mode score is non-finite.",call.=FALSE)
    hessian <- stats::optimHess(at,conditional,gr=gradientFunction)
    curvature <- .mmes_laplace_curvature(hessian)
    if(!curvature$valid)
      stop("Conditional mode Hessian is not positive definite; Laplace integration is undefined.",call.=FALSE)
    if(max(abs(gradient))<=settings$scoreTol)
      return(list(par=at,value=conditional(at),hessian=hessian,
        score=max(abs(gradient)),iterations=iteration))
    if(iteration>remaining) break
    step <- as.numeric(curvature$covariance %*% gradient)
    value <- conditional(at)
    fraction <- 1
    repeat{
      candidate <- at-fraction*step
      candidateValue <- conditional(candidate)
      if(is.finite(candidateValue) &&
         (candidateValue<=value-1e-4*fraction*sum(step*gradient) ||
          (candidateValue<=value+16*.Machine$double.eps*(1+abs(value)) &&
           max(abs(gradientFunction(candidate)))<max(abs(gradient))))) break
      fraction <- fraction/2
      if(max(abs(fraction*step)/(1+abs(at)))<settings$stepTol)
        stop("Conditional mode stalled before meeting the requested score tolerance.",call.=FALSE)
    }
    at <- candidate
  }
  stop("Conditional mode failed the requested numerical score tolerance.",call.=FALSE)
}

.spec_evaluate <- function(at,context,details=FALSE){
  p <- context$p
  q <- context$q
  n <- context$n
  widths <- context$widths
  variancePosition <- context$variancePosition
  covariancePositions <- context$covariancePositions
  records <- context$records
  spec <- context$spec
  data <- context$data
  design <- context$design
  randomDesign <- context$randomDesign
  zeros <- context$zeros
  y <- context$y
  family <- context$family
  exact <- context$exact
  control <- context$control
  splitMode <- function(mode){
    ends <- cumsum(widths)
    lapply(seq_along(widths),function(i)
      mode[seq_len(widths[i])+(if(i==1L) 0L else ends[i-1L])])
  }
  coefficients <- at[seq_len(p)]
  variance <- if(length(variancePosition)) unname(exp(at[variancePosition])) else 1
  if(!is.finite(variance) || variance<=0) return(if(details) NULL else Inf)
  covariances <- lapply(seq_along(records),function(i)
    .spec_covariance(records[[i]],at[covariancePositions[[i]]]))
  eta <- .spec_predictor(spec,coefficients,data,design,randomDesign,zeros)
  if(any(!is.finite(eta))) return(if(details) NULL else Inf)
  modeScore <- 0
  if(!q){
    value <- -.spec_logdensity(y,eta,family,variance)
    mode <- numeric()
  }else if(exact){
    z <- do.call(cbind,randomDesign)
    g <- as.matrix(Matrix::bdiag(covariances))
    v <- z %*% g %*% t(z)+diag(variance,n)
    factor <- chol(v)
    solved <- backsolve(factor,forwardsolve(t(factor),y-eta))
    value <- .5*(n*log(2*pi)+2*sum(log(diag(factor)))+sum((y-eta)*solved))
    mode <- as.numeric(g %*% crossprod(z,solved))
  }else{
    g <- as.matrix(Matrix::bdiag(covariances))
    factor <- chol(g)
    precision <- chol2inv(factor)
    conditional <- function(b){
      predictor <- .spec_predictor(spec,coefficients,data,design,randomDesign,splitMode(b))
      if(any(!is.finite(predictor))) return(Inf)
      -.spec_logdensity(y,predictor,family,variance)+
        .5*sum(b*(precision %*% b))
    }
    conditionalGradient <- function(b){
      predictor <- .spec_predictor(spec,coefficients,data,design,randomDesign,splitMode(b))
      mu <- family$linkinv(predictor)
      responseScore <- switch(family$family,
        gaussian=(mu-y)/variance,
        poisson=1-y/mu,
        binomial=(mu-y)/(mu*(1-mu)),
        Gamma=(mu-y)/(variance*mu^2))*family$mu.eta(predictor)
      derivatives <- vapply(seq_len(q),function(i){
        step <- .Machine$double.eps^(1/3)*(1+abs(b[i]))
        upper <- lower <- b
        upper[i] <- upper[i]+step
        lower[i] <- lower[i]-step
        (.spec_predictor(spec,coefficients,data,design,randomDesign,splitMode(upper))-
          .spec_predictor(spec,coefficients,data,design,randomDesign,splitMode(lower)))/(2*step)
      },numeric(n))
      as.numeric(crossprod(derivatives,responseScore)+precision %*% b)
    }
    inner <- .spec_mode(numeric(q),conditional,conditionalGradient,
      control$task$nonlinear$mode)
    mode <- inner$par
    modeScore <- inner$score
    value <- inner$value+sum(log(diag(factor)))+
      .5*as.numeric(determinant(inner$hessian,logarithm=TRUE)$modulus)
  }
  if(details) list(value=value,mode=splitMode(mode),covariances=covariances,
    eta=.spec_predictor(spec,coefficients,data,design,randomDesign,splitMode(mode)),
    variance=variance,modeScore=modeScore) else value
}

.mmes_specification_fit <- function(fixed,random=NULL,rcov=NULL,data,
    family=stats::gaussian(),parameters=list(),predictors=list(),start=NULL,
    control,specificationMethod=NULL,weights=NULL,W=NULL,vcc=NULL,covPar=NULL,
    familyParameters=list(),
    contrasts=NULL,naMethodX="exclude",naMethodY="exclude",
    naMethodRandom="exclude",naMethodR="exclude",...){
  .spec_dots(...)
  if(!is.list(familyParameters) || length(familyParameters))
    stop("Specification fits do not support nonempty familyParameters.",call.=FALSE)
  if(!inherits(control,"sommer_mmes_control"))
    stop("control must be an internal mmes control object.",call.=FALSE)
  if(!is.null(weights) || !is.null(W) || !is.null(vcc) ||
     !is.null(covPar) || !is.null(contrasts))
    stop("Specification fits do not support weights, W, vcc, covPar, or custom contrasts.",call.=FALSE)
  defaults <- .mmes_control_build()
  for(module in c("representation","solve","logdet","score","curvature",
      "update","proposal","output")){
    if(!identical(control[[module]],defaults[[module]]))
      stop("Specification fitting does not support non-default ",module," controls.",call.=FALSE)
  }
  if(control$formulation!="henderson")
    stop("Specification fitting does not support a requested direct formulation.",call.=FALSE)
  modeSettings <- control$task$nonlinear$mode
  if(modeSettings$warmStart!="off" || modeSettings$starts!="deterministic")
    stop("Specification modes require no warm starts and deterministic starts.",call.=FALSE)
  public <- attr(control,"publicControl")
  if(!is.null(public)){
    settings <- attr(public,"modelControls")
    if(!is.null(settings)){
      if(settings$optimizer$gradient=="analytic")
        stop("Specification expressions support numerical gradients, not requested analytic gradients.",call.=FALSE)
      if(!identical(settings$working,.mme_pql_defaults) ||
         !identical(settings$likelihood,
           list(backend="auto",trace="panels",blocks="off")) ||
         !identical(settings$inference,.mme_nonlinear_inference_defaults))
        stop("Specification fitting does not support working, likelihood backend, or post-inference controls.",call.=FALSE)
    }
  }
  integration <- "auto"
  if(!is.null(specificationMethod)){
    if(!is.list(specificationMethod) ||
       !specificationMethod$integration %in% c("auto","exact","laplace") ||
       !specificationMethod$approximation %in% c("auto","likelihood","pql"))
      stop("Invalid specificationMethod.",call.=FALSE)
    if(specificationMethod$approximation=="pql")
      stop("Specification likelihood models do not support PQL.",call.=FALSE)
    integration <- specificationMethod$integration
  }
  if(!inherits(fixed,"formula") || length(fixed)!=3L ||
     !is.symbol(fixed[[2L]]))
    stop("fixed must be a two-sided formula with a response column.",call.=FALSE)
  if(!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("data must be a non-empty data frame with unique column names.",call.=FALSE)
  if(!is.list(parameters) || !is.list(predictors))
    stop("parameters and predictors must be named lists.",call.=FALSE)
  .spec_named(parameters,"parameters")
  .spec_named(predictors,"predictors")
  if(!is.null(rcov) && (!inherits(rcov,"formula") || length(rcov)!=2L ||
     !identical(rcov[[2L]],as.name("units"))))
    stop("Specification fits support independent residuals: rcov=NULL or ~units.",call.=FALSE)
  if(!inherits(family,"family") ||
     !family$family %in% c("gaussian","binomial","poisson","Gamma") ||
     (family$family=="Gamma" && family$link!="log"))
    stop("Specification fitting supports Gaussian, Bernoulli binomial, Poisson, and Gamma-log families.",call.=FALSE)
  policies <- list(X=naMethodX,Y=naMethodY,Random=naMethodRandom,R=naMethodR)
  if(any(!vapply(policies,function(x) is.character(x) && length(x)==1L &&
      !is.na(x) && x %in% c("exclude","fail"),logical(1))))
    stop("Specification naMethod policies support only 'exclude' and 'fail'.",call.=FALSE)
  parameterVariables <- unlist(lapply(parameters,function(p){
    if(!inherits(p,"sommer_parameter")) return(character())
    all.vars(p$fixed)
  }),use.names=FALSE)
  randomVariables <- unique(c(all.vars(random),unlist(lapply(parameters,function(p){
    if(!inherits(p,"sommer_parameter")) return(character())
    all.vars(p$random)
  }),use.names=FALSE)))
  graphVariables <- unlist(lapply(predictors,function(graph){
    if(!inherits(graph,"sommer_nonlinear_graph")) return(character())
    unlist(lapply(graph$nodes,function(node)
      if(identical(node$op,"input")) node$column else character()),use.names=FALSE)
  }),use.names=FALSE)
  groups <- list(Y=as.character(fixed[[2L]]),
    X=unique(c(all.vars(fixed[[3L]]),parameterVariables,graphVariables)),
    Random=randomVariables,R=character())
  keep <- rep(TRUE,nrow(data))
  for(group in names(groups)){
    columns <- intersect(groups[[group]],names(data))
    missing <- if(length(columns)) !stats::complete.cases(data[,columns,drop=FALSE]) else
      rep(FALSE,nrow(data))
    if(any(missing) && policies[[group]]=="fail")
      stop("Missing values violate naMethod",group,"='fail'.",call.=FALSE)
    keep <- keep & !missing
  }
  included <- which(keep)
  excluded <- which(!keep)
  if(!length(included)) stop("No complete observations remain.",call.=FALSE)
  data <- data[included,,drop=FALSE]
  if(control$task$mode!="fit")
    stop("Specification fitting requires task$mode='fit'.",call.=FALSE)
  terms <- .spec_split(fixed[[3L]])
  graphs <- unique(unlist(lapply(terms,function(term)
    .spec_graph_names(term$expression)),use.names=FALSE))
  if(!setequal(graphs,names(predictors)))
    stop("predictors must name exactly the graph_term markers in fixed.",call.=FALSE)
  for(name in graphs){
    graph <- .nonlinear_graph_validate(predictors[[name]])
    if(family$family!="gaussian" || family$link!="identity" ||
       !identical(graph$nodes[[graph$output]]$op,"identity"))
      stop("Graph specification fits require Gaussian identity family and an identity output node.",call.=FALSE)
    if(length(graph$random) || !is.null(graph$joint) ||
       any(vapply(graph$parameterMap,function(x) identical(x$type,"random"),logical(1))))
      stop("Graph internal random maps are unsupported; declare namespaced parameter(random=) models.",call.=FALSE)
    for(coordinate in graph$parameters){
      key <- paste0(name,"::",coordinate)
      if(is.null(parameters[[key]]))
        parameters[[key]] <- parameter(link=graph$parameterMap[[coordinate]]$transform)
    }
    naturalMap <- lapply(graph$parameterMap,function(map){
      map$transform <- "identity"
      map
    })
    predictors[[name]] <- nonlinear_graph(graph$nodes,graph$output,
      graph$parameters,naturalMap)
  }
  if(any(!vapply(parameters,inherits,logical(1),"sommer_parameter")))
    stop("Every parameter declaration must be created by parameter().",call.=FALSE)
  if(any(names(parameters) %in% names(data)))
    stop("Parameter names must not shadow data columns.",call.=FALSE)
  nonlinear <- vapply(terms,function(term)
    length(intersect(all.vars(term$expression),names(parameters)))>0L ||
      length(.spec_graph_names(term$expression))>0L,logical(1))
  if(!any(nonlinear)){
    if(length(parameters)) stop("Declared parameters are unused in fixed.",call.=FALSE)
    if(!is.null(random))
      stop("Ordinary mixed models belong to the baseline fitting engine.",call.=FALSE)
    if(!is.null(start))
      stop("Starts are not used by fixed-only ordinary models.",call.=FALSE)
    base <- if(identical(family$family,"gaussian") && identical(family$link,"identity"))
      stats::lm(fixed,data=data,na.action=stats::na.fail,model=TRUE,x=TRUE,y=TRUE) else
      stats::glm(fixed,data=data,family=family,na.action=stats::na.fail,
        model=TRUE,x=TRUE,y=TRUE)
    if(anyNA(stats::coef(base))) stop("Fixed design is rank deficient.",call.=FALSE)
    n <- stats::nobs(base)
    gaussian <- identical(family$family,"gaussian")
    reml <- gaussian && !inherits(base,"glm") && control$task$objective=="reml"
    variance <- if(gaussian) sum(stats::residuals(base,type="response")^2)/
      (n-if(reml) length(stats::coef(base)) else 0L) else NULL
    likelihood <- if(inherits(base,"glm")) stats::logLik(base) else
      stats::logLik(base,REML=reml)
    if(family$family=="Gamma"){
      if(any(base$y<=0)) stop("Gamma responses must be positive.",call.=FALSE)
      phi <- stats::nlminb(log(mean(((base$y-stats::fitted(base))/
          stats::fitted(base))^2)),function(logPhi)
        -.spec_logdensity(base$y,base$linear.predictors,family,exp(logPhi)),
        control=list(iter.max=control$convergence$maxIters,
          rel.tol=control$convergence$likelihoodTol,
          x.tol=control$convergence$parameterTol))
      if(phi$convergence!=0L)
        stop("Gamma dispersion ML optimization did not converge.",call.=FALSE)
      variance <- exp(phi$par)
      likelihood <- structure(-phi$objective,df=length(stats::coef(base))+1L)
    }
    result <- list(baseFit=base,b=stats::coef(base),fixedFormula=fixed,
      fitted.values=stats::fitted(base),linear.predictors=if(inherits(base,"glm"))
        base$linear.predictors else stats::fitted(base),
      y=base$y,trainingData=data,family=family,control=control,
      llik=as.numeric(likelihood),nobs=n,npar=attr(likelihood,"df"),
      dispersion=if(gaussian || family$family=="Gamma") variance else 1,
      convergence=if(inherits(base,"glm")) base$converged else TRUE,
      likelihoodMethod=if(reml) "exact Gaussian REML" else "exact fixed-only ML",
      derivativeMethod="base lm/glm",obsInfo=list(included=included,excluded=excluded))
    class(result) <- c("mmes.specification","mmes")
    return(result)
  }
  if(control$task$objective!="ml")
    stop("Nonlinear specification models require task$objective='ml'.",call.=FALSE)
  if(64*as.double(nrow(data))^2>control$execution$denseMemoryMB*1048576)
    stop("Observation workspace exceeds execution$denseMemoryMB.",call.=FALSE)
  if(length(parameters)>30L || length(terms)>100L)
    stop("Specification exceeds 30 parameters or 100 additive terms.",call.=FALSE)
  used <- unique(c(unlist(lapply(terms[nonlinear],function(term)
    intersect(all.vars(term$expression),names(parameters))),use.names=FALSE),
    unlist(lapply(graphs,function(name)
      paste0(name,"::",predictors[[name]]$parameters)),use.names=FALSE)))
  if(!setequal(used,names(parameters))) stop("Declared parameters are unused in fixed.",call.=FALSE)
  expressions <- terms[nonlinear]
  for(term in expressions)
    .spec_validate_expression(term$expression,c(names(data),names(parameters)))
  inputs <- setdiff(unique(unlist(lapply(expressions,function(term)
    all.vars(term$expression)),use.names=FALSE)),names(parameters))
  rhs <- if(length(terms[!nonlinear])){
    ordinary <- terms[!nonlinear]
    out <- if(ordinary[[1L]]$sign==1) ordinary[[1L]]$expression else
      call("-",ordinary[[1L]]$expression)
    if(length(ordinary)>1L) for(term in ordinary[-1L])
      out <- call(if(term$sign==1) "+" else "-",out,term$expression)
    out
  }else 1
  ordinary <- stats::as.formula(call("~",rhs),env=environment(fixed))
  fixedDesign <- .spec_matrix(ordinary,data)
  design <- list(fixed=fixedDesign,parameters=list())
  positions <- seq_len(ncol(fixedDesign$x))
  coefficientNames <- paste0("fixed:",colnames(fixedDesign$x))
  if(!length(positions)) coefficientNames <- character()
  records <- list()
  for(name in names(parameters)){
    p <- parameters[[name]]
    design$parameters[[name]] <- .spec_matrix(p$fixed,data)
    size <- ncol(design$parameters[[name]]$x)
    if(!size) stop("Parameter regressions require at least one fixed coefficient.",call.=FALSE)
    p$positions <- length(coefficientNames)+seq_len(size)
    coefficientNames <- c(coefficientNames,paste0(name,":",
      colnames(design$parameters[[name]]$x)))
    p$randomIndex <- integer()
    if(!is.null(p$random)){
      records[[length(records)+1L]] <- .spec_random(p$random,data,control)
      p$randomIndex <- length(records)
    }
    parameters[[name]] <- p
  }
  additiveIndex <- integer()
  if(!is.null(random)){
    records[[length(records)+1L]] <- .spec_random(random,data,control)
    additiveIndex <- length(records)
  }
  y <- data[[as.character(fixed[[2L]])]]
  if(!is.numeric(y) || is.complex(y) || !is.null(dim(y)) ||
     length(y)!=nrow(data) || any(!is.finite(y)))
    stop("Response must be a finite numeric vector without missing values.",call.=FALSE)
  if(family$family=="binomial" && any(!y %in% 0:1))
    stop("Binomial specification fits currently require Bernoulli 0/1 responses.",call.=FALSE)
  if(family$family=="poisson" && any(y<0 | y!=floor(y)))
    stop("Poisson responses must be non-negative integers.",call.=FALSE)
  if(family$family=="Gamma" && any(y<=0))
    stop("Gamma responses must be positive.",call.=FALSE)
  n <- length(y)
  p <- length(coefficientNames)
  widths <- vapply(records,function(record) ncol(record$z),integer(1))
  q <- sum(widths)
  exact <- family$family=="gaussian" && family$link=="identity" &&
    all(vapply(parameters,function(p) !length(p$randomIndex),logical(1)))
  if(q && !exact && integration=="exact")
    stop("Exact integration is unavailable for parameter random effects or non-Gaussian random models.",call.=FALSE)
  if(p>60L || (!exact && q>512L) || n<=p ||
     8*(8*as.double(n)^2+8*as.double(q)^2+as.double(n)*(p+q))>
       control$execution$denseMemoryMB*1048576)
    stop("Specification exceeds observation, coefficient, mode, or dense memory limits.",call.=FALSE)
  spec <- list(parameters=parameters,predictors=predictors,expressions=expressions,
    inputs=inputs,fixedPositions=positions,additiveIndex=additiveIndex,control=control)
  zeros <- lapply(widths,numeric)
  randomDesign <- lapply(records,`[[`,"z")
  coefficients <- stats::setNames(numeric(p),coefficientNames)
  if(is.list(start)){
    if(any(!names(start) %in% c("parameters","coefficients")) ||
       length(start)!=1L || is.null(names(start)))
      stop("start must contain exactly parameters or coefficients.",call.=FALSE)
    explicit <- identical(names(start),"coefficients")
    start <- start[[1L]]
  }else explicit <- FALSE
  if(!is.numeric(start) || is.complex(start) || !is.null(dim(start)) ||
     any(!is.finite(start)))
    stop("Nonlinear fitting requires a complete finite named start vector.",call.=FALSE)
  .spec_named(start,"start")
  if(explicit){
    if(!setequal(names(start),coefficientNames))
      stop("Coefficient starts must name exactly: ",paste(coefficientNames,collapse=", "),".",call.=FALSE)
    coefficients[] <- start[coefficientNames]
  }else{
    if(!setequal(names(start),names(parameters)))
      stop("Natural parameter starts must name exactly: ",paste(names(parameters),collapse=", "),".",call.=FALSE)
    for(name in names(parameters)){
      declaration <- parameters[[name]]
      value <- start[[name]]
      if(declaration$link=="log"){
        if(value<=0) stop("Log-linked parameter starts must be positive.",call.=FALSE)
        value <- log(value)
      }
      x <- design$parameters[[name]]$x
      intercept <- match("(Intercept)",colnames(x))
      if(is.na(intercept))
        stop("Natural parameter starts require intercepts; use explicit coefficient starts.",call.=FALSE)
      coefficients[declaration$positions[intercept]] <- value
    }
    initial <- .spec_predictor(spec,coefficients,data,design,randomDesign,zeros)
    if(length(positions)){
      if(family$family=="gaussian" && family$link=="identity")
        coefficients[positions] <- stats::lm.fit(fixedDesign$x,y-initial)$coefficients
      else coefficients[positions] <- stats::glm.fit(fixedDesign$x,y,
        offset=initial,family=family)$coefficients
    }
  }
  if(any(!is.finite(coefficients))) stop("Initial mean design is rank deficient.",call.=FALSE)
  initialEta <- .spec_predictor(spec,coefficients,data,design,randomDesign,zeros)
  if(any(!is.finite(initialEta))) stop("Starts produce a non-finite predictor.",call.=FALSE)
  jacobian <- vapply(seq_len(p),function(i){
    upper <- lower <- coefficients
    step <- 1e-5*(1+abs(coefficients[i]))
    upper[i] <- upper[i]+step
    lower[i] <- lower[i]-step
    (.spec_predictor(spec,upper,data,design,randomDesign,zeros)-
      .spec_predictor(spec,lower,data,design,randomDesign,zeros))/(2*step)
  },numeric(n))
  if(any(!is.finite(jacobian)) || qr(jacobian,tol=1e-8)$rank<p)
    stop("Mean specification is locally rank deficient at the starts.",call.=FALSE)
  covariancePositions <- list()
  initial <- coefficients
  for(i in seq_along(records)){
    coordinates <- records[[i]]$mapping$map$initial
    covariancePositions[[i]] <- length(initial)+seq_along(coordinates)
    if(length(coordinates))
      initial <- c(initial,stats::setNames(coordinates,paste0("covariance",i,":",names(coordinates))))
  }
  variancePosition <- integer()
  if(family$family %in% c("gaussian","Gamma")){
    initialMean <- family$linkinv(initialEta)
    variance <- if(family$family=="Gamma")
      mean(((y-initialMean)/initialMean)^2) else mean((y-initialMean)^2)
    if(!is.finite(variance) || variance<=0)
      stop("Starting residual variance must be positive.",call.=FALSE)
    variancePosition <- length(initial)+1L
    initial <- c(initial,`residual:logVariance`=log(variance))
  }
  context <- list(p=p,q=q,n=n,widths=widths,variancePosition=variancePosition,
    covariancePositions=covariancePositions,records=records,spec=spec,data=data,
    design=design,randomDesign=randomDesign,zeros=zeros,y=y,family=family,
    exact=exact,control=control)
  evaluate <- function(at,details=FALSE) .spec_evaluate(at,context,details)
  if(!is.finite(evaluate(initial))) stop("Starts give an invalid likelihood.",call.=FALSE)
  optimizer <- stats::nlminb(initial,evaluate,
    gradient=function(at) .spec_gradient(at,evaluate),
    control=list(iter.max=control$convergence$maxIters,
      eval.max=20L*control$convergence$maxIters,
      rel.tol=control$convergence$likelihoodTol,
      x.tol=control$convergence$parameterTol))
  optimum <- evaluate(optimizer$par,TRUE)
  if(is.null(optimum) || !is.finite(optimum$value))
    stop("Specification optimization did not produce a finite likelihood.",call.=FALSE)
  if(optimizer$convergence!=0L)
    warning("Specification optimizer did not converge; inspect optimizer.",call.=FALSE)
  result <- list(b=optimizer$par[seq_len(p)],fixedFormula=fixed,randomFormula=random,
    family=family,trainingData=data,control=control,specification=spec,
    design=design,randomRecords=records,randomMode=optimum$mode,
    randomCovariances=optimum$covariances,optimizer=optimizer,
    likelihoodContext=context,modeScore=optimum$modeScore,
    linear.predictors=optimum$eta,fitted.values=family$linkinv(optimum$eta),
    y=y,nobs=n,npar=length(optimizer$par),llik=-optimum$value,
    dispersion=optimum$variance,convergence=optimizer$convergence==0L,
    likelihoodMethod=if(!q || exact) "exact ML" else "Laplace ML",
    derivativeMethod="numerical finite differences (outer score, mode score and conditional Hessian)",
    obsInfo=list(included=included,excluded=excluded))
  class(result) <- c("mmes.specification","mmes")
  result
}

coef.mmes.specification <- function(object,...){
  .spec_dots(...)
  object$b
}

fitted.mmes.specification <- function(object,type=c("response","link"),...){
  .spec_dots(...)
  if(match.arg(type)=="link") object$linear.predictors else object$fitted.values
}

residuals.mmes.specification <- function(object,type=c("response","pearson"),...){
  .spec_dots(...)
  residual <- object$y-object$fitted.values
  if(match.arg(type)=="pearson")
    residual <- residual/sqrt(object$family$variance(object$fitted.values)*object$dispersion)
  residual
}

predict.mmes.specification <- function(object,newdata=NULL,
    type=c("response","link"),target=c("conditional","zero","population"),...){
  .spec_dots(...)
  type <- match.arg(type)
  target <- match.arg(target)
  if(!is.null(object$baseFit)){
    if(is.null(newdata)) return(fitted.mmes.specification(object,type))
    return(if(inherits(object$baseFit,"glm"))
      stats::predict(object$baseFit,newdata=newdata,type=type) else
      stats::predict(object$baseFit,newdata=newdata))
  }
  if(is.null(newdata)) newdata <- object$trainingData
  if(!is.data.frame(newdata) || !nrow(newdata))
    stop("newdata must be a non-empty data frame.",call.=FALSE)
  design <- list(fixed=.spec_matrix(NULL,newdata,object$design$fixed),
    parameters=lapply(object$design$parameters,function(template)
      .spec_matrix(NULL,newdata,template)))
  records <- object$randomRecords
  if(target=="zero") random <- lapply(records,function(record)
    matrix(0,nrow(newdata),ncol(record$z))) else
    random <- lapply(records,.spec_random_design,newdata=newdata,
      training=object$trainingData,control=.mmes_fitting_control(object))
  mode <- if(target=="conditional") object$randomMode else
    lapply(records,function(record) numeric(ncol(record$z)))
  eta <- .spec_predictor(object$specification,object$b,newdata,design,random,mode)
  if(target!="population" || !length(records))
    return(if(type=="link") eta else object$family$linkinv(eta))
  # Integrate the low-dimensional predictor coordinates, not all group effects.
  active <- which(vapply(random,function(z) any(z!=0),logical(1)))
  if(length(active)>2L)
    stop("Population quadrature currently supports at most two random predictor coordinates.",call.=FALSE)
  variance <- lapply(active,function(i)
    rowSums((random[[i]] %*% object$randomCovariances[[i]])*random[[i]]))
  settings <- .mmes_fitting_control(object)$task$nonlinear$population
  integrate <- function(order){
    rule <- .laplace_normal_quadrature(order)
    grid <- expand.grid(rep(list(seq_along(rule$nodes)),length(active)))
    answer <- numeric(nrow(newdata))
    for(row in seq_len(nrow(grid))){
      shifts <- vector("list",length(records))
      for(i in seq_along(records)) shifts[[i]] <- numeric(nrow(newdata))
      weight <- 1
      for(j in seq_along(active)){
        shifts[[active[j]]] <- sqrt(variance[[j]])*rule$nodes[grid[row,j]]
        weight <- weight*rule$weights[grid[row,j]]
      }
      shifted <- object$specification
      value <- as.numeric(design$fixed$x %*% object$b[shifted$fixedPositions])+
        design$fixed$offset
      bindings <- new.env(parent=baseenv())
      for(name in shifted$inputs) bindings[[name]] <- newdata[[name]]
      values <- list()
      for(name in names(shifted$parameters)){
        p <- shifted$parameters[[name]]
        parameterEta <- as.numeric(design$parameters[[name]]$x %*% object$b[p$positions])
        if(length(p$randomIndex)) parameterEta <- parameterEta+shifts[[p$randomIndex]]
        values[[name]] <- if(p$link=="log") exp(parameterEta) else parameterEta
        bindings[[name]] <- values[[name]]
      }
      bindings$graph_term <- function(name){
        graph <- shifted$predictors[[name]]
        values <- do.call(cbind,values[paste0(name,"::",graph$parameters)])
        colnames(values) <- graph$parameters
        as.numeric(.nonlinear_graph_values(graph,newdata,values,
          .mmes_fitting_control(object)))
      }
      for(term in shifted$expressions)
        value <- value+term$sign*eval(term$expression,bindings)
      if(length(shifted$additiveIndex)) value <- value+shifts[[shifted$additiveIndex]]
      answer <- answer+weight*(if(type=="link") value else object$family$linkinv(value))
    }
    answer
  }
  order <- settings$order
  previous <- integrate(order)
  repeat{
    nextOrder <- min(2L*order,settings$maxOrder)
    if(nextOrder==order) stop("Population quadrature requires maxOrder greater than order.",call.=FALSE)
    current <- integrate(nextOrder)
    if(all(is.finite(current)) &&
       all(abs(current-previous)<=settings$tol*(1+abs(current)))) return(current)
    if(nextOrder>=settings$maxOrder)
      stop("Population quadrature did not converge at maxOrder.",call.=FALSE)
    order <- nextOrder
    previous <- current
  }
}

logLik.mmes.specification <- function(object,...){
  .spec_dots(...)
  structure(object$llik,class="logLik",nobs=object$nobs,df=object$npar)
}

nobs.mmes.specification <- function(object,...){
  .spec_dots(...)
  object$nobs
}

vcov.mmes.specification <- function(object,...){
  .spec_dots(...)
  if(!is.null(object$baseFit)){
    covariance <- stats::vcov(object$baseFit)
    if(object$family$family=="gaussian"){
      covariance <- covariance*object$dispersion/
        (sum(stats::residuals(object$baseFit,type="response")^2)/
          stats::df.residual(object$baseFit))
    }
    if(object$family$family=="Gamma")
      covariance <- covariance*object$dispersion/summary(object$baseFit)$dispersion
    return(covariance)
  }
  hessian <- stats::optimHess(object$optimizer$par,function(at)
    .spec_evaluate(at,object$likelihoodContext))
  curvature <- .mmes_laplace_curvature(hessian)
  if(!curvature$valid)
    stop("Coefficient covariance unavailable: ",curvature$reason,".",call.=FALSE)
  size <- length(object$b)
  result <- curvature$covariance[seq_len(size),seq_len(size),drop=FALSE]
  dimnames(result) <- list(names(object$b),names(object$b))
  result
}

simulate.mmes.specification <- function(object,nsim=1,seed=NULL,...){
  .spec_dots(...)
  stop("Simulation is not implemented for specification fits.",call.=FALSE)
}

summary.mmes.specification <- function(object,...){
  .spec_dots(...)
  structure(list(coefficients=object$b,logLik=object$llik,nobs=object$nobs,
    convergence=object$convergence,likelihoodMethod=object$likelihoodMethod,
    derivativeMethod=object$derivativeMethod,dispersion=object$dispersion),
    class="summary.mmes.specification")
}

print.summary.mmes.specification <- function(x,...){
  .spec_dots(...)
  cat("Model specification:",x$likelihoodMethod,"\nConverged:",x$convergence,
    "\nlogLik:",format(x$logLik),"\n")
  print(x$coefficients)
  cat("Derivatives:",x$derivativeMethod,"\n")
  invisible(x)
}

print.mmes.specification <- function(x,...){
  print(summary.mmes.specification(x),...)
  invisible(x)
}
