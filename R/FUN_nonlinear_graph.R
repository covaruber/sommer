#' Declare a bounded nonlinear predictor graph
#'
#' Create a serializable directed acyclic graph for evaluation and likelihood
#' mixed-model fitting with fixed or group-random graph coordinates.
#' @param nodes Ordered list of named \code{input}, \code{affine},
#'   \code{identity}, or \code{tanh} nodes. Each node may refer only to earlier
#'   nodes. Input nodes use \code{column}, optional \code{center}, and positive
#'   \code{scale}. Affine nodes use \code{inputs}, an integer \code{weights}
#'   matrix, and an integer \code{bias} vector. Weight rows follow concatenated
#'   input-node order; columns follow output units. Parameter indices are
#'   one-based, and repeated indices share coordinates. Zero bias indices omit
#'   an intercept.
#' @param output Name of one scalar output node.
#' @param parameters Unique parameter-coordinate names. A \code{log}
#'   transform treats supplied values as log coordinates.
#' @param parameterMap Optional named list mapping coordinates to
#'   \code{list(type="fixed")} or \code{list(type="random", group="column")}.
#'   Each coordinate may also set \code{transform="log"}.
#' @param random Optional named list of one-term random formulas, one for each
#'   independently distributed random graph coordinate. The formula declares
#'   its covariance using existing \code{vsm()} syntax.
#' @param joint Optional list with \code{parameters} and one random formula
#'   containing the internal \code{.parameter} axis to declare a joint
#'   covariance block, following \code{nonlinear_curve()}.
#' @details Use fixed graph coordinates with
#'   \code{mmes(y ~ graph_term("net"), predictors=list(net=graph))} for ML. Supply complete named
#'   natural-parameter starts through
#'   \code{mmes(start=list(parameters=...))}. Starts use natural parameter
#'   values; the standalone evaluator continues to use working coordinates. For random
#'   coordinates, mark them in \code{parameterMap} and declare their
#'   independent or joint covariance with \code{random} or \code{joint}.
#'   Gaussian deviations are on the declared working-coordinate scale
#'   (including log coordinates), and are integrated with an exact conditional
#'   Hessian and a Laplace likelihood. Outer Laplace scores for random graph
#'   fits use central finite differences of that objective and are reported in
#'   fit diagnostics. Conditional and zero-deviation predictions, prior-
#'   integrated population predictions, and simulation are supported.
#'   Supported response families are Gaussian identity, binomial
#'   logit/probit/cloglog (including grouped trials), Poisson, NB2 and ZIP log,
#'   Gamma log, ordinal logit/probit and reference-category multinomial,
#'   including \code{familym()} mixtures. The graph supplies a scalar
#'   link predictor per expanded category row; use explicit category-specific
#'   regressions or covariance maps when required. Free ordinal cutpoint
#'   locations must not alias nonlinear mean parameters.
#'   Use \code{tanh} only in hidden layers; fits require an identity output.
#'   Hidden-unit permutations and scaling can make raw weights non-unique;
#'   convergence does not establish a global optimum. Explicit
#'   \code{post_nonlinear()}, \code{vcov()}, confidence intervals and prediction
#'   uncertainty require locally identified, stationary, interior fits with
#'   checked positive curvature. Singular weight symmetries are rejected;
#'   conditional-mode covariance and approximate empirical Bayes PEV are
#'   different quantities. Namespaced \code{parameter()} overrides select
#'   the bounded dense reference route, which rejects internal graph random
#'   maps and reduced-rank covariance instead of adding a ridge.
#' @return A versioned, serializable \code{sommer_nonlinear_graph}.
#' @examples
#' graph <- nonlinear_graph(list(
#'   x=list(op="input", column="x"),
#'   mean=list(op="affine", inputs="x", weights=matrix(1L), bias=2L),
#'   output=list(op="identity", input="mean")
#' ), output="output", parameters=c("slope", "intercept"))
#' nonlinear_graph_eval(graph, data.frame(x=c(0, 1, 2)),
#'   c(slope=0.5, intercept=1))
#' fitData <- data.frame(x=seq_len(30), group=factor(rep(1:5, each=6)))
#' set.seed(2024)
#' fitData$y <- 0.4*fitData$x + rep(c(-0.2, 0.1, 0.3, -0.1, -0.1), each=6) +
#'   rnorm(nrow(fitData), sd=0.1)
#' fitGraph <- nonlinear_graph(list(
#'   x=list(op="input", column="x"),
#'   mean=list(op="affine", inputs="x", weights=matrix(1L), bias=2L),
#'   output=list(op="identity", input="mean")
#' ), output="output", parameters=c("slope", "intercept"))
#' fit <- mmes(y~0+graph_term("net"), random=~group, data=fitData,
#'   predictors=list(net=fitGraph), method=mmes_method(objective="ml"),
#'   start=list(parameters=c(slope=0.4, intercept=0)))
#' @export
nonlinear_graph <- function(nodes, output, parameters, parameterMap=NULL,
    random=list(), joint=NULL){
  if(!is.list(nodes) || !length(nodes) ||
     !is.character(output) || length(output)!=1L || is.na(output) ||
     !nzchar(output) || !is.character(parameters) || !length(parameters) ||
     anyNA(parameters) || any(!nzchar(parameters)) || anyDuplicated(parameters)){
    stop("nodes, output and unique parameter names are required.",call.=FALSE)
  }
  if(is.null(parameterMap)){
    parameterMap <- stats::setNames(lapply(parameters,function(name)
      list(type="fixed",transform="identity")),parameters)
  }
  if(!is.list(parameterMap) || is.null(names(parameterMap)) ||
     anyNA(names(parameterMap)) || anyDuplicated(names(parameterMap)) ||
     !setequal(names(parameterMap),parameters)){
    stop("parameterMap must define every parameter exactly once.",call.=FALSE)
  }
  parameterMap <- parameterMap[parameters]
  parameterMap <- lapply(seq_along(parameterMap),function(i){
    entry <- parameterMap[[i]]
    name <- parameters[[i]]
    if(!is.list(entry) || !is.character(entry$type) || length(entry$type)!=1L ||
       is.na(entry$type) || !entry$type %in% c("fixed","random") ||
       !is.null(entry$transform) &&
         (!is.character(entry$transform) || length(entry$transform)!=1L ||
           is.na(entry$transform) ||
           !entry$transform %in% c("identity","log"))){
      stop("Invalid parameter map for '",name,"'.",call.=FALSE)
    }
    transform <- if(is.null(entry$transform)) "identity" else entry$transform
    if(identical(entry$type,"fixed")){
      if(!is.null(entry$group)) stop("Fixed parameters cannot name a group.",call.=FALSE)
      list(type="fixed",transform=transform)
    }else{
      if(!is.character(entry$group) || length(entry$group)!=1L ||
         is.na(entry$group) || !nzchar(entry$group))
        stop("Random parameters require one group-column name.",call.=FALSE)
      list(type="random",group=entry$group,transform=transform)
    }
  })
  names(parameterMap) <- parameters
  randomCoordinates <- parameters[vapply(parameterMap,function(map)
    identical(map$type,"random"),logical(1))]
  if(!is.list(random) ||
     (length(random) && (is.null(names(random)) || anyNA(names(random)) ||
       anyDuplicated(names(random)) ||
       any(!names(random) %in% randomCoordinates) ||
       any(!vapply(random,function(formula) inherits(formula,"formula") &&
         length(formula)==2L && length(.mmes_split_plus(formula[[2L]]))==1L,
         logical(1)))))){
    stop("random must contain uniquely named one-term formulas for independent random graph coordinates.",
      call.=FALSE)
  }
  if(!is.null(joint)){
    if(!is.list(joint) || length(joint)!=2L || anyDuplicated(names(joint)) ||
       !setequal(names(joint),c("parameters","random")) ||
       !is.character(joint$parameters) || length(joint$parameters)<2L ||
       anyNA(joint$parameters) || anyDuplicated(joint$parameters) ||
       any(!joint$parameters %in% randomCoordinates) ||
       any(joint$parameters %in% names(random)) ||
       !inherits(joint$random,"formula") || length(joint$random)!=2L ||
       length(.mmes_split_plus(joint$random[[2L]]))!=1L ||
       !".parameter" %in% all.vars(joint$random)){
      stop("joint requires unique random graph parameters disjoint from random and one random formula containing .parameter.",
        call.=FALSE)
    }
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
  hasCovarianceMap <- length(random)>0L || !is.null(joint)
  if(hasCovarianceMap &&
     !setequal(c(names(random),if(!is.null(joint)) joint$parameters),
       randomCoordinates))
    stop("Every random graph coordinate must be declared exactly once in random or joint.",
      call.=FALSE)
  for(parameter in names(random)){
    group <- parameterMap[[parameter]]$group
    if(!group %in% all.vars(random[[parameter]]))
      stop("Random formula for graph parameter '",parameter,
        "' must include its declared group column '",group,"'.",call.=FALSE)
  }
  if(!is.null(joint)){
    groups <- vapply(joint$parameters,function(parameter)
      parameterMap[[parameter]]$group,character(1))
    if(length(unique(groups))!=1L || !groups[1L] %in% all.vars(joint$random))
      stop("Joint random graph coordinates must declare the same group column and the joint formula must include it.",
        call.=FALSE)
  }
  widths <- integer()
  compiled <- vector("list",length(nodes))
  names(compiled) <- names(nodes)
  if(is.null(names(nodes)) || anyNA(names(nodes)) || any(!nzchar(names(nodes))) ||
     anyDuplicated(names(nodes)))
    stop("nodes must have unique non-empty names.",call.=FALSE)
  for(i in seq_along(nodes)){
    node <- nodes[[i]]
    name <- names(nodes)[[i]]
    if(!is.list(node) || !is.character(node$op) || length(node$op)!=1L ||
       is.na(node$op) || !nzchar(node$op))
      stop("Each graph node requires one operation.",call.=FALSE)
    op <- node$op
    if(op=="input"){
      if(!is.character(node$column) || length(node$column)!=1L ||
         is.na(node$column) || !nzchar(node$column))
        stop("Input node '",name,"' requires one column.",call.=FALSE)
      center <- if(is.null(node$center)) 0 else node$center
      scale <- if(is.null(node$scale)) 1 else node$scale
      if(!is.numeric(center) || length(center)!=1L || !is.finite(center) ||
         !is.numeric(scale) || length(scale)!=1L || !is.finite(scale) || scale<=0)
        stop("Input transforms require finite center and positive scale.",call.=FALSE)
      compiled[[i]] <- list(op=op,column=node$column,center=as.numeric(center),
        scale=as.numeric(scale),inputs=character(),width=1L)
    }else if(op=="affine"){
      inputs <- node$inputs
      if(!is.character(inputs) || !length(inputs) || anyNA(inputs) ||
         any(!inputs %in% names(widths)))
        stop("Affine node '",name,"' must reference earlier graph nodes.",call.=FALSE)
      inputWidth <- sum(widths[inputs])
      weights <- node$weights
      bias <- node$bias
      if(!is.matrix(weights) || !is.numeric(weights) && !is.integer(weights) ||
         nrow(weights)!=inputWidth || !ncol(weights) ||
         anyNA(weights) || any(weights<1) || any(weights>length(parameters)) ||
         !is.numeric(bias) && !is.integer(bias) ||
         length(bias)!=ncol(weights) || anyNA(bias) || any(bias<0) ||
         any(bias>length(parameters)) ||
         any(weights!=as.integer(weights)) || any(bias!=as.integer(bias)))
        stop("Affine node '",name,"' has invalid parameter indices or dimensions.",
          call.=FALSE)
      compiled[[i]] <- list(op=op,inputs=inputs,
        weights=matrix(as.integer(weights),nrow=nrow(weights)),
        bias=as.integer(bias),width=as.integer(ncol(weights)))
    }else if(op %in% c("identity","tanh")){
      input <- if(is.null(node$input)) node$inputs else node$input
      if(!is.character(input) || length(input)!=1L || is.na(input) ||
         !input %in% names(widths))
        stop("Activation node '",name,"' must reference an earlier graph node.",
          call.=FALSE)
      compiled[[i]] <- list(op=op,inputs=input,width=widths[[input]])
    }else{
      stop("Unsupported nonlinear graph operation: ",op,".",call.=FALSE)
    }
    widths[[name]] <- compiled[[i]]$width
  }
  if(!output %in% names(widths) || widths[[output]]!=1L)
    stop("Graph output must name one scalar node.",call.=FALSE)
  result <- list(version=if(hasCovarianceMap) 2L else 1L,
    nodes=compiled,output=output,parameters=unname(parameters),
    parameterMap=parameterMap)
  if(length(random)) result$random <- random
  if(!is.null(joint)) result$joint <- joint
  structure(result,class="sommer_nonlinear_graph")
}

#' Evaluate a nonlinear predictor graph and derivative products
#'
#' Evaluate a declared graph without fitting. Return native reverse-mode
#' gradients and directional derivatives through order three.
#' @param graph A graph created by \code{nonlinear_graph()}.
#' @param data Data frame containing graph input and random-group columns.
#' @param parameters Named parameter vector, or observation-by-parameter matrix.
#'   Random coordinates in a matrix must be constant within each group.
#' @param directions Optional parameter-by-three matrix shared across rows, or
#'   observation-by-parameter-by-three array for group-specific directions.
#'   Fixed-coordinate directions must be constant across rows; random-coordinate
#'   directions must be constant within their declared groups. Omission sets
#'   all directions to zero.
#' @param control Validated \code{mmes_control()} object for workspace limits.
#' @return A list with an output vector and observation-by-parameter gradient,
#'   Hessian-vector product and third-order directional gradient matrices.
#'   Directional arrays have dimensions observations-by-3,
#'   observations-by-3-by-3 and observations-by-3-by-3-by-3. Directional
#'   outputs are zero when \code{directions} is omitted.
#' @export
nonlinear_graph_eval <- function(graph,data,parameters,directions=NULL,
    control=mmes_control()){
  graph <- .nonlinear_graph_validate(graph)
  if(!is.data.frame(data) || !nrow(data))
    stop("data must be a non-empty data frame.",call.=FALSE)
  if(!inherits(control,"sommer_mmes_control"))
    stop("control must be created by mmes_control().",call.=FALSE)
  control <- do.call(.mmes_control_build,.mmes_engine_arguments(control))
  required <- unique(c(vapply(graph$nodes,function(node)
      if(identical(node$op,"input")) node$column else character(1),
    character(1)),vapply(graph$parameterMap,function(map)
      if(identical(map$type,"random")) map$group else character(1),
    character(1))))
  required <- required[nzchar(required)]
  if(any(!required %in% names(data)))
    stop("Graph data is missing columns: ",
      paste(setdiff(required,names(data)),collapse=", "),".",call.=FALSE)
  requiredBytes <- 8*as.double(nrow(data))*
    (sum(vapply(graph$nodes,function(node) node$width,numeric(1)))+
       7*length(graph$parameters)+
       sum(vapply(graph$nodes,function(node) identical(node$op,"input"),
         logical(1)))+16)
  if(requiredBytes>control$execution$denseMemoryMB*1048576)
    stop("Nonlinear graph workspace exceeds execution$denseMemoryMB.",call.=FALSE)
  inputNodes <- which(vapply(graph$nodes,function(node)
    identical(node$op,"input"),logical(1)))
  inputs <- do.call(cbind,lapply(graph$nodes[inputNodes],function(node){
    value <- data[[node$column]]
    if(!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
       length(value)!=nrow(data) || any(!is.finite(value)))
      stop("Graph input '",node$column,"' must contain finite numeric values.",
        call.=FALSE)
    as.numeric(value)
  }))
  if(is.null(dim(inputs))) inputs <- matrix(inputs,ncol=1L)

  colnames(inputs) <- names(graph$nodes)[vapply(graph$nodes,function(node)
    identical(node$op,"input"),logical(1))]
  if(!is.matrix(parameters) && is.numeric(parameters) && is.null(dim(parameters))){
    if(is.null(names(parameters)) || anyNA(names(parameters)) ||
       anyDuplicated(names(parameters)) ||
       !setequal(names(parameters),graph$parameters))
      stop("parameters must be named with every graph coordinate.",call.=FALSE)
    parameterValues <- matrix(unname(parameters[graph$parameters]),
      nrow=nrow(data),ncol=length(graph$parameters),byrow=TRUE)
  }else if(is.matrix(parameters) && is.numeric(parameters) &&
      nrow(parameters)==nrow(data) && ncol(parameters)==length(graph$parameters) &&
      !is.null(colnames(parameters)) && setequal(colnames(parameters),graph$parameters)){
    parameterValues <- unname(parameters[,graph$parameters,drop=FALSE])
  }else{
    stop("parameters must be a complete named vector or observation-by-parameter matrix.",
      call.=FALSE)
  }
  if(any(!is.finite(parameterValues)))
    stop("Graph parameters must be finite.",call.=FALSE)
  for(j in seq_along(graph$parameterMap)){
    map <- graph$parameterMap[[j]]
    if(identical(map$transform,"log") && any(parameterValues[,j]>log(.Machine$double.xmax)))
      stop("Log-transformed graph parameter exceeds numeric range.",call.=FALSE)
    if(identical(map$type,"fixed") && any(parameterValues[,j]!=parameterValues[1L,j]))
      stop("Fixed graph parameters must be constant across observations.",call.=FALSE)
    if(identical(map$type,"random")){
      group <- data[[map$group]]
      if(!is.atomic(group) || is.list(group) || !is.null(dim(group)) ||
         length(group)!=nrow(data) ||
         !(is.factor(group) || is.character(group) || is.numeric(group) ||
           is.logical(group)))
        stop("Random graph group labels must be one atomic vector.",call.=FALSE)
      if(is.numeric(group) && any(!is.finite(group)))
        stop("Random graph group labels must be complete and finite.",call.=FALSE)
      if(anyNA(group))
        stop("Random graph group labels must not be missing.",call.=FALSE)
      for(level in unique(group)){
        rows <- which(group==level)
        if(any(parameterValues[rows,j]!=parameterValues[rows[1L],j]))
          stop("Random graph parameters must be constant within group '",
            map$group,"'.",call.=FALSE)
      }
    }
  }
  logParameters <- which(vapply(graph$parameterMap,function(map)
    identical(map$transform,"log"),logical(1)))
  if(length(logParameters)){
    natural <- exp(parameterValues[,logParameters,drop=FALSE])
    if(any(!is.finite(natural)) || any(natural<=0))
      stop("Log-transformed graph parameters exceed numeric range.",call.=FALSE)
  }
  parameterCount <- length(graph$parameters)
  if(is.null(directions)){
    directionValues <- matrix(0,parameterCount,3L*nrow(data))
  }else if(is.matrix(directions) && is.numeric(directions) &&
      identical(dim(directions),c(parameterCount,3L))){
    if(any(!is.finite(directions)))
      stop("directions must be finite.",call.=FALSE)
    directionValues <- matrix(rep(directions,nrow(data)),
      nrow=parameterCount,ncol=3L*nrow(data))
  }else if(is.array(directions) && is.numeric(directions) &&
      identical(dim(directions),c(nrow(data),parameterCount,3L))){
    if(any(!is.finite(directions)))
      stop("directions must be finite.",call.=FALSE)
    directionValues <- matrix(0,parameterCount,3L*nrow(data))
    for(i in seq_len(nrow(data))){
      for(d in seq_len(3L))
        directionValues[,3L*(i-1L)+d] <- directions[i,,d]
    }
  }else{
    stop("directions must be a parameter-by-three matrix or observation-by-parameter-by-three array.",
      call.=FALSE)
  }
  for(j in seq_along(graph$parameterMap)){
    map <- graph$parameterMap[[j]]
    for(d in seq_len(3L)){
      directionColumn <- seq(d,by=3L,length.out=nrow(data))
      if(identical(map$type,"fixed")){
        if(any(directionValues[j,directionColumn]!=directionValues[j,directionColumn[1L]]))
          stop("Fixed parameter directions must be constant across observations.",
            call.=FALSE)
      }else{
        group <- data[[map$group]]
        for(level in unique(group)){
          rows <- which(group==level)
          selected <- 3L*(rows-1L)+d
          if(any(directionValues[j,selected]!=directionValues[j,selected[1L]]))
            stop("Random parameter directions must be constant within group '",
              map$group,"'.",call.=FALSE)
        }
      }
    }
  }
  compiled <- lapply(seq_along(graph$nodes),function(i){
    node <- graph$nodes[[i]]
    if(identical(node$op,"input"))
      node$columnIndex <- match(names(graph$nodes)[[i]],colnames(inputs))-1L
    node$inputIndices <- match(node$inputs,names(graph$nodes))-1L
    node$inputWidths <- vapply(node$inputs,function(name)
      graph$nodes[[name]]$width,integer(1))
    node
  })
  names(compiled) <- names(graph$nodes)
  parameterTransforms <- as.integer(vapply(graph$parameterMap,function(map)
    identical(map$transform,"log"),logical(1)))
  result <- sommer_nonlinear_graph_eval_cpp(compiled,
    match(graph$output,names(compiled))-1L,
    inputs,parameterValues,parameterTransforms,directionValues,
    control$execution$denseMemoryMB)
  colnames(result$gradient) <- graph$parameters
  colnames(result$hessianVector) <- graph$parameters
  colnames(result$thirdDirectionalGradient) <- graph$parameters
  result
}

.nonlinear_graph_validate <- function(graph){
  if(!inherits(graph,"sommer_nonlinear_graph") || !is.list(graph) ||
     !graph$version %in% 1:2 || is.null(graph$nodes) ||
     is.null(graph$output) || is.null(graph$parameters) ||
     is.null(graph$parameterMap))
    stop("graph must be a versioned object created by nonlinear_graph().",
      call.=FALSE)
  rebuilt <- nonlinear_graph(graph$nodes,graph$output,graph$parameters,
    graph$parameterMap,if(is.null(graph$random)) list() else graph$random,
    graph$joint)
  if(!identical(unclass(graph),unclass(rebuilt)))
    stop("The nonlinear graph layout has been modified.",call.=FALSE)
  graph
}

.nonlinear_graph_native_layout <- function(graph,data,allowMissing=FALSE){
  graph <- .nonlinear_graph_validate(graph)
  inputNodes <- which(vapply(graph$nodes,function(node)
    identical(node$op,"input"),logical(1)))
  columns <- vapply(graph$nodes[inputNodes],`[[`,character(1),"column")
  if(!length(columns) || any(!columns %in% names(data)))
    stop("Graph data is missing columns: ",
      paste(setdiff(columns,names(data)),collapse=", "),".",call.=FALSE)
  inputs <- do.call(cbind,lapply(columns,function(column){
    value <- data[[column]]
    if(!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
       length(value)!=nrow(data) || any(is.infinite(value)) ||
       (!allowMissing && anyNA(value)))
      stop("Graph input '",column,"' must contain finite numeric values.",
        call.=FALSE)
    as.numeric(value)
  }))
  if(is.null(dim(inputs))) inputs <- matrix(inputs,ncol=1L)
  compiled <- lapply(seq_along(graph$nodes),function(i){
    node <- graph$nodes[[i]]
    if(identical(node$op,"input"))
      node$columnIndex <- match(names(graph$nodes)[[i]],
        names(graph$nodes)[inputNodes])-1L
    node$inputIndices <- match(node$inputs,names(graph$nodes))-1L
    node$inputWidths <- vapply(node$inputs,function(name)
      graph$nodes[[name]]$width,integer(1))
    node
  })
  names(compiled) <- names(graph$nodes)
  list(graph=graph,nodes=compiled,
    output=match(graph$output,names(compiled))-1L,inputs=inputs,
    transformed=as.integer(vapply(graph$parameterMap,function(map)
      identical(map$transform,"log"),logical(1))))
}

.nonlinear_graph_first <- function(graph,data,parameters,
    control=mmes_control()){
  layout <- .nonlinear_graph_native_layout(graph,data)
  if(!is.numeric(parameters) || is.complex(parameters) ||
     !is.null(dim(parameters)) || is.null(names(parameters)) ||
     anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
     !setequal(names(parameters),layout$graph$parameters) ||
     any(!is.finite(parameters))){
    stop("parameters must be a finite named vector containing every graph coordinate.",
      call.=FALSE)
  }
  value <- sommer_nonlinear_graph_first_cpp(layout$nodes,layout$output,
    layout$inputs,as.numeric(parameters[layout$graph$parameters]),
    layout$transformed,control$execution$denseMemoryMB)
  if(!isTRUE(value$valid))
    stop("Nonlinear graph evaluation failed: ",value$reason,".",call.=FALSE)
  colnames(value$jacobian) <- layout$graph$parameters
  value
}

.nonlinear_graph_values <- function(graph,data,parameters,
    control=mmes_control()){
  layout <- .nonlinear_graph_native_layout(graph,data)
  if(!is.matrix(parameters) || !is.numeric(parameters) ||
     nrow(parameters)!=nrow(data) ||
     ncol(parameters)!=length(layout$graph$parameters) ||
     is.null(colnames(parameters)) ||
     !setequal(colnames(parameters),layout$graph$parameters) ||
     any(!is.finite(parameters)))
    stop("parameters must be a finite observation-by-parameter matrix with every graph coordinate.",
      call.=FALSE)
  sommer_nonlinear_graph_values_cpp(layout$nodes,layout$output,
    layout$inputs,unname(parameters[,layout$graph$parameters,drop=FALSE]),
    layout$transformed,control$execution$denseMemoryMB)
}
