.laplace_is_categorical_mixed <- function(object){
  !is.null(object$categorical)
}

.laplace_has_categorical <- function(family){
  inherits(family, "sommer_familym") &&
    any(vapply(family$families, function(component)
      component$family %in% c("ordinal", "multinomial"), logical(1)))
}

.laplace_mixed_expand <- function(data, family, components){
  index <- .laplace_family_index(family, data, nrow(data))
  dimensions <- vapply(components, function(component)
    if(component$family == "multinomial") length(component$predictorLevels) else 1L,
    integer(1))
  widths <- dimensions[index]
  widths[is.na(widths)] <- 1L
  rowMap <- rep(seq_len(nrow(data)), widths)
  coordinate <- sequence(widths)
  expanded <- data[rowMap, , drop=FALSE]
  labels <- unlist(lapply(seq_along(components), function(position){
    component <- components[[position]]
    paste(names(components)[position],
      if(component$family == "multinomial") component$predictorLevels else "location",
      sep=":")
  }), use.names=FALSE)
  values <- vapply(seq_along(rowMap), function(row){
    position <- index[rowMap[row]]
    if(is.na(position)) return(NA_character_)
    component <- components[[position]]
    paste(names(components)[position],
      if(component$family == "multinomial") component$predictorLevels[coordinate[row]]
      else "location", sep=":")
  }, character(1))
  expanded$.familyCoordinate <- factor(values, levels=labels)
  list(data=expanded, rowMap=rowMap, coordinate=coordinate)
}

.laplace_mixed_design <- function(basis, index, coordinate, map){
  design <- Matrix::Matrix(0, nrow=nrow(basis), ncol=nrow(map), sparse=TRUE)
  for(column in seq_len(nrow(map))){
    rows <- if(is.na(map$trait[column])) !index %in% map$nominalTraits[[column]] else
      index == map$trait[column] & coordinate == map$coordinate[column]
    rows[is.na(rows)] <- FALSE
    design[rows, column] <- basis[rows, map$basis[column]]
  }
  colnames(design) <- map$name
  design
}

.laplace_prepare_mixed <- function(fixed, random, rcov, data, family,
                                    naMethodY, mmesArgs, familyParameters){
  if(is.null(data) || !is.data.frame(data)){
    stop("Categorical familym() requires a data frame with numeric category codes.",
      call.=FALSE)
  }
  reserved <- c(".familyCoordinate", ".familyResponse")
  if(any(reserved %in% names(data)) || any(reserved %in% all.vars(fixed)) ||
     ".familyResponse" %in% all.vars(random)){
    stop("Categorical familym() reserves .familyCoordinate and .familyResponse; use .familyCoordinate only in explicit random-effect terms.",
      call.=FALSE)
  }
  if(".nominalCategory" %in% all.vars(random)){
    stop("Categorical familym() random expansion is ambiguous with .nominalCategory; specify shared effects or explicit .familyCoordinate covariance terms.",
      call.=FALSE)
  }
  frame <- stats::model.frame(fixed, data=data, na.action=stats::na.pass)
  response <- stats::model.response(frame)
  if(!is.numeric(response) || is.object(response) || !is.null(dim(response))){
    stop("A mixed categorical response must use plain numeric category codes (contiguous integers starting at one) and scalar observations in one long-format column.",
      call.=FALSE)
  }
  model <- .laplace_family_specs(family, familyParameters)
  index <- .laplace_family_index(family, data, nrow(data))
  components <- model$specifications
  for(position in seq_along(components)){
    component <- components[[position]]
    if(!component$family %in% c("ordinal", "multinomial")) next
    value <- response[which(index == position)]
    nominal <- .laplace_nominal_response(value, component$reference)
    component$levels <- nominal$levels
    if(component$family == "multinomial"){
      component <- c(component, nominal[setdiff(names(nominal), names(component))])
      component$reference <- nominal$reference
      response[which(index == position)] <- nominal$response
    }else if(!is.null(component$thresholds) &&
             length(component$thresholds) != length(component$levels) - 1L){
      stop("Each ordinal trait requires one fewer cutpoint than categories.", call.=FALSE)
    }
    components[[position]] <- component
  }
  dimensions <- vapply(components, function(component)
    if(component$family == "multinomial") length(component$predictorLevels) else 1L,
    integer(1))
  width <- max(dimensions)
  bytes <- 2 * as.double(width) * as.double(utils::object.size(data)) +
    8 * nrow(data) * width * 24
  if(!is.finite(bytes) || bytes > mmesArgs$engine$denseMemoryMB * 1048576 ||
     as.double(nrow(data)) * width > .Machine$integer.max){
    stop("Categorical familym() preparation exceeds execution$denseMemoryMB or sparse index limits.",
      call.=FALSE)
  }
  expanded <- .laplace_mixed_expand(data, family, components)
  if(!is.null(mmesArgs$observationMask))
    mmesArgs$observationMask <- mmesArgs$observationMask[expanded$rowMap]
  expanded$data$.familyResponse <- response[expanded$rowMap]
  numericFixed <- stats::formula(stats::terms(fixed, data=data))
  numericFixed[[2L]] <- as.name(".familyResponse")
  prepared <- .mmes_prepare_model(numericFixed, random, rcov, expanded$data,
    family, naMethodY, mmesArgs, .ordinalPrepared=TRUE)
  included <- prepared$obsInfo$included
  kept <- which(included)
  retainedOriginal <- unique(expanded$rowMap[kept])
  retainedCounts <- tabulate(expanded$rowMap[kept], nbins=nrow(data))
  fullCounts <- tabulate(expanded$rowMap, nbins=nrow(data))
  if(any(retainedCounts > 0L & retainedCounts != fullCounts)){
    stop("Categorical filtering must retain or exclude every predictor coordinate of an observation together.",
      call.=FALSE)
  }
  retainedIndex <- index[expanded$rowMap[kept]]
  basis <- prepared$X
  if(ncol(basis) && Matrix::rankMatrix(basis)[1L] != ncol(basis))
    stop("Categorical preparation requires a full-rank original fixed design.",
      call.=FALSE)
  memory <- mmesArgs$engine$denseMemoryMB * 1048576
  nominalTraits <- which(vapply(components, function(component)
    component$family == "multinomial", logical(1)))
  width <- ncol(basis) * (1 + sum(vapply(components[nominalTraits],
    function(component) length(component$predictorLevels), integer(1))))
  if(8 * as.double(nrow(basis)) * max(1, width) * 8 > memory){
    stop("Categorical familym() design expansion exceeds execution$denseMemoryMB.",
      call.=FALSE)
  }
  basis <- as.matrix(basis)
  basisModel <- stats::model.matrix(attr(prepared$modelFrame, "terms"),
    prepared$retainedFrame, contrasts.arg=mmesArgs$contrasts)
  map <- data.frame(basis=seq_len(ncol(basis)), trait=NA_integer_,
    coordinate=NA_integer_, name=colnames(basis))
  map$nominalTraits <- rep(list(nominalTraits), nrow(map))
  for(position in nominalTraits){
    for(coordinate in seq_along(components[[position]]$predictorLevels)){
      addition <- data.frame(basis=seq_len(ncol(basis)), trait=position,
        coordinate=coordinate,
        name=paste(names(components)[position],
          components[[position]]$predictorLevels[coordinate], colnames(basis), sep=":"))
      addition$nominalTraits <- rep(list(nominalTraits), nrow(addition))
      map <- rbind(map, addition)
    }
  }
  design <- .laplace_mixed_design(basis, retainedIndex,
    expanded$coordinate[kept], map)
  ordinalTraits <- which(vapply(components, function(component)
    component$family == "ordinal" && component$free, logical(1)))
  locations <- sapply(ordinalTraits, function(position)
    as.numeric(retainedIndex == position))
  if(!length(ordinalTraits)) locations <- matrix(0, nrow(design), 0)
  locations <- matrix(locations, nrow(design))
  # Put estimated cutpoint locations first in the rank selection, so only
  # fixed columns aliased with these locations are removed.
  augmented <- cbind(locations, as.matrix(design))
  rank <- qr(augmented, LAPACK=FALSE, tol=1e-9)
  selected <- sort(rank$pivot[seq_len(rank$rank)])
  columns <- selected[selected > ncol(locations)] - ncol(locations)
  design <- design[, columns, drop=FALSE]
  map <- map[columns, , drop=FALSE]
  for(position in seq_along(components)){
    component <- components[[position]]
    rows <- which(retainedIndex == position)
    if(!length(rows)) stop("Each declared family must have retained observations.", call.=FALSE)
    if(component$family %in% c("ordinal", "multinomial")){
      observed <- expanded$data$.familyResponse[kept[rows]]
      counts <- tabulate(observed, nbins=length(component$levels))
      if((component$free || component$family == "multinomial") && any(counts == 0L)){
        stop("Estimated categorical parameters require retained observations in every category.",
          call.=FALSE)
      }
      component$counts <- counts
      if(component$family == "ordinal" && is.null(component$thresholds)){
        component$thresholds <- model$families[[position]]$linkfun(
          head(cumsum(counts) / sum(counts), -1L))
      }
      if(component$family == "multinomial"){
        component$nominalGroups <- matrix(rows,
          ncol=length(component$predictorLevels), byrow=TRUE)
      }
    }
    components[[position]] <- component
  }
  prepared$X <- design
  prepared$categorical <- list(components=components, map=map,
    basisColumns=colnames(basis), expandedData=expanded$data,
    contrasts=attr(basisModel, "contrasts"),
    assignment=attr(basisModel, "assign")[map$basis],
    numericFixed=numericFixed, originalFixed=fixed,
    originalData=data[retainedOriginal, , drop=FALSE],
    originalRows=retainedOriginal,
    rowMap=match(expanded$rowMap[kept], retainedOriginal),
    coordinate=expanded$coordinate[kept],
    originalIndex=index[retainedOriginal], included=included,
    originalObsInfo=prepared$obsInfo[match(seq_len(nrow(data)), expanded$rowMap), , drop=FALSE])
  prepared
}

.laplace_mixed_output <- function(object, eta, type="response", index=NULL,
                                   rowMap=NULL, coordinate=NULL){
  metadata <- object$categorical
  if(is.null(index)) index <- object$laplaceFamilyIndex
  if(is.null(rowMap)) rowMap <- metadata$rowMap
  if(is.null(coordinate)) coordinate <- metadata$coordinate
  specifications <- object$familySpecifications
  result <- lapply(seq_along(specifications), function(position){
    rows <- which(index == position)
    specification <- specifications[[position]]
    if(!length(rows)){
      if(specification$family == "multinomial" ||
         (specification$family == "ordinal" && type == "response")){
        labels <- if(type == "link") specification$predictorLevels else specification$levels
        return(matrix(numeric(), 0L, length(labels), dimnames=list(NULL, labels)))
      }
      return(numeric())
    }
    if(specification$family == "multinomial"){
      originals <- unique(rowMap[rows])
      logits <- matrix(eta[rows], nrow=length(originals), byrow=TRUE,
        dimnames=list(as.character(originals), specification$predictorLevels))
      if(type == "link") return(logits)
      return(.laplace_nominal_probabilities(logits, specification))
    }
    if(specification$family == "ordinal"){
      if(type == "link") return(eta[rows])
      return(.laplace_ordinal_distribution(eta[rows], specification,
        specification$levels, object$laplacePostPEVSettings$denseMemoryMB)$probability)
    }
    if(type == "link") return(eta[rows])
    object$laplaceFamilies[[position]]$linkinv(eta[rows])
  })
  names(result) <- names(specifications)
  result
}

.laplace_numeric_gradient <- function(parameters, objective){
  vapply(seq_along(parameters), function(position){
    step <- 1e-5 * max(1, abs(parameters[position]))
    upper <- lower <- parameters
    upper[position] <- upper[position] + step
    lower[position] <- lower[position] - step
    (objective(upper) - objective(lower)) / (2 * step)
  }, numeric(1))
}

.laplace_mixed_residuals <- function(object, type){
  gaussianCovariance <- if(!is.null(object$laplaceResidual))
    .laplace_residual_matrix(object) else NULL
  gaussianStandardized <- if(!is.null(gaussianCovariance))
    as.numeric(forwardsolve(t(chol(gaussianCovariance)),
      object$y[object$laplaceResidual$rows]-
        object$linear.predictors[object$laplaceResidual$rows])) else NULL
  result <- lapply(seq_along(object$familySpecifications), function(position){
    specification <- object$familySpecifications[[position]]
    rows <- which(object$laplaceFamilyIndex == position)
    local <- object
    local$laplaceResidual <- NULL
    if(specification$family=="gaussian" && !is.null(gaussianCovariance)){
      selected <- match(rows,object$laplaceResidual$rows)
      if(type=="deviance") return(gaussianStandardized[selected])
      local$laplaceResidual <- list(rows=seq_along(rows))
      local$laplaceResidualCovariance <- gaussianCovariance[selected,selected,drop=FALSE]
    }
    local$categorical <- NULL
    local$family <- object$laplaceFamilies[[position]]
    local$laplaceFamilies <- list(local$family)
    local$familySpecifications <- list(specification)
    local$laplaceFamilyIndex <- rep(1L, length(rows))
    local$y <- object$y[rows]
    local$linear.predictors <- object$linear.predictors[rows]
    local$fitted.values <- object$fitted.values[rows]
    local$priorWeights <- object$priorWeights[rows]
    local$dispersion <- object$dispersion[position]
    if(specification$family == "ordinal"){
      local$ordinal <- list(levels=specification$levels)
    }else if(specification$family == "multinomial"){
      local$nominal <- specification
      local$linear.predictors <- matrix(local$linear.predictors,
        ncol=length(specification$predictorLevels), byrow=TRUE)
      local$category.probabilities <- .laplace_nominal_probabilities(
        local$linear.predictors, specification)
      response <- local$y[seq.int(1L, length(local$y),
        by=length(specification$predictorLevels))]
      local$y <- factor(specification$nativeLevels[response], levels=specification$levels)
    }else if(!is.null(specification$zeroDesign)){
      local$familySpecifications[[1L]]$zeroDesign <- specification$zeroDesign[rows, , drop=FALSE]
      local$familySpecifications[[1L]]$zeroOffset <- specification$zeroOffset[rows]
    }
    residuals(local, type=type)
  })
  names(result) <- names(object$familySpecifications)
  result
}

.laplace_mixed_prediction_design <- function(object, newdata, includeRandom=TRUE,
                                               population=FALSE){
  metadata <- object$categorical
  if(is.null(newdata)){
    return(list(fixed=object$fixedDesign, random=if(includeRandom || population)
      object$randomDesign else object$randomDesign * 0, offset=object$offset,
      index=object$laplaceFamilyIndex, rowMap=metadata$rowMap,
      coordinate=metadata$coordinate))
  }
  if(!is.data.frame(newdata) || !nrow(newdata)){
    stop("Categorical familym() prediction requires nonempty newdata.", call.=FALSE)
  }
  width <- max(vapply(object$familySpecifications, function(specification)
    if(specification$family == "multinomial") length(specification$predictorLevels) else 1L,
    integer(1)))
  bytes <- 2 * as.double(width) * as.double(utils::object.size(newdata)) +
    8 * nrow(newdata) * width * (24 + ncol(object$fixedDesign))
  if(!is.finite(bytes) || bytes > object$laplacePostPEVSettings$denseMemoryMB * 1048576)
    stop("Categorical prediction preparation exceeds execution$denseMemoryMB.", call.=FALSE)
  expanded <- .laplace_mixed_expand(newdata, object$family, object$familySpecifications)
  frame <- stats::model.frame(object$fixedTerms, data=expanded$data,
    na.action=stats::na.pass, xlev=object$fixedXlevels)
  basis <- stats::model.matrix(object$fixedTerms, frame, contrasts.arg=object$contrasts)
  if(!setequal(colnames(basis), metadata$basisColumns) || any(!is.finite(basis))){
    stop("Categorical prediction fixed design differs from the fitted model or contains missing/nonfinite predictors.",
      call.=FALSE)
  }
  basis <- basis[, metadata$basisColumns, drop=FALSE]
  index <- .laplace_family_index(object$family, expanded$data, nrow(expanded$data))
  if(anyNA(index)) stop("Every prediction row must name a declared trait.", call.=FALSE)
  fixed <- .laplace_mixed_design(basis, index, expanded$coordinate, metadata$map)
  offset <- stats::model.offset(frame)
  if(is.null(offset)) offset <- numeric(nrow(frame))
  if(any(!is.finite(offset))) stop("Prediction offsets must be finite.", call.=FALSE)
  random <- object$randomDesign[rep(1L, nrow(fixed)), , drop=FALSE] * 0
  if(includeRandom || population){
    stacked <- expanded$data
    for(variable in setdiff(names(object$randomFactorLevels), ".familyCoordinate")){
      if(!variable %in% names(stacked)) next
      values <- as.character(stacked[[variable]])
      known <- object$randomFactorLevels[[variable]]
      unseen <- setdiff(unique(values), known)
      if(length(unseen)){
        if(!population){
          stop("Conditional categorical prediction encountered unseen random-effect levels; use re.form=NA or supported population prediction.",
            call.=FALSE)
        }
        .laplace_population_new_levels(object,variable)
        values[values %in% unseen] <- known[1L]
      }
      stacked[[variable]] <- factor(values, levels=known)
    }
    stacked$.familyResponse <- 1
    training <- metadata$expandedData[metadata$included, , drop=FALSE]
    training$.familyResponse <- 1
    bytes <- 2 * as.double(utils::object.size(stacked)) +
      2 * as.double(utils::object.size(training)) +
      64 * (nrow(stacked) + nrow(training))
    if(bytes > object$laplacePostPEVSettings$denseMemoryMB * 1048576)
      stop("Categorical random prediction preparation exceeds execution$denseMemoryMB.",
        call.=FALSE)
    common <- intersect(names(stacked), names(training))
    stacked <- rbind(stacked[, common, drop=FALSE], training[, common, drop=FALSE])
    randomFixed <- metadata$numericFixed
    randomFixed[[3L]] <- 1
    randomArguments <- object$predictionSettings$mmesArgs
    randomArguments$contrasts <- NULL
    prepared <- .mmes_prepare_model(randomFixed, object$randomFormula,
      NULL, stacked, stats::gaussian(), object$predictionSettings$naMethodY,
      randomArguments, .ordinalPrepared=TRUE)
    if(!all(prepared$obsInfo$included) ||
       length(prepared$Z) != length(object$randomBlockColumns)){
      stop("Categorical prediction requires complete, unambiguous random-effect designs.",
        call.=FALSE)
    }
    blocks <- lapply(seq_along(prepared$Z), function(term){
      block <- prepared$Z[[term]][seq_len(nrow(fixed)), , drop=FALSE]
      columns <- object$randomBlockColumns[[term]]
      if(is.null(columns) || is.null(colnames(block))){
        if(ncol(block) != object$randomDesignBlockWidths[term])
          stop("Categorical random-effect expansion is ambiguous; specify shared effects or explicit .familyCoordinate covariance terms.",
            call.=FALSE)
        return(block)
      }
      positions <- match(colnames(block), columns)
      entries <- Matrix::summary(block)
      if(anyNA(positions[entries$j]))
        stop("Categorical prediction random-effect columns are absent from the fitted model.",
          call.=FALSE)
      Matrix::sparseMatrix(i=entries$i, j=positions[entries$j], x=entries$x,
        dims=c(nrow(block), length(columns)), dimnames=list(NULL, columns))
    })
    random <- Matrix::Matrix(do.call(cbind, blocks), sparse=TRUE)
  }
  list(fixed=fixed, random=random, offset=offset, index=index,
    rowMap=expanded$rowMap, coordinate=expanded$coordinate)
}

.laplace_mixed_population_output <- function(object, prediction, eta){
  result <- .laplace_mixed_output(object, eta, index=prediction$index,
    rowMap=prediction$rowMap, coordinate=prediction$coordinate)
  variance <- .laplace_mixed_population_variance(object, prediction$random)
  for(position in seq_along(object$familySpecifications)){
    specification <- object$familySpecifications[[position]]
    rows <- which(prediction$index == position)
    if(!length(rows)) next
    local <- object
    local$categorical <- NULL
    local$family <- object$laplaceFamilies[[position]]
    local$laplaceFamilies <- list(local$family)
    local$familySpecifications <- list(specification)
    if(specification$family == "multinomial"){
      d <- length(specification$predictorLevels)
      n <- length(rows)/d
      ordering <- as.vector(matrix(rows, ncol=d, byrow=TRUE))
      local$nominal <- specification
      result[[position]] <- .laplace_nominal_population_prediction(local,
        list(n=n, fixed=prediction$fixed[ordering, , drop=FALSE],
          random=prediction$random[ordering, , drop=FALSE],
          offset=prediction$offset[ordering]))
    }else if(specification$family == "ordinal"){
      local$ordinal <- list(levels=specification$levels)
      result[[position]] <- .laplace_ordinal_output(local, eta[rows], "response",
        randomVariance=variance[rows])
    }else{
      result[[position]] <- .laplace_population_mean(local, eta[rows],
        rep(1L, length(rows)), variance[rows])
    }
  }
  result
}

.laplace_mixed_population_variance <- function(object, design){
  .laplace_population_variance(design,object$covStruct,object$randomRelationships,
    object$randomBlockWidths,object$laplacePostPEVSettings$denseMemoryMB)
}

.laplace_mixed_predict <- function(object, arguments){
  allowed <- c("newdata", "type", "re.form", "se.fit", "cov.fit", "uncertainty")
  if(length(arguments) && (is.null(names(arguments)) || anyDuplicated(names(arguments)) ||
      anyNA(names(arguments)) || any(!names(arguments) %in% allowed))){
    stop("Categorical familym() prediction supports newdata, type, re.form, se.fit, cov.fit and uncertainty.",
      call.=FALSE)
  }
  type <- if(is.null(arguments$type)) "response" else
    match.arg(arguments$type, c("response", "probabilities", "link"))
  if(type == "probabilities") type <- "response"
  se <- if(is.null(arguments$se.fit)) FALSE else .mme_flag(arguments$se.fit, "se.fit")
  covariance <- if(is.null(arguments$cov.fit)) FALSE else .mme_flag(arguments$cov.fit, "cov.fit")
  uncertainty <- if(is.null(arguments$uncertainty)) "conditional" else
    match.arg(arguments$uncertainty, c("conditional", "full"))
  if(!se && (covariance || uncertainty == "full"))
    stop("cov.fit=TRUE and uncertainty='full' require se.fit=TRUE.", call.=FALSE)
  randomForm <- arguments$re.form
  population <- identical(randomForm, "population")
  if(!is.null(randomForm) && !identical(randomForm, NA) && !population)
    stop("re.form must be NULL, NA or 'population'.", call.=FALSE)
  if(population && type != "response"){
    stop("Mixed categorical population prediction supports response probabilities/means only.",
      call.=FALSE)
  }
  prediction <- .laplace_mixed_prediction_design(object, arguments$newdata,
    includeRandom=is.null(randomForm), population=population)
  if(!is.null(arguments$newdata)){
    expanded <- .laplace_mixed_expand(arguments$newdata, object$family,
      object$familySpecifications)
    object$familySpecifications <- .laplace_prediction_specs(object,
      expanded$data, prediction$index)
    object$laplaceFamilies <- .laplace_update_families(object$laplaceFamilies,
      object$familySpecifications, prediction$index)
  }
  eta <- prediction$offset + as.numeric(prediction$fixed %*% object$b)
  if(is.null(randomForm)) eta <- eta + as.numeric(prediction$random %*% object$randomMode)
  if(population){
    value <- .laplace_mixed_population_output(object,prediction,eta)
    if(!se) return(value)
    return(.laplace_population_delta(object,value,function(local)
      .laplace_mixed_population_output(local,prediction,
        prediction$offset+as.numeric(prediction$fixed %*% local$b)),covariance))
  }
  result <- .laplace_mixed_output(object, eta, type, prediction$index,
    prediction$rowMap, prediction$coordinate)
  if(!se) return(result)
  .laplace_mixed_prediction_uncertainty(object, prediction, result, type,
    uncertainty, covariance)
}

.laplace_mixed_output_derivative <- function(object, prediction, eta, type){
  result <- .laplace_mixed_output(object, eta, type, prediction$index,
    prediction$rowMap, prediction$coordinate)
  total <- sum(lengths(result))
  derivative <- Matrix::Matrix(0, nrow=total, ncol=length(eta), sparse=TRUE)
  start <- 0L
  for(position in seq_along(result)){
    rows <- which(prediction$index == position)
    specification <- object$familySpecifications[[position]]
    if(!length(rows)) next
    n <- length(rows)
    if(specification$family == "multinomial"){
      d <- length(specification$predictorLevels)
      n <- n/d
      groups <- matrix(rows, ncol=d, byrow=TRUE)
      if(type == "link"){
        derivative[cbind(start + seq_len(n*d), as.vector(groups))] <- 1
      }else{
        p <- result[[position]]
        for(category in seq_len(ncol(p))) for(coordinate in seq_len(d)){
          outputRows <- start + (category - 1L)*n + seq_len(n)
          derivative[cbind(outputRows, groups[, coordinate])] <-
            p[, category] * (as.numeric(specification$levels[category] ==
              specification$predictorLevels[coordinate]) -
                p[, match(specification$predictorLevels[coordinate], specification$levels)])
        }
      }
    }else if(specification$family == "ordinal" && type == "response"){
      distribution <- .laplace_ordinal_distribution(eta[rows], specification,
        specification$levels, object$laplacePostPEVSettings$denseMemoryMB)
      derivative[cbind(start + seq_along(distribution$derivative),
        rep(rows, ncol(distribution$derivative)))] <- as.numeric(distribution$derivative)
    }else{
      derivative[cbind(start + seq_len(n), rows)] <- if(type == "link") 1 else
        object$laplaceFamilies[[position]]$mu.eta(eta[rows])
    }
    start <- start + length(result[[position]])
  }
  derivative
}

.laplace_mixed_prediction_uncertainty <- function(object, prediction, fitted,
                                                   type, mode, covariance){
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  training <- object
  training$familySpecifications <- object$categorical$components
  uncertainty <- .mmes_laplace_uncertainty(training)
  if(!uncertainty$valid)
    stop("Mixed categorical prediction uncertainty is unavailable: ",
      uncertainty$reason, ".", call.=FALSE)
  parameters <- object$optimizer$par
  outputCount <- sum(lengths(fitted))
  q <- ncol(prediction$random)
  memory <- object$laplacePostPEVSettings$denseMemoryMB * 1048576
  bytes <- 8 * (6 * as.double(outputCount)^2 +
    4 * as.double(outputCount) * (length(parameters) + q))
  if(!is.finite(bytes) || bytes > memory)
    stop("Mixed categorical prediction uncertainty exceeds execution$denseMemoryMB.",
      call.=FALSE)
  context <- .laplace_objective_context(training,
    memoryMB=object$laplacePostPEVSettings$denseMemoryMB)
  evaluate <- function(proposal){
    local <- object
    local$familySpecifications <- .laplace_proposal_specs(object$familySpecifications,
      proposal, object$laplaceParameterLayout)
    local$laplaceFamilies <- .laplace_update_families(object$laplaceFamilies,
      local$familySpecifications, prediction$index)
    value <- .laplace_constrained_evaluate(proposal, object$laplaceParameterLayout,
      function(expanded) sommer_laplace_objective_eval(context, expanded, gradient=FALSE,
        maxModeIterations=object$laplacePostPEVSettings$modeMaxit,
        tolerance=object$laplacePostPEVSettings$modeTol))
    if(!isTRUE(value$valid))
      stop("Mixed prediction uncertainty encountered invalid parameters: ", value$reason,
        call.=FALSE)
    eta <- prediction$offset + as.numeric(prediction$fixed %*%
      proposal[seq_len(length(object$b))]) +
      as.numeric(prediction$random %*% value$mode)
    as.numeric(unlist(.laplace_mixed_output(local, eta, type, prediction$index,
      prediction$rowMap, prediction$coordinate), use.names=FALSE))
  }
  selected <- seq_along(parameters)
  parameterCovariance <- uncertainty$parameterCovariance
  if(mode == "conditional"){
    counts <- vapply(object$familySpecifications, .laplace_family_count, integer(1))
    selected <- seq_len(length(object$b))
    for(position in seq_along(counts)){
      if(object$familySpecifications[[position]]$family != "ordinal") next
      selected <- c(selected, length(object$b) +
        object$laplaceParameterLayout$covarianceCount +
        sum(head(counts, position-1L)) + seq_len(counts[position]))
    }
    parameterCovariance <- solve(uncertainty$hessian[selected, selected, drop=FALSE])
  }
  jacobian <- matrix(0, outputCount, length(selected))
  for(column in seq_along(selected)){
    position <- selected[column]
    step <- 1e-5 * max(1, abs(parameters[position]))
    upper <- lower <- parameters
    upper[position] <- upper[position] + step
    lower[position] <- lower[position] - step
    jacobian[, column] <- (evaluate(upper) - evaluate(lower))/(2*step)
  }
  resultCovariance <- jacobian %*% parameterCovariance %*% t(jacobian)
  eta <- prediction$offset + as.numeric(prediction$fixed %*% object$b) +
    as.numeric(prediction$random %*% object$randomMode)
  derivative <- .laplace_mixed_output_derivative(object, prediction, eta, type)
  randomDerivative <- derivative %*% prediction$random
  solveMode <- .laplace_mode_solver(training)
  solved <- solveMode(Matrix::t(randomDerivative))
  resultCovariance <- resultCovariance + as.matrix(randomDerivative %*% solved)
  resultCovariance <- (resultCovariance + t(resultCovariance))/2
  variances <- diag(resultCovariance)
  if(any(!is.finite(resultCovariance)) || any(variances < 0))
    stop("Mixed prediction covariance is invalid; no uncertainty was returned.", call.=FALSE)
  standardErrors <- fitted
  start <- 0L
  for(position in seq_along(fitted)){
    values <- sqrt(variances[start + seq_along(fitted[[position]])])
    if(is.matrix(fitted[[position]])){
      standardErrors[[position]][] <- values
    }else standardErrors[[position]] <- values
    start <- start + length(values)
  }
  labels <- .laplace_prediction_labels(fitted)
  dimnames(resultCovariance) <- list(labels, labels)
  result <- list(fit=fitted, se.fit=standardErrors, uncertainty=mode)
  if(covariance) result$cov.fit <- resultCovariance
  result
}

.laplace_mixed_simulate <- function(object, nsim, includeRandomEffects,
                                     maxMemoryMB, randomFactors=NULL){
  if(includeRandomEffects && is.null(randomFactors))
    randomFactors <- .laplace_random_effect_sampling_factors(object, maxMemoryMB)
  simulations <- matrix(NA_real_, object$nobs, nsim)
  residualFactor <- if(!is.null(object$laplaceResidual))
    chol(.laplace_residual_matrix(object,maxMemoryMB)) else NULL
  for(simulation in seq_len(nsim)){
    eta <- object$linear.predictors
    if(includeRandomEffects){
      eta <- object$offset + as.numeric(object$fixedDesign %*% object$b) +
        as.numeric(object$randomDesign %*% .laplace_draw_random_effects(randomFactors))
    }
    output <- .laplace_mixed_output(object, eta)
    for(position in seq_along(output)){
      rows <- which(object$laplaceFamilyIndex == position)
      originals <- unique(object$categorical$rowMap[rows])
      specification <- object$familySpecifications[[position]]
      if(!is.null(residualFactor) && specification$family=="gaussian") next
      if(specification$family %in% c("ordinal", "multinomial")){
        p <- output[[position]]
        uniform <- stats::runif(nrow(p))
        cumulative <- p[, 1L]
        draw <- rep(1L, nrow(p))
        for(category in 2:ncol(p)){
          draw <- draw + as.integer(uniform > cumulative)
          cumulative <- cumulative + p[, category]
        }
      }else{
        local <- object
        local$laplaceResidual <- NULL
        local$categorical <- NULL
        local$family <- object$laplaceFamilies[[position]]
        local$laplaceFamilies <- list(local$family)
        if(!is.null(specification$zeroDesign)){
          specification$zeroDesign <- specification$zeroDesign[rows, , drop=FALSE]
          specification$zeroOffset <- specification$zeroOffset[rows]
        }
        local$familySpecifications <- list(specification)
        local$laplaceFamilyIndex <- rep(1L, length(rows))
        local$nobs <- length(rows)
        local$linear.predictors <- eta[rows]
        local$fitted.values <- output[[position]]
        local$priorWeights <- object$priorWeights[rows]
        draw <- as.numeric(.laplace_simulate_response(local, 1L, maxMemoryMB=maxMemoryMB))
      }
      simulations[originals, simulation] <- draw
    }
    if(!is.null(residualFactor)){
      rows <- object$laplaceResidual$rows
      simulations[object$categorical$rowMap[rows],simulation] <- eta[rows]+
        as.numeric(t(residualFactor) %*% stats::rnorm(length(rows)))
    }
  }
  simulations
}
