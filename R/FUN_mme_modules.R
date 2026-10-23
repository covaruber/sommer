.mme_module <- function(value, defaults, name){
  if(!is.list(value) || (length(value) &&
     (is.null(names(value)) || anyNA(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value))))){
    stop(name, " must be a list with unique, non-empty names.", call.=FALSE)
  }
  unknown <- setdiff(names(value), names(defaults))
  if(length(unknown)){
    stop("Unknown ", name, " settings: ", paste(unknown, collapse=", "), ".", call.=FALSE)
  }
  utils::modifyList(defaults, value, keep.null=TRUE)
}

.mme_choice <- function(value, choices, name){
  if(!is.character(value) || length(value) != 1L || is.na(value) ||
     !tolower(value) %in% choices){
    stop(name, " must be one of ", paste(shQuote(choices), collapse=", "), ".", call.=FALSE)
  }
  tolower(value)
}

.mme_flag <- function(value, name){
  if(!is.logical(value) || length(value) != 1L || is.na(value)){
    stop(name, " must be a single TRUE/FALSE value.", call.=FALSE)
  }
  value
}

.mme_positive <- function(value, name, integer=FALSE){
  if(!is.numeric(value) || is.complex(value) || length(value) != 1L ||
     !is.finite(value) || value <= 0 ||
     (integer && (value > .Machine$integer.max || value != floor(value)))){
    stop(name, " must be one positive ", if(integer) "integer." else "finite number.", call.=FALSE)
  }
  if(integer) as.integer(value) else value
}

.mme_blend_weight <- function(value){
  if(is.null(value)) return(value)
  if(!is.numeric(value) || is.complex(value) || !length(value)){
    stop("update$blendWeight must be a non-empty numeric scalar or vector.", call.=FALSE)
  }
  if(any(!is.finite(value)) || any(value < 0)){
    stop("update$blendWeight must contain finite non-negative values.", call.=FALSE)
  }
  if(length(value) > 1L && any(value > 1)){
    stop("update$blendWeight vectors must contain values between 0 and 1.", call.=FALSE)
  }
  if(length(value) == 1L && value > 1 && value != floor(value)){
    stop("update$blendWeight above 1 must be a whole-number decay iteration count.", call.=FALSE)
  }
  value
}

.mme_mode_paths <- c("representation$sparseFactors", "representation$matrixFree",
  "solve$blockSchur", "solve$blockBanded", "solve$factorSchurMultiFactor",
  "solve$dFirstML", "score$covarianceChecks", "score$traceBatch",
  "curvature$reduced", "update$commonScale", "output$selectedInverse",
  "diagnostics$profile")

.mme_pcg_paths <- c(setdiff(unname(.mme_pcg_mapping()),
  c("representation$compression", "execution$seed")),
  "representation$matrixFree", "execution$seed")

.mme_pcg_solve_names <- names(.mme_pcg_mapping())[
  startsWith(.mme_pcg_mapping(), "solve$")]

.mme_pql_defaults <- list(maxit=20L, tol=1e-5, estimateTheta=FALSE, warmStart=TRUE)

.mme_control_paths <- function(control){
  unlist(lapply(setdiff(names(control), "formulation"), function(module)
    paste0(module, "$", names(control[[module]]))), use.names=FALSE)
}

.mme_path_value <- function(control, path){
  parts <- strsplit(path, "$", fixed=TRUE)[[1L]]
  Reduce(function(x, part) x[[part]], parts, control)
}

# A setting is requested when the user supplied it and it changes the
# default method; experimental toggles only request work when forced.
.mme_requested <- function(control, defaults, supplied, path){
  if(!path %in% supplied) return(FALSE)
  value <- .mme_path_value(control, path)
  if(path %in% .mme_mode_paths) return(identical(value, "force"))
  if(path == "update$blendWeight") return(!is.null(value) && any(value > 0))
  !identical(value, .mme_path_value(defaults, path))
}

# Single compatibility table for model-independent option combinations.
# Each rule names a context and the settings that are either unused or
# unsupported there; model-dependent eligibility is checked during fitting.
.mmes_compatibility_rules <- function(){
  accelerated <- function(x) x$proposal$extrapolation != "none" ||
    x$curvature$method != "exact-ai"
  emLinked <- c("update$recovery", "update$blendWeight", "update$commonScale",
    "curvature$reduced")
  list(
    list(id="direct-engine",
      when=function(x) x$formulation == "direct" && x$task$mode == "fit",
      message="formulation='direct' inverts the observation covariance V directly; Henderson solvers, representations, curvature variants, PCG estimators and Henderson outputs do not apply.",
      settings=function(x) setdiff(.mme_control_paths(x), c("task$mode",
        "task$objective", "task$pql", "update$recovery", "update$blendWeight",
        "proposal$stepWeight", "convergence$maxIters", "convergence$likelihoodTol",
        "convergence$parameterTol", "output$setup",
        "representation$fixedRank",
        "execution$denseMemoryMB", "execution$memoryPolicy",
        "execution$blasThreads", "execution$inverseTol",
        paste0("diagnostics$", names(x$diagnostics)))),
      bad=list(`output$inverse`=function(value) value == 1L),
      hint="Remove them or use formulation='henderson'; the direct engine supports output$inverse 0 or 2."),
    list(id="solve-mode",
      when=function(x) x$task$mode == "solve",
      message="task$mode='solve' runs one Henderson PCG solve at the covariance parameters supplied in covPar; formulation='direct', estimation modules and non-PCG solvers are not used.",
      settings=function(x) c("formulation", "task$objective", "task$pql",
        "solve$method", "solve$nystromRank", "solve$recycleRank",
        "solve$recycleMethod", "solve$recycleState", "solve$exportRecycle",
        "solve$rhsMethod", "solve$blockSize", "solve$preconditioner",
        grep("^(representation|logdet|score|curvature|update|proposal)\\$",
          .mme_control_paths(x), value=TRUE),
        "convergence$maxIters", "convergence$likelihoodTol",
        "convergence$parameterTol", "convergence$accuracy",
        "convergence$maxTraceProbes", "convergence$maxLanczosSteps",
        "convergence$logDetTol", "convergence$scoreTol", "output$inverse",
        "output$covariance", "output$coefficientMatrix",
        "output$selectedInverse", "output$panelSize", "execution$gpu",
        "execution$gpuMinSize"),
      bad=list(`solve$method`=function(value) !value %in% c("auto", "pcg"),
        `solve$preconditioner`=function(value) value == "nystrom"),
      hint="Remove them, or use task$mode='fit' to estimate covariance parameters."),
    list(id="deterministic-solver",
      when=function(x) x$formulation == "henderson" && x$task$mode == "fit" &&
        x$solve$method != "pcg",
      message="PCG preconditioning, recycling, stochastic log-determinant and trace-estimation settings only apply when solve$method='pcg' ('auto' selects 'ldlt' or 'cholmod').",
      settings=function(x) .mme_pcg_paths,
      hint="Remove them or set solve$method='pcg'."),
    list(id="pcg-estimation",
      when=function(x) x$formulation == "henderson" && x$task$mode == "fit" &&
        x$solve$method == "pcg",
      message="solve$method='pcg' uses stochastic AI-REML; EM/PXEM recovery, blending, common-scale profiling, reduced AI, Aitken, hybrid AI, saddle-point equations, factor scores and ML need a deterministic factorization.",
      settings=function(x) c(emLinked, "proposal$extrapolation", "curvature$method",
        "representation$equations", "representation$factorScores", "task$objective"),
      hint="Remove them or use solve$method 'auto', 'ldlt' or 'cholmod'."),
    list(id="accelerated-proposal",
      when=function(x) x$formulation == "henderson" && x$task$mode == "fit" &&
        accelerated(x),
      message="Aitken extrapolation and hybrid AI replace the EM/PXEM recovery path; recovery, blending, common-scale profiling, reduced AI and saddle-point equations are not combined with them.",
      settings=function(x) c(emLinked, "representation$equations",
        if(x$curvature$method == "hybrid-ai") "representation$factorScores"),
      hint="Remove them or use proposal$extrapolation='none' with curvature$method='exact-ai'."),
    list(id="saddle-point",
      when=function(x) x$formulation == "henderson" && x$task$mode == "fit" &&
        x$representation$equations == "saddle-point",
      message="representation$equations='saddle-point' uses its own deterministic AI path without EM/PXEM recovery, blending, common-scale profiling, reduced AI or factor scores.",
      settings=function(x) c(emLinked, "representation$factorScores"),
      hint="Remove them or use representation$equations='normal'."),
    list(id="factor-scores",
      when=function(x) x$formulation == "henderson" && x$task$mode == "fit" &&
        x$representation$factorScores != "none",
      message="representation$factorScores augments the Henderson equations for REML without C inverse, EM/PXEM recovery, compression or reduced covariance outputs.",
      settings=function(x) c(emLinked, "task$objective", "output$inverse",
        "output$covariance", "output$coefficientMatrix", "representation$compression"),
      hint="Remove them or use representation$factorScores='none'."),
    list(id="blending",
      when=function(x) !is.null(x$update$blendWeight) && any(x$update$blendWeight > 0),
      message="Positive update$blendWeight blends EM/PXEM directions and cannot be combined with forced common-scale profiling.",
      settings=function(x) "update$commonScale",
      hint="Remove update$commonScale or use update$blendWeight=0."),
    list(id="coefficient-matrix",
      when=function(x) x$update$method != "conditional-fa" &&
        x$representation$equations != "conditional" &&
        (x$representation$equations != "saddle-point" || x$output$inverse != 0L),
      message="output$coefficientMatrix=FALSE requires representation$equations='saddle-point' and output$inverse=0.",
      settings=function(x) "output$coefficientMatrix",
      hint="Remove it or set those two settings."),
    list(id="conditional-fa",
      when=function(x) x$update$method == "conditional-fa" ||
        x$representation$equations == "conditional",
      message="Conditional FA uses sparse latent equations with AI or EM updates, without PCG, augmentation, rotation, covariance uncertainty or a retained public coefficient matrix.",
      settings=function(x) setdiff(.mme_control_paths(x), c(
        "task$mode", "task$objective", "update$method", "solve$method",
        "representation$equations",
        "convergence$maxIters", "convergence$likelihoodTol", "convergence$parameterTol",
        "output$covariance", "output$coefficientMatrix", "output$inverse",
        "output$panelSize", "output$setup", "execution$denseMemoryMB",
        if(x$update$method == "ai") "execution$inverseTol",
        "execution$blasThreads", "diagnostics$verbose", "diagnostics$dateWarning")),
      hint="Remove unrelated method requests or use update$method='ai'.")
  )
}

.mmes_check_compatibility <- function(control, supplied,
                                      defaults=.mmes_control_build()){
  problems <- character()
  for(rule in .mmes_compatibility_rules()){
    if(!isTRUE(rule$when(control))) next
    paths <- unique(rule$settings(control))
    offending <- paths[vapply(paths, function(path){
      if(!is.null(rule$bad[[path]])){
        return(path %in% supplied && isTRUE(rule$bad[[path]](.mme_path_value(control, path))))
      }
      if(path == "formulation"){
        return(path %in% supplied && control$formulation != defaults$formulation)
      }
      .mme_requested(control, defaults, supplied, path)
    }, logical(1))]
    if(length(offending)){
      problems <- c(problems, paste0("* ", rule$message, " Offending: ",
        paste(offending, collapse=", "), ". ", rule$hint))
    }
  }
  if(length(problems)){
    stop("Incompatible mmes_control() settings:\n", paste(problems, collapse="\n"),
         call.=FALSE)
  }
  invisible(TRUE)
}

#' Modular Mixed Model Numerical Controls
#'
#' Assigns numerical requests to mathematical modules. Model-dependent eligibility
#' is checked during fitting. No new combinations of numerical methods are implied.
#' @param formulation Equation domain: \code{"henderson"} or \code{"direct"}.
#' @param task Named list selecting \code{mode}, \code{objective} and outer \code{pql} controls.
#' @param representation Named list of equation, compression and covariance representation requests.
#' @param solve Named list of linear solver, preconditioning, reuse and elimination requests.
#' @param logdet Named list controlling stochastic determinant scaling and Lanczos depth.
#' @param score Named list controlling trace estimation and covariance checks.
#' @param curvature Named list selecting exact/hybrid AI and reduced-direction assembly.
#' @param update Named list controlling EM/PXEM recovery, blending and common-scale profiling.
#' @param proposal Named list controlling step scaling and Aitken extrapolation.
#' @param convergence Named list of inner iteration, convergence and stochastic accuracy controls.
#' @param output Named list of coefficient/covariance uncertainty and setup output requests.
#'   For Henderson fits, \code{inverse=0} skips fit-time PEVs and final inverse
#'   work; \code{inverse=1} computes random-effect PEVs without the full inverse;
#'   \code{inverse=2} computes PEVs and the full inverse. With mode 0, use
#'   \code{\link{postPEV}} after fitting when \code{coefficientMatrix=TRUE}.
#' @param execution Named list of memory, hardware, threading and numerical precision policies.
#' @param diagnostics Named list of ordinary verbosity, package-age messages and phase profiling.
#' @return A \code{sommer_mmes_control} object with strictly validated module names.
#' @section Compatibility:
#' Settings are checked against one compatibility table before fitting. A supplied
#' setting that conflicts with another, or that has no effect on the selected path,
#' stops with a message naming the offending module settings. Examples include PCG
#' estimator settings with a deterministic solver, Henderson-only settings with
#' \code{formulation="direct"}, and estimation modules with \code{task$mode="solve"}.
#' Rules that depend on the data, family, \code{vcc} or rotation are checked by
#' \code{\link{mmes}}.
#' @export
mmes_control <- function(formulation="henderson", task=list(),
                         representation=list(), solve=list(), logdet=list(),
                         score=list(), curvature=list(), update=list(),
                         proposal=list(), convergence=list(), output=list(),
                         execution=list(), diagnostics=list()){
  modules <- list(task=task, representation=representation, solve=solve,
    logdet=logdet, score=score, curvature=curvature, update=update,
    proposal=proposal, convergence=convergence, output=output,
    execution=execution, diagnostics=diagnostics)
  supplied <- c(if(!missing(formulation)) "formulation",
    unlist(lapply(names(modules), function(module){
      value <- modules[[module]]
      if(is.list(value) && length(names(value))) paste0(module, "$", names(value))
    }), use.names=FALSE))
  control <- do.call(.mmes_control_build, c(list(formulation=formulation), modules))
  .mmes_check_compatibility(control, supplied)
  attr(control, "supplied") <- supplied
  control
}

.mmes_control_build <- function(formulation="henderson", task=list(),
                         representation=list(), solve=list(), logdet=list(),
                         score=list(), curvature=list(), update=list(),
                         proposal=list(), convergence=list(), output=list(),
                         execution=list(), diagnostics=list()){
  base <- .mmes_engine_control()
  formulation <- .mme_choice(formulation, c("henderson", "direct"), "formulation")
  task <- .mme_module(task, list(mode="fit", objective="reml",
    pql=.mme_pql_defaults), "task")
  task$mode <- .mme_choice(task$mode, c("fit", "solve"), "task$mode")
  task$objective <- .mme_choice(task$objective, c("reml", "ml"), "task$objective")
  task$pql <- .mme_module(task$pql, .mme_pql_defaults, "task$pql")
  task$pql$maxit <- .mme_positive(task$pql$maxit, "task$pql$maxit", TRUE)
  task$pql$tol <- .mme_positive(task$pql$tol, "task$pql$tol")
  for(name in c("estimateTheta", "warmStart")){
    task$pql[[name]] <- .mme_flag(task$pql[[name]], paste0("task$pql$", name))
  }
  representation <- .mme_module(representation, list(
    compression="off", factorScores="none", fixedRank="error",
    equations=if(base$experimental$saddlePoint == "force") "saddle-point" else "normal",
    sparseFactors=base$experimental$sparseFactors,
    matrixFree=base$experimental$matrixFreePcg), "representation")
  representation$compression <- .mme_compression_control(representation$compression)
  representation$fixedRank <- .mme_choice(representation$fixedRank,
    c("error", "drop"), "representation$fixedRank")
  representation$factorScores <- .mme_choice(representation$factorScores,
    c("none", "fixed-shape", "profile"), "representation$factorScores")
  representation$equations <- .mme_choice(representation$equations,
    c("normal", "saddle-point", "conditional"), "representation$equations")
  pcgSolve <- .mme_pcg_solve_names
  solve <- .mme_module(solve, c(list(method="auto"), base$pcg[pcgSolve],
    base$experimental[c("blockSchur", "blockBanded", "factorSchurMultiFactor", "dFirstML")]), "solve")
  solve$method <- .mme_choice(solve$method, c("auto", "ldlt", "cholmod", "pcg"), "solve$method")
  logdet <- .mme_module(logdet, list(scaling=base$pcg$logDetScaling,
    lanczosSteps=base$pcg$lanczosSteps), "logdet")
  score <- .mme_module(score, list(traceProbes=base$pcg$traceProbes,
    estimator=base$pcg$traceEstimator, rank=base$pcg$traceRank,
    covarianceChecks=base$experimental$pdCertificates,
    traceBatch=base$experimental$traceBatch), "score")
  curvature <- .mme_module(curvature, list(method="exact-ai",
    reduced=base$experimental$reducedAI), "curvature")
  curvature$method <- .mme_choice(curvature$method,
    c("exact-ai", "hybrid-ai"), "curvature$method")
  update <- .mme_module(update, list(method="ai", recovery="em", blendWeight=NULL,
    commonScale=base$experimental$commonScale), "update")
  update$method <- .mme_choice(update$method, c("ai", "em", "conditional-fa"), "update$method")
  if(update$method == "em" && representation$equations != "conditional"){
    stop("update$method='em' requires representation$equations='conditional'; ordinary AI workflows use update$recovery for EM recovery.",call.=FALSE)
  }
  if(update$method == "conditional-fa" && representation$equations == "saddle-point"){
    stop("update$method='conditional-fa' cannot use saddle-point equations.",call.=FALSE)
  }
  update$recovery <- .mme_choice(update$recovery, c("em", "pxem"), "update$recovery")
  update["blendWeight"] <- list(.mme_blend_weight(update$blendWeight))
  proposal <- .mme_module(proposal, list(stepWeight=NULL, extrapolation="none"), "proposal")
  proposal$extrapolation <- .mme_choice(proposal$extrapolation,
    c("none", "aitken"), "proposal$extrapolation")
  if(curvature$method == "hybrid-ai" && proposal$extrapolation != "none"){
    stop("Hybrid AI and Aitken extrapolation are not a supported combination.", call.=FALSE)
  }
  convergence <- .mme_module(convergence, list(maxIters=30L,
    likelihoodTol=1e-4, parameterTol=1e-4, accuracy=base$pcg$accuracy,
    maxTraceProbes=NULL, maxLanczosSteps=NULL,
    logDetTol=base$pcg$logDetTol, scoreTol=base$pcg$scoreTol), "convergence")
  convergence$maxIters <- .mme_positive(convergence$maxIters, "convergence$maxIters", TRUE)
  for(name in c("likelihoodTol", "parameterTol")){
    convergence[[name]] <- .mme_positive(convergence[[name]], paste0("convergence$", name))
  }
  output <- .mme_module(output, list(inverse=0L, covariance="full",
    coefficientMatrix=TRUE, setup=FALSE,
    selectedInverse=base$experimental$cholmodFinalSelected,
    panelSize=base$experimental$panelCi2), "output")
  if(!is.numeric(output$inverse) || length(output$inverse) != 1L ||
     is.na(output$inverse) || !output$inverse %in% 0:2){
    stop("output$inverse must be 0, 1, or 2.", call.=FALSE)
  }
  output$inverse <- as.integer(output$inverse)
  output$covariance <- .mme_covariance_information(output$covariance)
  output$panelSize <- .mme_positive(output$panelSize, "output$panelSize", TRUE)
  for(name in c("coefficientMatrix", "setup")){
    output[[name]] <- .mme_flag(output[[name]], paste0("output$", name))
  }
  execution <- .mme_module(execution, list(denseMemoryMB=base$denseMemoryMB,
    memoryPolicy=base$memoryPolicy, gpu=base$gpu, gpuMinSize=base$gpuMinSize,
    blasThreads=NULL, inverseTol=1e-6, seed=base$pcg$seed), "execution")
  execution$inverseTol <- .mme_positive(execution$inverseTol, "execution$inverseTol")
  diagnostics <- .mme_module(diagnostics, list(verbose=TRUE, dateWarning=TRUE,
    profile=base$experimental$profile), "diagnostics")
  for(name in c("verbose", "dateWarning")){
    diagnostics[[name]] <- .mme_flag(diagnostics[[name]], paste0("diagnostics$", name))
  }
  modes <- list(representation=c("sparseFactors", "matrixFree"),
    solve=c("blockSchur", "blockBanded", "factorSchurMultiFactor", "dFirstML"),
    score=c("covarianceChecks", "traceBatch"), curvature="reduced",
    update="commonScale", output="selectedInverse", diagnostics="profile")
  control <- list(formulation=formulation, task=task, representation=representation,
    solve=solve, logdet=logdet, score=score, curvature=curvature, update=update,
    proposal=proposal, convergence=convergence, output=output, execution=execution,
    diagnostics=diagnostics)
  for(module in names(modes)) for(name in modes[[module]]){
    control[[module]][[name]] <- .mme_experimental_mode(control[[module]][[name]])
  }
  compiled <- .mmes_compile_control(control)
  control$solve[pcgSolve] <- compiled$pcg[pcgSolve]
  control$logdet$scaling <- compiled$pcg$logDetScaling
  control$score$estimator <- compiled$pcg$traceEstimator
  control$score$rank <- compiled$pcg$traceRank
  control$score$traceProbes <- compiled$pcg$traceProbes
  control$logdet$lanczosSteps <- compiled$pcg$lanczosSteps
  control$convergence$accuracy <- compiled$pcg$accuracy
  control$convergence$maxTraceProbes <- compiled$pcg$maxTraceProbes
  control$convergence$maxLanczosSteps <- compiled$pcg$maxLanczosSteps
  control$execution[c("gpu", "gpuMinSize", "denseMemoryMB", "memoryPolicy")] <-
    compiled[c("gpu", "gpuMinSize", "denseMemoryMB", "memoryPolicy")]
  control$execution$seed <- compiled$pcg$seed
  class(control) <- "sommer_mmes_control"
  attr(control, "native") <- compiled
  control
}

.mmes_compile_control <- function(control){
  mapping <- .mme_pcg_mapping()
  mapping <- mapping[names(mapping) != "compression"]
  pcg <- lapply(mapping, function(path) .mme_path_value(control, path))
  experimental <- c(control$solve[c("blockSchur", "blockBanded", "factorSchurMultiFactor", "dFirstML")],
    list(sparseFactors=control$representation$sparseFactors,
      matrixFreePcg=control$representation$matrixFree,
      pdCertificates=control$score$covarianceChecks, traceBatch=control$score$traceBatch,
      reducedAI=control$curvature$reduced, commonScale=control$update$commonScale,
      saddlePoint=if(control$representation$equations == "saddle-point") "force" else "off",
      cholmodFinalSelected=control$output$selectedInverse,
      panelCi2=control$output$panelSize, profile=control$diagnostics$profile))
  acceleration <- if(control$proposal$extrapolation == "aitken") "aitken"
    else if(control$curvature$method == "hybrid-ai") "hybrid-ai" else "none"
  engine <- .mmes_engine_control(solver=control$solve$method,
    experimental=experimental, pcg=pcg, acceleration=acceleration,
    factorScoreAugmentation=control$representation$factorScores,
    compression=control$representation$compression,
    covarianceInformation=control$output$covariance,
    retainCoefficientMatrix=control$output$coefficientMatrix,
    denseMemoryMB=control$execution$denseMemoryMB,
    memoryPolicy=control$execution$memoryPolicy, gpu=control$execution$gpu,
    gpuMinSize=control$execution$gpuMinSize)
  engine$updateMethod <- control$update$method
  engine$conditionalEquations <- control$representation$equations == "conditional"
  engine$conditionalInverseTol <- control$execution$inverseTol
  engine$fixedRank <- control$representation$fixedRank
  class(engine) <- "sommer_mmes_engine"
  engine
}

.mmes_control_args <- function(control){
  list(nIters=control$convergence$maxIters,
    tolParConvLL=control$convergence$likelihoodTol,
    tolParConvNorm=control$convergence$parameterTol,
    tolParInv=control$execution$inverseTol, returnParam=control$output$setup,
    dateWarning=control$diagnostics$dateWarning, verbose=control$diagnostics$verbose,
    stepWeight=control$proposal$stepWeight, blend=control$update$recovery,
    blendWeight=control$update$blendWeight,
    henderson=control$formulation == "henderson", computeCi=control$output$inverse,
    solver=control$solve$method, REML=control$task$objective == "reml",
    pqlControl=control$task$pql, solveOnly=control$task$mode == "solve",
    blasThreads=control$execution$blasThreads, engine=attr(control, "native"))
}

print.sommer_mmes_control <- function(x, ...){
  cat("-- sommer modular numerical controls --\n")
  cat("Formulation:", x$formulation, "\n")
  for(module in setdiff(names(x), "formulation")){
    values <- vapply(x[[module]], function(value){
      if(is.null(value)) "NULL" else if(is.list(value)) "<list>"
      else if(length(value) > 4L) paste0("<", length(value), " values>")
      else paste(value, collapse=",")
    }, character(1))
    cat(module, ": ", paste(paste0(names(values), "=", values), collapse=" | "), "\n", sep="")
  }
  invisible(x)
}

mmes <- function(fixed, random, rcov, data, W, weights=NULL,
                 naMethodX="exclude", naMethodY="exclude",
                 naMethodRandom="exclude", naMethodR="exclude",
                 contrasts=NULL, family=stats::gaussian(), vcc=NULL,
                 covPar=NULL, control=mmes_control(), .internal=NULL, ...){
  if(!inherits(control, "sommer_mmes_control")){
    stop("control must be created by mmes_control().", call.=FALSE)
  }
  legacy <- list(...)
  if(length(legacy)){
    if(is.null(names(legacy)) || anyNA(names(legacy)) ||
       any(!nzchar(names(legacy))) || anyDuplicated(names(legacy))){
      stop("Legacy mmes() arguments in ... must have unique, non-empty names.",
           call.=FALSE)
    }
    unknown <- setdiff(names(legacy), c("getPEV", "tolParInv", "verbose"))
    if(length(unknown)){
      stop("unused argument(s): ", paste(unknown, collapse=", "), call.=FALSE)
    }
  }
  supplied <- attr(control, "supplied")
  if(length(legacy)){
    if(is.null(supplied)) supplied <- c("formulation", .mme_control_paths(control))
    replacements <- c(getPEV="output$inverse", tolParInv="execution$inverseTol",
      verbose="diagnostics$verbose")
    if("getPEV" %in% names(legacy)){
      getPEV <- .mme_flag(legacy$getPEV, "getPEV")
      # getPEV only disabled PEV; TRUE did not request an inverse by itself.
      if(!getPEV) control$output$inverse <- 0L
    }
    if("tolParInv" %in% names(legacy)){
      control$execution["inverseTol"] <- legacy["tolParInv"]
    }
    if("verbose" %in% names(legacy)){
      control$diagnostics["verbose"] <- legacy["verbose"]
    }
    supplied <- unique(c(supplied, unname(replacements[names(legacy)])))
    control <- do.call(mmes_control, unclass(control))
    warning("Passing ", paste(names(legacy), collapse=", "),
      " through mmes(...) is deprecated. Start using control=mmes_control() with ",
      paste(paste0(names(legacy), " -> ", replacements[names(legacy)]),
        collapse=", "),
      ". Legacy arguments override the corresponding control settings.",
      if("getPEV" %in% names(legacy))
        " getPEV=FALSE sets output$inverse=0, while getPEV=TRUE leaves it unchanged.",
      call.=FALSE)
  }
  control <- do.call(.mmes_control_build, unclass(control))
  if(is.null(supplied)) supplied <- c("formulation", .mme_control_paths(control))
  .mmes_check_compatibility(control, supplied)
  attr(control, "supplied") <- supplied
  .mmes_check_family_controls(control, supplied, family, .internal)
  call <- match.call()
  model <- as.list(call)[-1L]
  model[c("control", ".internal", names(legacy))] <- NULL
  args <- lapply(model, eval, envir=parent.frame())
  args <- c(args, .mmes_control_args(control), list(.internal=.internal))
  result <- do.call(".mmes_fit", args, envir=environment(mmes))
  result$control <- control
  result$engineControl <- NULL
  if(is.list(result$engineDiagnostics)){
    result$engineDiagnostics$controlObject <- NULL
  }
  result$executionPlan <- .mmes_execution_plan(control, result)
  result$call <- call
  result
}

.mmes_check_family_controls <- function(control, supplied, family, internal){
  if(!is.null(internal) || !"task$pql" %in% supplied ||
     identical(control$task$pql, .mme_pql_defaults)) return(invisible(TRUE))
  if(inherits(family, "family") && identical(family$family, "gaussian") &&
     identical(family$link, "identity")){
    stop("Incompatible mmes() settings:\n* task$pql only applies to non-Gaussian families fitted by PQL. Offending: task$pql. Remove it or supply a non-Gaussian family.",
         call.=FALSE)
  }
  invisible(TRUE)
}

.mmes_execution_plan <- function(control, result){
  diagnostics <- result$engineDiagnostics
  solver <- if(control$task$mode == "solve") "pcg"
    else if(control$formulation == "direct") "direct"
    else if(!is.null(result$solver)) result$solver
    else diagnostics$control$solver
  if(is.null(solver)) solver <- NA_character_
  selected <- function(name){
    diagnostic <- diagnostics[[name]]
    if(is.null(diagnostic)){
      return(list(active=NA, reason="Not reported by this path; inspect model-specific diagnostics."))
    }
    list(active=if(is.null(diagnostic$active)) NA else diagnostic$active,
      reason=if(is.null(diagnostic$fallbackReason)) "" else diagnostic$fallbackReason)
  }
  list(version=1L, task=control$task$mode, objective=control$task$objective,
    formulation=control$formulation, solver=list(requested=control$solve$method,
      selected=solver),
    representation=control$representation,
    curvature=control$curvature, proposal=control$proposal,
    admission=list(compression=selected("compression"),
      reducedAI=selected("reducedAI"), commonScale=selected("commonScale"),
      hybridAI=selected("hybridAI"), saddlePoint=selected("saddlePoint"),
      conditionalFA=selected("conditionalFA")),
    update=control$update,
    output=control$output,
    status=if(control$output$setup && control$task$mode == "fit") "assembled" else "evaluated")
}
