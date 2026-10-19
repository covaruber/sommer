#!/usr/bin/env Rscript

# Usage: Rscript benchmarks/validate_pcg_reuse.R LIBRARY OUTPUT.rds [REPEATS] [SOLVER]
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2L) stop("Expected LIBRARY OUTPUT.rds [REPEATS].")
.libPaths(c(normalizePath(args[[1L]], mustWork=TRUE), .libPaths()))
suppressPackageStartupMessages(library(sommer))
repeats <- if(length(args) >= 3L) as.integer(args[[3L]]) else 3L
if(length(repeats) != 1L || is.na(repeats) || repeats < 1L) stop("Invalid repeat count.")
solver <- if(length(args) >= 4L) match.arg(args[[4L]], c("pcg", "ldlt", "cholmod", "auto")) else "pcg"

set.seed(617)
dataset <- expand.grid(id=factor(seq_len(40L)), environment=factor(seq_len(3L)),
  replicate=seq_len(300L))
signal <- rnorm(nlevels(dataset$id), sd=1.2)
dataset$response <- signal[as.integer(dataset$id)] +
  as.integer(dataset$environment) / 3 + rnorm(nrow(dataset), sd=0.8)
variants <- list(
  baseline=list(),
  groups=list(representation=list(compression="groups")),
  classes=list(representation=list(compression="classes")),
  classes_krylov=list(representation=list(compression="classes"),
    solve=list(recycleMethod="krylov", recycleRank=8L)))
if(solver != "pcg") variants$classes_krylov <- NULL
run <- function(modules) {
  suppressMessages(mmes(response~environment, random=~id,
    rcov=~vsm(dsm(environment), ism(units)), data=dataset,
    control=do.call(mmes_control, utils::modifyList(list(
      solve=list(method=solver, tol=1e-10),
      score=list(traceProbes=16L, estimator="control-variate", rank=1000L),
      logdet=list(lanczosSteps=60L), convergence=list(maxIters=30L),
      update=list(blendWeight=0),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)), modules))))
}
reference <- run(variants$baseline)
records <- list()
fits <- list()
for(repetition in seq_len(repeats)){
  # Rotate execution order to limit systematic warm-cache bias.
  order <- ((seq_along(variants) + repetition - 2L) %% length(variants)) + 1L
  for(index in order){
    name <- names(variants)[index]
    gc()
    elapsed <- system.time(fit <- run(variants[[index]]))[["elapsed"]]
    fields <- c("bu", "theta", "llik", "dLu", "InfMat", "theta_se")
    errors <- vapply(fields, function(field) {
      candidate <- as.numeric(unlist(fit[[field]]))
      baseline <- as.numeric(unlist(reference[[field]]))
      if(length(candidate) != length(baseline)) return(Inf)
      if(!length(candidate)) return(0)
      max(abs(candidate - baseline)) / max(1, max(abs(baseline)))
    }, numeric(1L))
    diagnostics <- fit$engineDiagnostics$pcg
    compression <- fit$engineDiagnostics$compression
    records[[length(records) + 1L]] <- data.frame(
      variant=name, repetition=repetition, elapsed=elapsed,
      convergence=isTRUE(fit$convergence), maxRelativeFieldError=max(errors),
      trueRelativeResidual=if(solver == "pcg") fit$pcgSolveMaxRelativeResidual else NA_real_,
      operatorVectors=if(solver == "pcg") diagnostics$operatorVectorEquivalent else NA_real_,
      iterations=fit$engineDiagnostics$pcgBatchIterations,
      groups=compression$groups,
      retainedStatisticBytes=compression$retainedBytes,
      statisticSetupSeconds=compression$setupSeconds,
      recycleSetupSeconds=if(solver == "pcg") diagnostics$seconds$preconditionerAndCoarseSetup else 0)
    fits[[paste(name, repetition, sep="_")]] <- list(
      convergence=fit$convergence, errors=errors, diagnostics=diagnostics)
    print(records[[length(records)]])
  }
}
measurements <- do.call(rbind, records)
metadata <- list(session=sessionInfo(), package=packageDescription("sommer"),
  nativeHash=tools::md5sum(system.file("libs", .Platform$r_arch, "sommer.so", package="sommer")),
  observations=nrow(dataset), coefficientCount=length(reference$bu), seed=617L, solver=solver,
  description="Whole Gaussian REML fits, including grouping and model setup; matched full-rank trace control variates and SLQ budgets.")
saveRDS(list(measurements=measurements, fits=fits, metadata=metadata), args[[2L]])
write.csv(measurements, sub("\\.rds$", ".csv", args[[2L]]), row.names=FALSE)
if(any(!is.finite(measurements$maxRelativeFieldError)) ||
   any(measurements$maxRelativeFieldError > 1e-6) ||
   any(measurements$trueRelativeResidual > 1e-10, na.rm=TRUE)){
  stop("Accuracy-matched benchmark failed; inspect saved measurements.")
}
print(aggregate(elapsed ~ variant, measurements, median))
