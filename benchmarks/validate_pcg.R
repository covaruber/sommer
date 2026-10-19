arguments <- commandArgs(TRUE)
if(!length(arguments) %in% c(2L, 3L)){
  stop("Usage: Rscript benchmarks/validate_pcg.R LIBRARY OUTPUT.rds [VARIANT,...]")
}
if(!grepl("\\.rds$", arguments[2L])){
  stop("The validation output path must end in .rds.")
}
.libPaths(c(normalizePath(arguments[1L], mustWork=TRUE), .libPaths()))
suppressPackageStartupMessages(library(sommer))

variants <- list(
  baseline=list(),
  block_jacobi=list(solve=list(preconditioner="block-jacobi")),
  nystrom=list(solve=list(preconditioner="nystrom", nystromRank=12L)),
  warm=list(solve=list(warmStart=TRUE)),
  recycle=list(solve=list(warmStart=TRUE, recycleRank=6L)),
  block_cg=list(solve=list(rhsMethod="block", blockSize=4L)),
  scaled_slq=list(logdet=list(scaling="diagonal")),
  control_variate=list(score=list(estimator="control-variate", rank=12L))
)
if(length(arguments) == 3L){
  selected <- strsplit(arguments[3L], ",", fixed=TRUE)[[1L]]
  if(!length(selected) || any(!selected %in% names(variants))){
    stop("Unknown benchmark variant. Choose from: ", paste(names(variants), collapse=", "))
  }
  variants <- variants[selected]
}
rows <- list()
fits <- list()
for(relationship in c("none", "sparse", "dense")){
  workload <- generate_mme_workload(n_ind=40L, n_levels=3L,
    structure="us", relationship=relationship, seed=721L)
  for(repetition in seq_len(3L)){
    order <- if(repetition %% 2L) names(variants) else rev(names(variants))
    for(name in order){
      control <- do.call(mmes_control, utils::modifyList(list(
        solve=list(method="pcg", tol=1e-10), score=list(traceProbes=32L),
        logdet=list(lanczosSteps=40L), update=list(blendWeight=0),
        convergence=list(maxIters=40L,
          likelihoodTol=1e-5, parameterTol=1e-4),
        diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)),
        variants[[name]]))
      elapsed <- system.time(fit <- suppressMessages(mmes(
        workload$fixed, random=workload$random, rcov=workload$rcov,
        data=workload$data, control=control)))[["elapsed"]]
      inputs <- suppressMessages(mmes(workload$fixed,
        random=workload$random, rcov=workload$rcov, data=workload$data,
        control=mmes_control(output=list(setup=TRUE),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
      descriptors <- sommer:::.mmes_known_covpar(inputs$covStruct,
        names(fit$covStruct), fit)
      for(term in seq_along(descriptors)){
        descriptors[[term]]$free <- inputs$covStruct[[term]]$free
      }
      reference <- sommer:::ai_mme_sp2(
        inputs$X, inputs$Z, inputs$Zind, inputs$Ai, inputs$yvar, inputs$W,
        inputs$useH, inputs$residualBlock, inputs$residualIndex,
        1L, 1e-4, 1e-4, inputs$tolParInv, descriptors, 0, 1, FALSE, 0L,
        "ldlt", 1e-10, 0L, 8L, 20L, TRUE, inputs$responsePrepared,
        inputs$preparedMean, inputs$preparedSd, inputs$preparedIntercept)
      if(length(reference$bu) != length(fit$bu) ||
         !length(reference$dLu) || length(reference$dLu) != length(fit$dLu)){
        stop("Direct-reference output dimensions differ: ", relationship, "/", name)
      }
      coefficientError <- max(abs(as.numeric(reference$bu) -
        as.numeric(fit$bu))) / max(1, max(abs(reference$bu)))
      scoreError <- max(abs(as.numeric(reference$dLu) -
        as.numeric(fit$dLu))) / max(1, max(abs(reference$dLu)))
      if(!is.finite(scoreError)){
        stop("Direct-reference score comparison is non-finite: ", relationship, "/", name)
      }
      if(!is.finite(coefficientError) || coefficientError > 1e-6){
        stop("PCG coefficient reference failed: ", relationship, "/", name,
             " relative error ", coefficientError)
      }
      key <- paste(relationship, repetition, name, sep="/")
      rows[[key]] <- data.frame(
        relationship=relationship, repetition=repetition, variant=name,
        elapsed=elapsed, converged=isTRUE(fit$convergence),
        relativeCoefficientError=coefficientError,
        logLik=tail(as.numeric(fit$llik), 1L),
        directLogLik=tail(as.numeric(reference$llik), 1L),
        absoluteLogLikError=abs(tail(as.numeric(reference$llik), 1L) -
          tail(as.numeric(fit$llik), 1L)),
        relativeScoreError=scoreError,
        logDetProbeSE=fit$pcgLogDetProbeStandardError,
        solveResidual=fit$pcgSolveMaxRelativeResidual
      )
      fits[[key]] <- list(diagnostics=fit$engineDiagnostics,
        parameters=fit$covPar, controls=control)
      cat(key, "seconds", elapsed, "converged", fit$convergence, "\n")
    }
  }
}
result <- list(rows=do.call(rbind, rows), fits=fits,
  session=sessionInfo(), package=packageDescription("sommer"),
  nativeHashes=tools::md5sum(list.files(system.file("libs", package="sommer"),
    pattern="\\.(so|dll|dylib)$", full.names=TRUE, recursive=TRUE)),
  driverSourceHashes=tools::md5sum(c("src/ai_mme_sp2.cpp", "src/mme_pcg.h",
    "src/mme_pcg_estimators.h")))
saveRDS(result, arguments[2L])
write.csv(result$rows, sub("\\.rds$", ".csv", arguments[2L]), row.names=FALSE)
cat("Saved matched-point validation. Compare speed only for converged fits with comparable direct likelihood and score errors.\n")
