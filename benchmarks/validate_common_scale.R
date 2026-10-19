# Usage: R_LIBS=candidate-lib Rscript benchmarks/validate_common_scale.R output.rds
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1L){
  stop("usage: validate_common_scale.R output.rds", call.=FALSE)
}
suppressPackageStartupMessages(library(sommer))
set.seed(913)
dataset <- expand.grid(id=factor(seq_len(24L)), level=factor(seq_len(3L)),
                       replicate=seq_len(3L))
effects <- matrix(rnorm(72L), 24L) %*%
  chol(0.4^abs(outer(seq_len(3L), seq_len(3L), "-")))
dataset$y <- as.integer(dataset$level) / 3 +
  effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
  rnorm(nrow(dataset), sd=0.8)
measurements <- list()
for(reml in c(TRUE, FALSE)){
  model <- function(mode, information){
    mmes(y~level, random=~vsm(ar1m(level), ism(id), sigma2=0.2),
      rcov=~vsm(ism(units), sigma2=0.3), data=dataset,
      control=mmes_control(task=list(objective=if(reml) "reml" else "ml"),
        solve=list(method="ldlt"), update=list(commonScale=mode),
        convergence=list(maxIters=80L, likelihoodTol=1e-11, parameterTol=1e-9),
        output=list(inverse=2L, covariance=information),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE, profile="force")))
  }
  reference <- model("off", "full")
  for(mode in c("off", "force")) for(information in c("full", "none")){
    invisible(model(mode, information))
  }
  for(repetition in seq_len(5L)){
    modes <- if(repetition %% 2L) c("off", "force") else c("force", "off")
    for(mode in modes) for(information in c("full", "none")){
      invisible(gc())
      elapsed <- system.time(fit <- model(mode, information))[["elapsed"]]
      if(!isTRUE(fit$convergence) || !isTRUE(reference$convergence)){
        stop("Common-scale benchmark did not converge.", call.=FALSE)
      }
      for(field in c("covPar", "bu", "Ci", "uPevList")){
        if(!isTRUE(all.equal(fit[[field]], reference[[field]], tolerance=2e-5))){
          stop("Common-scale benchmark parity failed for ", field, call.=FALSE)
        }
      }
      if(abs(tail(c(fit$llik), 1L) - tail(c(reference$llik), 1L)) > 1e-8){
        stop("Common-scale benchmark likelihood parity failed.", call.=FALSE)
      }
      diagnostics <- fit$engineDiagnostics
      sensitivities <- if(mode == "force") diagnostics$commonScale$sensitivityColumns else
        diagnostics$profile$counters$sensitivityRHSColumns
      measurements[[length(measurements) + 1L]] <- data.frame(
        repetition=repetition, likelihood=if(reml) "REML" else "ML",
        mode=mode, information=information, elapsed=elapsed,
        sensitivityColumns=sensitivities,
        coefficientEvaluations=diagnostics$coefficientEvaluations,
        acceptedIterations=ncol(fit$llik),
        finalLogLik=tail(c(fit$llik), 1L))
    }
  }
}
measurements <- do.call(rbind, measurements)
saveRDS(list(measurements=measurements, session=sessionInfo()), args[1L])
write.csv(measurements, sub("\\.rds$", ".csv", args[1L]), row.names=FALSE)
print(aggregate(cbind(elapsed, sensitivityColumns, coefficientEvaluations,
  acceptedIterations) ~ likelihood + mode + information,
  data=measurements, FUN=median), row.names=FALSE)
cat("COMMON_SCALE_VALIDATION_OK\n")
