# Run in a fresh process with the desired installed sommer library selected.
# Example: R_LIBS=/path/to/library OMP_NUM_THREADS=1 Rscript benchmarks/benchmark_hybrid_ai.R --quick
library(sommer)

quick <- "--quick" %in% commandArgs(trailingOnly=TRUE)
replications <- if(quick) 3L else 5L
levels <- if(quick) 6L else 12L
individuals <- 30L
output <- Sys.getenv("SOMMER_HYBRID_RESULTS", unset="")

fixture <- function(relationship){
  parameters <- if(relationship == "none") levels else 6L
  size <- if(!quick && relationship == "none") 1000L else individuals
  workload <- generate_mme_workload(structure="us", n_levels=parameters,
    n_ind=size, relationship=relationship, seed=391L)
  dataset <- workload$data
  covariance <- 0.5^abs(outer(seq_len(parameters), seq_len(parameters), "-"))
  effects <- matrix(rnorm(size * parameters), size) %*% chol(covariance)
  dataset$response <- as.integer(dataset$level) / 3 +
    effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
    rnorm(nrow(dataset), sd=0.7)
  workload$data <- dataset
  workload
}

rows <- list()
errors <- list()
for(relationship in c("none", "sparse", "dense")){
  workload <- fixture(relationship)
  for(reml in c(TRUE, FALSE)){
    fitModel <- function(mode){
      mmes(workload$fixed,
  random=workload$random,
  rcov=workload$rcov,
  data=workload$data,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(proposal=list(extrapolation=if(mode == "aitken") "aitken" else "none"),
  curvature=list(method=if(mode == "hybrid-ai") "hybrid-ai" else "exact-ai"),
  diagnostics=list(profile="force")
)), list(task=list(objective=if(reml) "reml" else "ml"), convergence=list(maxIters=150L, likelihoodTol=1e-10, parameterTol=1e-5), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
    }
    invisible(fitModel("none"))
    invisible(fitModel("hybrid-ai"))
    for(replication in seq_len(replications)){
      order <- if(replication %% 2L) c("none", "hybrid-ai") else c("hybrid-ai", "none")
      fits <- list()
      for(mode in order){
        started <- proc.time()[["elapsed"]]
        fit <- fitModel(mode)
        elapsed <- proc.time()[["elapsed"]] - started
        fits[[mode]] <- fit
        diagnostic <- fit$engineDiagnostics$hybridAI
        phases <- as.list(fit$engineDiagnostics$profile$seconds)
        names(phases) <- paste0("phase_", names(phases))
        rows[[length(rows) + 1L]] <- as.data.frame(c(list(
          relationship=relationship, individuals=workload$dimensions$n_ind,
          levels=workload$dimensions$n_levels,
          REML=reml, replication=replication,
          acceleration=mode, elapsed=elapsed, iterations=length(fit$llik),
          converged=fit$convergence, likelihood=tail(c(fit$llik), 1L),
          sensitivityColumns=diagnostic$sensitivityRHSColumns,
          exactRefreshes=diagnostic$exactRefreshes,
          bfgsUpdates=diagnostic$bfgsUpdates), phases))
      }
      relativeError <- function(field){
        first <- unlist(fits$none[[field]], use.names=FALSE)
        second <- unlist(fits[["hybrid-ai"]][[field]], use.names=FALSE)
        max(abs(first - second) / pmax(1, abs(first), abs(second)))
      }
      errors[[length(errors) + 1L]] <- data.frame(relationship=relationship,
        REML=reml, replication=replication,
        likelihood=abs(tail(c(fits$none$llik), 1L) - tail(c(fits[["hybrid-ai"]]$llik), 1L)),
        parameters=relativeError("covPar"), coefficients=relativeError("bu"),
        uncertainty=relativeError("theta_se"))
    }
    message("Completed ", relationship, if(reml) " REML" else " ML")
  }
}
result <- list(rows=do.call(rbind, rows), errors=do.call(rbind, errors),
  session=sessionInfo(), dimensions=c(levels=levels, individuals=individuals,
    identityIndividuals=if(quick) individuals else 1000L),
  threads=Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")))
print(aggregate(cbind(elapsed, iterations, sensitivityColumns) ~ relationship + REML + acceleration,
  data=result$rows, FUN=median), row.names=FALSE)
print(result$errors, row.names=FALSE)
if(nzchar(output)) saveRDS(result, output)
if(!all(result$rows$converged)) stop("Some benchmark fits did not converge.")
if(any(result$errors$likelihood > 1e-6 | result$errors$parameters > 1e-4 |
       result$errors$coefficients > 1e-4 | result$errors$uncertainty > 1e-4)){
  stop("Benchmark final outputs did not meet parity tolerances.")
}
