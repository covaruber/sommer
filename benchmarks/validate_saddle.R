# Run against isolated installations with identical thread settings.
# Baseline: R_LIBS=/baseline Rscript benchmarks/validate_saddle.R --baseline --out=/tmp/baseline.rds
# Candidate: R_LIBS=/candidate Rscript benchmarks/validate_saddle.R --reference=/tmp/baseline.rds --out=/tmp/candidate.rds
suppressPackageStartupMessages(library(sommer))
arguments <- commandArgs(trailingOnly=TRUE)
value <- function(prefix){
  found <- arguments[startsWith(arguments, prefix)]
  if(length(found)) substring(found[1L], nchar(prefix) + 1L) else ""
}
baselineOnly <- "--baseline" %in% arguments
output <- value("--out=")
referencePath <- value("--reference=")
replications <- 3L
fields <- c("llik", "bu", "C", "theta", "covPar", "covParWorking", "InfMat",
            "dLu", "theta_se", "Ci", "uList", "uPevList", "monitor",
            "normMonitor", "convergence")
fixtures <- c("iid", "sparse", "dense", "residual-ar1", "diagonal-weight", "sparse-weight")
defaultOutputs <- list()
rows <- list()
errors <- list()
relativeError <- function(first, second){
  first <- if(inherits(first, "Matrix")) as.numeric(first) else unlist(first, use.names=FALSE)
  second <- if(inherits(second, "Matrix")) as.numeric(second) else unlist(second, use.names=FALSE)
  if(length(first) != length(second)) stop("Compared outputs have different lengths.")
  if(!length(first)) return(0)
  max(abs(first - second) / pmax(1, abs(first), abs(second)))
}
for(fixture in fixtures){
  set.seed(892)
  dataset <- expand.grid(id=factor(seq_len(30L)), level=factor(seq_len(6L)),
                         replicate=seq_len(2L))
  effects <- rnorm(30L)
  dataset$response <- as.integer(dataset$level) / 4 +
    effects[as.integer(dataset$id)] + rnorm(nrow(dataset), sd=0.6)
  random <- ~id
  if(fixture %in% c("sparse", "dense")){
    precision <- if(fixture == "dense") diag(30L) + matrix(0.02, 30L, 30L) else
      Matrix::bandSparse(30L, k=c(-1L, 0L, 1L),
        diagonals=list(rep(-0.1, 29L), rep(1, 30L), rep(-0.1, 29L)))
    dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
    attr(precision, "inverse") <- TRUE
    random <- ~vsm(ism(id), Gu=precision)
  }
  rcov <- if(fixture == "residual-ar1") ~vsm(ar1m(level, rho=0.2), ism(units)) else ~units
  weight <- NULL
  if(fixture == "diagonal-weight") weight <- diag(seq(0.7, 1.3, length.out=nrow(dataset)))
  if(fixture == "sparse-weight"){
    weight <- Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
      diagonals=list(rep(0.05, nrow(dataset)-1L), rep(1.2, nrow(dataset)),
                     rep(0.05, nrow(dataset)-1L)))
  }
  for(solver in c("auto", "ldlt", "cholmod")) for(reml in c(TRUE, FALSE)) for(mode in c(0L, 1L, 2L)){
    key <- paste(fixture, solver, reml, mode, sep="/")
    fitModel <- function(saddle){
      control <- if(baselineOnly) mmes_control(solver=solver, experimental=list(profile="force")) else
        mmes_control(solver=solver, experimental=list(profile="force", saddlePoint=if(saddle) "force" else "off"))
      args <- list(fixed=response~level, random=random, rcov=rcov, data=dataset,
        REML=reml, computeCi=mode, nIters=80L, tolParConvLL=1e-10, tolParConvNorm=1e-8,
        engine=control, verbose=FALSE, dateWarning=FALSE)
      if(!is.null(weight)) args$W <- weight
      do.call(mmes, args)
    }
    plain <- fitModel(FALSE)
    defaultOutputs[[key]] <- plain[fields]
    if(baselineOnly) next
    candidate <- fitModel(TRUE)
    if(!plain$convergence || !candidate$convergence) stop("Fit did not converge: ", key)
    if(mode != 0L){
      errors[[key]] <- vapply(fields, function(field){
        if(mode == 1L && field == "Ci"){
          return(relativeError(diag(plain$Ci), diag(candidate$Ci)))
        }
        relativeError(plain[[field]], candidate[[field]])
      }, numeric(1))
      if(mode == 1L){
        inverse <- solve(as.matrix(candidate$C)) * candidate$Cscale
        stored <- summary(candidate$Ci)
        error <- max(abs(stored$x - inverse[cbind(stored$i, stored$j)]) /
                       pmax(1, abs(stored$x)))
        if(error > 1e-7) stop("Selected inverse entries are inaccurate: ", key)
      }
      next
    }
    for(replication in seq_len(replications)){
      fits <- list()
      for(saddle in if(replication %% 2L) c(FALSE, TRUE) else c(TRUE, FALSE)){
        started <- proc.time()[["elapsed"]]
        fit <- fitModel(saddle)
        elapsed <- proc.time()[["elapsed"]] - started
        fits[[if(saddle) "saddle" else "default"]] <- fit
        phase <- fit$engineDiagnostics$profile
        rows[[length(rows)+1L]] <- data.frame(fixture=fixture, solver=solver, REML=reml,
          replication=replication, saddle=saddle, elapsed=elapsed, converged=fit$convergence,
          iterations=length(fit$llik), assembly=phase$seconds[["assembly"]],
          factorization=phase$seconds[["factorization"]],
          traces=phase$seconds[["scoreTraces"]], final=phase$seconds[["final"]])
      }
      errors[[paste(key, replication, sep="/")]] <- vapply(fields, function(field){
        relativeError(fits$default[[field]], fits$saddle[[field]])
      }, numeric(1))
    }
    message("Completed ", key)
  }
}
result <- list(defaultOutputs=defaultOutputs,
  rows=if(length(rows)) do.call(rbind, rows) else NULL,
  errors=if(length(errors)) do.call(rbind, errors) else NULL, session=sessionInfo(),
  threads=Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")))
if(nzchar(referencePath)){
  reference <- readRDS(referencePath)$defaultOutputs
  result$defaultIdentical <- vapply(names(reference), function(key){
    identical(reference[[key]], defaultOutputs[[key]])
  }, logical(1))
  if(!all(result$defaultIdentical)) stop("Default output changed against preserved baseline.")
}
if(nzchar(output)) saveRDS(result, output)
if(!baselineOnly){
  print(aggregate(cbind(elapsed, iterations, assembly, factorization, traces, final) ~ fixture + solver + REML + saddle,
    data=result$rows, FUN=median), row.names=FALSE)
  print(apply(result$errors, 2L, max))
  if(!all(result$rows$converged) || max(result$errors) > 1e-6) stop("Saddle validation gate failed.")
}
