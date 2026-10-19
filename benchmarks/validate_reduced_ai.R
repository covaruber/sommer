# Usage: Rscript benchmarks/validate_reduced_ai.R baseline-lib candidate-lib output.rds
args <- commandArgs(trailingOnly=TRUE)
if(length(args) && args[1L] == "--child"){
  suppressPackageStartupMessages(library(sommer, lib.loc=args[2L]))
  set.seed(331)
  dataset <- expand.grid(id=factor(seq_len(36L)), level=factor(seq_len(3L)),
                         replicate=seq_len(3L))
  dataset$h <- factor(rep(seq_len(18L), length.out=nrow(dataset)))
  dataset$y <- rnorm(nrow(dataset)) + rep(rnorm(36L), 9L) +
    as.integer(dataset$level) / 3
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "dLu", "InfMat",
              "theta_se", "Ci", "uPevList", "convergence")
  signatures <- list()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(tied in c(FALSE, TRUE)) for(fixed in c(FALSE, TRUE)){
      random <- if(fixed) ~vsm(ar1m(level, rho=0.3, fixed=TRUE), ism(id))+h
        else ~id+h
      fit <- mmes(y~level, random=random, data=dataset,
        vcc=if(tied) data.frame(parameter=if(fixed) c(1L, 3L) else 1:2,
          group=1L, scale=c(2, 1)) else NULL,
        control=mmes_control(convergence=list(maxIters=5L),
          task=list(objective=if(reml) "reml" else "ml"), solve=list(method=solver),
          output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      key <- paste(solver, reml, tied, fixed, sep="/")
      signatures[[key]] <- fit[fields]
    }
  }
  measurements <- list()
  if(identical(args[4L], "candidate")){
    terms <- 24L
    for(term in seq_len(terms)) dataset[[paste0("id", term)]] <- dataset$id
    random <- reformulate(paste0("id", seq_len(terms)))
    model <- function(mode, information){
      mmes(y~level, random=random, data=dataset,
        vcc=data.frame(parameter=seq_len(terms), group=1L),
        control=mmes_control(solve=list(method="ldlt"),
          curvature=list(reduced=mode), output=list(covariance=information),
          convergence=list(maxIters=15L, likelihoodTol=1e-10, parameterTol=1e-8),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE, profile="force")))
    }
    reference <- model("off", "full")
    for(mode in c("off", "force")) for(information in c("full", "none")){
      invisible(model(mode, information))
    }
    for(repetition in seq_len(5L)){
      order <- if(repetition %% 2L) c("off", "force") else c("force", "off")
      for(mode in order) for(information in c("full", "none")){
        invisible(gc())
        elapsed <- system.time(fit <- model(mode, information))[["elapsed"]]
        for(field in c("theta", "bu", "llik", "convergence")){
          if(!isTRUE(all.equal(fit[[field]], reference[[field]], tolerance=1e-7))){
            stop("Reduced-AI benchmark parity failed for ", field, call.=FALSE)
          }
        }
        if(!isTRUE(fit$convergence)) stop("Benchmark did not converge.", call.=FALSE)
        profile <- fit$engineDiagnostics$profile
        counters <- profile$counters
        measurements[[length(measurements) + 1L]] <- data.frame(
          repetition=repetition, mode=mode, information=information,
          elapsed=elapsed,
          sensitivities=unname(profile$seconds[["sensitivities"]]),
          aiAssembly=unname(profile$seconds[["aiAssembly"]]),
          sensitivityColumns=counters$sensitivityRHSColumns,
          fullParameters=terms + 1L, independentParameters=2L,
          acceptedIterations=ncol(fit$llik))
      }
    }
  }
  saveRDS(list(signatures=signatures,
    measurements=if(length(measurements)) do.call(rbind, measurements) else NULL,
    session=sessionInfo()), args[3L])
  quit(status=0L)
}

if(length(args) != 3L){
  stop("usage: validate_reduced_ai.R baseline-lib candidate-lib output.rds", call.=FALSE)
}
script <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value=TRUE)[1L]))
output <- args[3L]
libraries <- c(baseline=normalizePath(args[1L], mustWork=TRUE),
               candidate=normalizePath(args[2L], mustWork=TRUE))
results <- list()
for(mode in names(libraries)){
  child <- paste0(output, ".", mode, ".rds")
  log <- paste0(child, ".log")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), "--child", shQuote(libraries[[mode]]), shQuote(child), mode),
    env=c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
          "BLIS_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1"),
    stdout=log, stderr=log)
  if(status != 0L) stop("Validation child failed; see ", log, call.=FALSE)
  results[[mode]] <- readRDS(child)
}
reference <- results$baseline$signatures
candidate <- results$candidate$signatures
if(!identical(reference, candidate)){
  stop("Default-off outputs changed against the preserved baseline: ",
    paste(all.equal(reference, candidate, tolerance=0), collapse="; "), call.=FALSE)
}
measurements <- results$candidate$measurements
saveRDS(list(defaultParity=TRUE, libraries=libraries, results=results), output)
write.csv(measurements, sub("\\.rds$", ".csv", output), row.names=FALSE)
print(aggregate(cbind(elapsed, sensitivities, aiAssembly, sensitivityColumns) ~
  mode + information, data=measurements, FUN=median), row.names=FALSE)
cat("DEFAULT_PARITY_OK:", length(reference), "fits,",
    length(reference[[1L]]), "fields, bit-identical\n")
cat("REDUCED_AI_VALIDATION_OK\n")
