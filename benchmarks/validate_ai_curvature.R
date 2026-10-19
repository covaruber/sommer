# Compare AI curvature and stationary-AR derivative implementations.
# Usage: Rscript benchmarks/validate_ai_curvature.R baseline candidate output.rds
args <- commandArgs(trailingOnly = TRUE)
if(length(args) && args[1L] == "--child"){
  library_path <- normalizePath(args[2L], mustWork = TRUE)
  output <- args[3L]
  suppressPackageStartupMessages(library(sommer, lib.loc = library_path))
  workload <- generate_mme_workload(structure = "ar2", n_levels = 1000L,
                                    n_ind = 6L, relationship = "none", seed = 391L)
  fit_model <- function(){
    mmes(workload$fixed,
  random=workload$random,
  rcov=workload$rcov,
  data=workload$data,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(solve=list(method="ldlt", blockSchur="off"),
  diagnostics=list(profile="force")
)), list(task=list(objective=if(TRUE) "reml" else "ml"), convergence=list(maxIters=4L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
  }
  invisible(fit_model())
  invisible(gc())
  elapsed <- system.time(fit <- fit_model())[["elapsed"]]
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")
  signature <- setNames(lapply(fields, function(field) {
    as.numeric(unlist(fit[[field]], use.names = FALSE))
  }), fields)
  saveRDS(list(signature = signature, elapsed = elapsed,
               profile = fit$engineDiagnostics$profile,
               convergence = fit$convergence,
               library = find.package("sommer")), output)
  quit(status = 0L)
}

if(length(args) < 3L){
  stop("usage: validate_ai_curvature.R baseline candidate output.rds", call. = FALSE)
}
if(!file.exists("/usr/bin/time")){
  stop("/usr/bin/time is required for peak RSS measurement", call. = FALSE)
}
libraries <- c(baseline = normalizePath(args[1L], mustWork = TRUE),
               candidate = normalizePath(args[2L], mustWork = TRUE))
output <- args[3L]
script <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)[1L]))
rscript <- file.path(R.home("bin"), "Rscript")
darwin <- identical(Sys.info()[["sysname"]], "Darwin")
runs <- paste0(output, ".runs")
dir.create(runs, recursive = TRUE, showWarnings = FALSE)
rows <- list()
errors <- list()
fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")
tolerances <- setNames(rep(1e-8, length(fields)), fields)
# Banded AR precision differentiation changes only roundoff-sensitive information and score.
tolerances[c("InfMat", "dLu")] <- 1e-6

for(repetition in 1:3){
  order <- if(repetition %% 2L) c("baseline", "candidate") else c("candidate", "baseline")
  pair <- list()
  for(mode in order){
    path <- file.path(runs, sprintf("rep-%d-%s.rds", repetition, mode))
    rss_path <- sub("\\.rds$", ".rss", path)
    log_path <- sub("\\.rds$", ".log", path)
    time_args <- if(darwin) {
      c("-l", "-o", shQuote(rss_path), shQuote(rscript), shQuote(script),
        "--child", shQuote(libraries[[mode]]), shQuote(path))
    } else {
      c("-f", "%M", "-o", shQuote(rss_path), shQuote(rscript), shQuote(script),
        "--child", shQuote(libraries[[mode]]), shQuote(path))
    }
    status <- system2("/usr/bin/time", time_args,
      env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1",
              "MKL_NUM_THREADS=1", "BLIS_NUM_THREADS=1",
              "VECLIB_MAXIMUM_THREADS=1"),
      stdout = log_path, stderr = log_path)
    if(status != 0L){ stop("Child failed; see ", log_path, call. = FALSE) }
    pair[[mode]] <- readRDS(path)
    profile <- pair[[mode]]$profile
    counters <- profile$counters
    rss_lines <- readLines(rss_path, warn = FALSE)
    if(darwin){
      rss_line <- grep("maximum resident set size", rss_lines, value = TRUE)
      if(length(rss_line) != 1L){ stop("Could not parse peak RSS; see ", rss_path, call. = FALSE) }
      peak_rss_kb <- as.numeric(sub("^\\s*([0-9]+).*", "\\1", rss_line)) / 1024
    } else {
      peak_rss_kb <- as.numeric(rss_lines[1L])
    }
    rows[[length(rows) + 1L]] <- data.frame(
      repetition = repetition, mode = mode,
      elapsed_seconds = pair[[mode]]$elapsed,
      ai_seconds = unname(profile$seconds[["aiAssembly"]]),
      random_curvature_seconds = counters$aiRandomCurvatureSeconds,
      trace_prep_seconds = unname(profile$seconds[["tracePrep"]]),
      peak_rss_kb = peak_rss_kb,
      converged = pair[[mode]]$convergence
    )
    cat(repetition, mode, pair[[mode]]$elapsed, "seconds; AI curvature",
        counters$aiRandomCurvatureSeconds, "seconds\n")
  }
  reference <- pair$baseline$signature
  candidate <- pair$candidate$signature
  difference <- vapply(fields, function(field){
    first <- reference[[field]]
    second <- candidate[[field]]
    if(length(first) != length(second)) return(Inf)
    if(!length(first)) return(0)
    if(any(!is.finite(first)) || any(!is.finite(second))) return(Inf)
    max(abs(second - first) / (1 + abs(first)))
  }, numeric(1))
  errors[[repetition]] <- difference
  if(any(difference > tolerances[names(difference)])){
    stop("Parity failed in repetition ", repetition, ": ",
         paste(names(difference)[difference > tolerances[names(difference)]],
               collapse = ", "),
         call. = FALSE)
  }
}

results <- do.call(rbind, rows)
summary <- aggregate(cbind(elapsed_seconds, ai_seconds, random_curvature_seconds,
                           peak_rss_kb) ~ mode, data = results, FUN = median)
paired <- merge(results[results$mode == "baseline", ],
                results[results$mode == "candidate", ],
                by = "repetition", suffixes = c("_baseline", "_candidate"))
paired$whole_fit_ratio <- paired$elapsed_seconds_baseline /
  paired$elapsed_seconds_candidate
paired$curvature_ratio <- paired$random_curvature_seconds_baseline /
  paired$random_curvature_seconds_candidate
timingUnstable <- any(vapply(
  paired[c("whole_fit_ratio", "curvature_ratio")],
  function(ratio){
    !all(is.finite(ratio)) || max(ratio) / min(ratio) > 1.2
  },
  logical(1)
))
saveRDS(list(results = results, paired = paired, summary = summary, errors = errors,
             libraries = libraries, session = sessionInfo(),
             note = "Four-iteration REML screening fit; timings are not converged-fit estimates."),
        output)
write.csv(results, sub("\\.rds$", ".csv", output), row.names = FALSE)
print(summary, row.names = FALSE)
print(paired[, c("repetition", "whole_fit_ratio", "curvature_ratio")], row.names = FALSE)
if(timingUnstable){
  cat("TIMING_UNSTABLE: paired ratios vary by more than 20%; do not report a speedup.\n")
} else {
  cat("Median paired whole-fit speedup:", median(paired$whole_fit_ratio), "\n")
}
cat("Maximum relative output difference:", max(unlist(errors)), "\n")
cat("VALIDATION_OK\n")
