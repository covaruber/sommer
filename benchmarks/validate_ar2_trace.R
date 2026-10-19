# Run from the package root with an explicitly installed candidate library.
# Rscript benchmarks/validate_ar2_trace.R library output.rds [--quick]
args <- commandArgs(trailingOnly = TRUE)
if(length(args) && args[1L] == "--child"){
  suppressPackageStartupMessages(library(sommer, lib.loc = args[2L]))
  output <- args[3L]
  threads <- as.integer(args[4L])
  mode <- args[5L]
  levels <- as.integer(args[6L])
  workload <- generate_mme_workload(structure = "ar2", n_levels = levels, n_ind = 6L, seed = 42L)
  fit_model <- function(){
    mmes(workload$fixed,
  random=workload$random,
  rcov=workload$rcov,
  data=workload$data,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(solve=list(method="ldlt", blockSchur="off"),
  diagnostics=list(profile="force"),
  score=list(traceBatch=mode)
)), list(convergence=list(maxIters=4L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
  }
  invisible(fit_model())
  invisible(gc())
  elapsed <- system.time(fit <- fit_model())[["elapsed"]]
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")
  signature <- setNames(lapply(fields, function(field) as.numeric(unlist(fit[[field]]))), fields)
  saveRDS(list(signature = signature, elapsed = elapsed,
               profile = fit$engineDiagnostics$profile, library = find.package("sommer")), output)
  quit(status = 0L)
}
if(length(args) < 2L) stop("usage: validate_ar2_trace.R library output.rds [--quick]")
if(!file.exists("/usr/bin/time")) stop("/usr/bin/time is required for peak RSS measurement")
library_path <- normalizePath(args[1L], mustWork = TRUE)
output <- args[2L]
quick <- "--quick" %in% args
levels <- if(quick) 40L else 1000L
repetitions <- if(quick) 1L else 3L
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
directory <- paste0(output, ".runs")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
rows <- list()
errors <- list()
for(threads in c(1L, 4L)){
  for(repetition in seq_len(repetitions)){
    modes <- if(repetition %% 2L) c("off", "force") else c("force", "off")
    pair <- list()
    for(mode in modes){
      path <- file.path(directory, sprintf("threads-%d-rep-%d-%s.rds", threads, repetition, mode))
      rss_path <- sub("\\.rds$", ".rss", path)
      log_path <- sub("\\.rds$", ".log", path)
      status <- system2("/usr/bin/time", c("-f", "%M", "-o", shQuote(rss_path),
        shQuote(file.path(R.home("bin"), "Rscript")), shQuote(script), "--child",
        shQuote(library_path), shQuote(path), threads, mode, levels),
        env = c(paste0("OMP_NUM_THREADS=", threads), "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
                "BLIS_NUM_THREADS=1"), stdout = log_path, stderr = log_path)
      if(status != 0L) stop("Child failed; see ", log_path)
      pair[[mode]] <- readRDS(path)
      counters <- pair[[mode]]$profile$counters
      rows[[length(rows) + 1L]] <- data.frame(
        levels = levels, threads = threads, repetition = repetition, mode = mode,
        elapsed_seconds = pair[[mode]]$elapsed,
        score_seconds = unname(pair[[mode]]$profile$seconds["scoreTraces"]),
        peak_rss_kb = as.numeric(readLines(rss_path, warn = FALSE)),
        fallback_columns = counters$traceFallbackSolves, batch_solves = counters$traceBatchSolves,
        panel_workspace_bytes = counters$traceBatchWorkspaceBytes)
      cat(threads, repetition, mode, pair[[mode]]$elapsed, "seconds\n")
    }
    reference <- pair$off$signature
    candidate <- pair$force$signature
    difference <- vapply(names(reference), function(field){
      if(length(reference[[field]]) != length(candidate[[field]])) return(Inf)
      if(!length(reference[[field]])) return(0)
      if(any(!is.finite(reference[[field]])) || any(!is.finite(candidate[[field]]))) return(Inf)
      max(abs(candidate[[field]] - reference[[field]]) / (1 + abs(reference[[field]])))
    }, numeric(1))
    if(any(difference > 1e-8)) stop("Parity failed: ", paste(names(difference)[difference > 1e-8], collapse = ", "))
    errors[[paste(threads, repetition, sep = "-")]] <- difference
  }
}
results <- do.call(rbind, rows)
summary <- aggregate(cbind(elapsed_seconds, score_seconds, peak_rss_kb) ~ threads + mode,
                     data = results, FUN = median)
saveRDS(list(results = results, summary = summary, errors = errors,
             library = library_path, session = sessionInfo(),
             memory_note = "Peak process RSS includes package loading, one warmup, and one timed whole fit."), output)
write.csv(results, sub("\\.rds$", ".csv", output), row.names = FALSE)
print(summary, row.names = FALSE)
cat("Maximum relative output difference:", max(unlist(errors)), "\n")
cat("VALIDATION_OK\n")
