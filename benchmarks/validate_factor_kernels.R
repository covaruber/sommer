# Rscript benchmarks/validate_factor_kernels.R library output.rds [--quick]
args <- commandArgs(trailingOnly = TRUE)
if(length(args) && args[1L] == "--child"){
  suppressPackageStartupMessages(library(sommer, lib.loc = args[2L]))
  output <- args[3L]
  model <- args[4L]
  levels <- as.integer(args[5L])
  mode <- args[6L]
  individuals <- if(model == "ante" && levels > 16L) 8L else if(model == "ante") 32L else 6L
  workload <- generate_mme_workload(structure = if(model %in% c("ante", "ar3")) "ar2" else model,
                                   n_levels = levels, n_ind = individuals, seed = 303L)
  if(model == "ar3") workload$random <- ~vsm(ar3m(level), ism(id))
  if(model == "ante"){
    set.seed(303L)
    covariance <- 0.4 ^ abs(outer(seq_len(levels), seq_len(levels), "-"))
    effects <- matrix(rnorm(individuals * levels), individuals, levels) %*% chol(covariance)
    dataset <- workload$data
    workload$data$response <- effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
      rnorm(nrow(dataset), sd = 0.5) + as.integer(dataset$level) / levels
    workload$random <- ~vsm(antem(level, order = 2L), ism(id))
  }
  fit_model <- function(){
    mmes(workload$fixed,
  random=workload$random,
  rcov=workload$rcov,
  data=workload$data,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(solve=list(method="ldlt", blockSchur="off"),
  diagnostics=list(profile="force"),
  score=list(traceBatch="force"),
  representation=list(sparseFactors=mode)
)), list(convergence=list(maxIters=if(model == "ante" && levels > 16L) 2L else 4L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
  }
  invisible(fit_model())
  invisible(gc())
  elapsed <- system.time(fit <- fit_model())[["elapsed"]]
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")
  signature <- setNames(lapply(fields, function(field) as.numeric(unlist(fit[[field]]))), fields)
  saveRDS(list(signature = signature, elapsed = elapsed, profile = fit$engineDiagnostics$profile,
               library = find.package("sommer")), output)
  quit(status = 0L)
}
if(length(args) < 2L) stop("usage: validate_factor_kernels.R library output.rds [--quick]")
if(!file.exists("/usr/bin/time")) stop("/usr/bin/time is required for peak RSS measurement")
library_path <- normalizePath(args[1L], mustWork = TRUE)
output <- args[2L]
quick <- "--quick" %in% args
repetitions <- if(quick) 1L else 3L
fixtures <- if(quick){
  data.frame(model = c("ar2", "ar3", "ante", "sar", "car"), levels = c(30L, 30L, 8L, 30L, 30L))
}else{
  data.frame(model = c("ar2", "ar2", "ar3", "ar3", "ante", "ante", "ante", "ante",
                       "sar", "sar", "car", "car"),
             levels = c(200L, 1000L, 200L, 1000L, 8L, 16L, 64L, 128L,
                        200L, 1000L, 200L, 1000L))
}
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
directory <- paste0(output, ".runs")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
rows <- list()
errors <- list()
for(fixture in seq_len(nrow(fixtures))){
  for(repetition in seq_len(repetitions)){
    modes <- if(repetition %% 2L) c("off", "force") else c("force", "off")
    pair <- list()
    for(mode in modes){
      path <- file.path(directory, sprintf("%s-%d-rep-%d-%s.rds", fixtures$model[fixture],
                                           fixtures$levels[fixture], repetition, mode))
      rss_path <- sub("\\.rds$", ".rss", path)
      log_path <- sub("\\.rds$", ".log", path)
      status <- system2("/usr/bin/time", c("-f", "%M", "-o", shQuote(rss_path),
        shQuote(file.path(R.home("bin"), "Rscript")), shQuote(script), "--child",
        shQuote(library_path), shQuote(path), fixtures$model[fixture], fixtures$levels[fixture], mode),
        env = c("OMP_NUM_THREADS=4", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1", "BLIS_NUM_THREADS=1"),
        stdout = log_path, stderr = log_path)
      if(status != 0L) stop("Child failed; see ", log_path)
      pair[[mode]] <- readRDS(path)
      phases <- pair[[mode]]$profile$seconds
      rows[[length(rows) + 1L]] <- data.frame(model = fixtures$model[fixture], levels = fixtures$levels[fixture],
        repetition = repetition, mode = mode, elapsed_seconds = pair[[mode]]$elapsed,
        peak_rss_kb = as.numeric(readLines(rss_path, warn = FALSE)),
        factor_eval_seconds = unname(phases["factorEval"]), trace_prep_seconds = unname(phases["tracePrep"]),
        ai_seconds = unname(phases["aiAssembly"]), proposal_seconds = unname(phases["proposal"]))
      cat(fixtures$model[fixture], fixtures$levels[fixture], repetition, mode, pair[[mode]]$elapsed, "seconds\n")
    }
    reference <- pair$off$signature
    candidate <- pair$force$signature
    difference <- vapply(names(reference), function(field){
      if(length(reference[[field]]) != length(candidate[[field]])) return(Inf)
      if(!length(reference[[field]])) return(0)
      if(any(!is.finite(reference[[field]])) || any(!is.finite(candidate[[field]]))) return(Inf)
      max(abs(candidate[[field]] - reference[[field]]) / (1 + abs(reference[[field]])))
    }, numeric(1))
    errors[[paste(fixtures$model[fixture], fixtures$levels[fixture], repetition, sep = "-")]] <- difference
    if(any(difference > 1e-6)) stop("Parity failed: ", paste(names(difference)[difference > 1e-6], collapse = ", "))
  }
}
results <- do.call(rbind, rows)
summary <- aggregate(cbind(elapsed_seconds, peak_rss_kb, factor_eval_seconds, trace_prep_seconds,
                           ai_seconds, proposal_seconds) ~ model + levels + mode, data = results, FUN = median)
saveRDS(list(results = results, summary = summary, errors = errors, library = library_path,
             session = sessionInfo(), memory_note = "Peak RSS includes package loading, warmup, and one timed whole fit."), output)
write.csv(results, sub("\\.rds$", ".csv", output), row.names = FALSE)
print(summary, row.names = FALSE)
cat("Maximum relative output difference:", max(unlist(errors)), "\n")
cat("VALIDATION_OK\n")
