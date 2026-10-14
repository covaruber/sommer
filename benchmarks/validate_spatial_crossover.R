# Rscript benchmarks/validate_spatial_crossover.R library output.rds [--quick]
args <- commandArgs(trailingOnly = TRUE)
make_weights <- function(model, levels, density, seed){
  set.seed(seed)
  weights <- matrix(0, levels, levels)
  if(model == "car"){
    ring <- cbind(seq_len(levels), c(seq_len(levels - 1L) + 1L, 1L))
    weights[ring] <- runif(levels, 0.5, 1.5)
    weights[cbind(ring[, 2L], ring[, 1L])] <- weights[ring]
    possible <- which(upper.tri(weights), arr.ind = TRUE)
    target <- max(levels, round(density * nrow(possible)))
    selected <- setdiff(seq_len(nrow(possible)), match(
      paste(pmin(ring[, 1L], ring[, 2L]), pmax(ring[, 1L], ring[, 2L])),
      paste(possible[, 1L], possible[, 2L])))
    if(target > levels){
      chosen <- possible[sample(selected, target - levels), , drop = FALSE]
      weights[chosen] <- runif(nrow(chosen), 0.5, 1.5)
      weights[cbind(chosen[, 2L], chosen[, 1L])] <- weights[chosen]
    }
  }else{
    cycle <- cbind(seq_len(levels), c(seq_len(levels - 1L) + 1L, 1L))
    weights[cycle] <- runif(levels, 0.5, 1.5)
    possible <- which(row(weights) != col(weights) & weights == 0, arr.ind = TRUE)
    target <- max(levels, round(density * levels * (levels - 1L)))
    additional <- max(0L, target - levels)
    if(additional > 0L){
      chosen <- possible[sample(seq_len(nrow(possible)), additional), , drop = FALSE]
      weights[chosen] <- runif(nrow(chosen), 0.5, 1.5)
    }
  }
  labels <- paste0("L", seq_len(levels))
  dimnames(weights) <- list(labels, labels)
  weights
}

run_child <- function(library, output, model, levels, density, mode){
  suppressPackageStartupMessages(library(sommer, lib.loc = library))
  idCount <- if(levels >= 1000L) 2L else 6L
  weights <- make_weights(model, levels, density, 303L)
  dataset <- expand.grid(id = factor(seq_len(idCount)),
                         level = factor(rownames(weights), levels = rownames(weights)),
                         replicate = seq_len(2L))
  set.seed(1303L)
  dataset$response <- rnorm(nrow(dataset)) + as.integer(dataset$id) / 7 +
    as.integer(dataset$level) / levels
  rho <- if(model == "sar") 0.05 / max(rowSums(abs(weights))) else 0.05
  shape <- if(model == "sar") sar(dataset$level, weights, rho = rho) else
    car(dataset$level, weights, rho = rho)
  fit_model <- function(){
    mmes(response ~ 1, random = ~vsm(shape, ism(id)), rcov = ~units, data = dataset,
         nIters = 3L, nThreads = 4L, verbose = FALSE, dateWarning = FALSE,
         engine = mmes_control(solver = "ldlt", experimental = list(
           profile = "force", blockSchur = "off", sparseFactors = mode)))
  }
  invisible(fit_model())
  invisible(gc())
  elapsed <- system.time(fit <- fit_model())[["elapsed"]]
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")
  signature <- setNames(lapply(fields, function(field) as.numeric(unlist(fit[[field]]))), fields)
  saveRDS(list(signature = signature, elapsed = elapsed, profile = fit$engineDiagnostics$profile,
               cacheReuses = fit$engineDiagnostics$sparseSpatialFactorizationReuses), output)
}

if(length(args) && args[1L] == "--child"){
  run_child(args[2L], args[3L], args[4L], as.integer(args[5L]), as.numeric(args[6L]), args[7L])
  quit(status = 0L)
}
if(length(args) < 2L) stop("usage: validate_spatial_crossover.R library output.rds [--quick]")
if(!file.exists("/usr/bin/time")) stop("/usr/bin/time is required for peak RSS measurement")
library_path <- normalizePath(args[1L], mustWork = TRUE)
output <- args[2L]
quick <- "--quick" %in% args
large <- "--large" %in% args
repetitions <- if(quick) 1L else 3L
fixtures <- if(large){
  data.frame(model = c("sar", "car"), levels = 1000L, density = 1.0)
}else if(quick){
  data.frame(model = c("sar", "sar", "car", "car"), levels = 32L,
             density = c(0.1, 1, 0.1, 1))
}else{
  expand.grid(model = c("sar", "car"), levels = c(60L, 150L, 300L),
              density = c(0.05, 0.25, 1), KEEP.OUT.ATTRS = FALSE,
              stringsAsFactors = FALSE)
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
      path <- file.path(directory, sprintf("%s-%d-d%.2f-rep-%d-%s.rds",
        fixtures$model[fixture], fixtures$levels[fixture], fixtures$density[fixture], repetition, mode))
      rssPath <- sub("\\.rds$", ".rss", path)
      logPath <- sub("\\.rds$", ".log", path)
      status <- system2("/usr/bin/time", c("-f", "%M", "-o", shQuote(rssPath),
        shQuote(file.path(R.home("bin"), "Rscript")), shQuote(script), "--child",
        shQuote(library_path), shQuote(path), fixtures$model[fixture], fixtures$levels[fixture],
        fixtures$density[fixture], mode),
        env = c("OMP_NUM_THREADS=4", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1", "BLIS_NUM_THREADS=1"),
        stdout = logPath, stderr = logPath)
      if(status != 0L) stop("Child failed; see ", logPath)
      pair[[mode]] <- readRDS(path)
      rows[[length(rows) + 1L]] <- data.frame(
        model = fixtures$model[fixture], levels = fixtures$levels[fixture],
        requested_density = fixtures$density[fixture], repetition = repetition, mode = mode,
        elapsed_seconds = pair[[mode]]$elapsed,
        factor_eval_seconds = unname(pair[[mode]]$profile$seconds["factorEval"]),
        peak_rss_kb = as.numeric(readLines(rssPath, warn = FALSE)),
        factor_cache_reuses = pair[[mode]]$cacheReuses)
      cat(fixtures$model[fixture], fixtures$levels[fixture], fixtures$density[fixture],
          repetition, mode, pair[[mode]]$elapsed, "seconds\n")
    }
    difference <- vapply(names(pair$off$signature), function(field){
      reference <- pair$off$signature[[field]]
      candidate <- pair$force$signature[[field]]
      if(length(reference) != length(candidate)) return(Inf)
      if(!length(reference)) return(0)
      if(any(!is.finite(reference)) || any(!is.finite(candidate))) return(Inf)
      max(abs(candidate - reference) / (1 + abs(reference)))
    }, numeric(1))
    errors[[paste(fixtures$model[fixture], fixtures$levels[fixture], fixtures$density[fixture], repetition)]] <- difference
    if(any(difference > 1e-6)) stop("Parity failed: ", paste(names(difference)[difference > 1e-6], collapse = ", "))
  }
}
results <- do.call(rbind, rows)
summary <- aggregate(cbind(elapsed_seconds, factor_eval_seconds, peak_rss_kb, factor_cache_reuses) ~
  model + levels + requested_density + mode, data = results, FUN = median)
saveRDS(list(results = results, summary = summary, errors = errors, library = library_path,
             session = sessionInfo(),
             memory_note = "Peak process RSS includes package loading, warmup and timed fit."), output)
write.csv(results, sub("\\.rds$", ".csv", output), row.names = FALSE)
print(summary, row.names = FALSE)
cat("Maximum relative output difference:", max(unlist(errors)), "\n")
cat("VALIDATION_OK\n")