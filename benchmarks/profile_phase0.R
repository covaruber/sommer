# Phase 0 profile: where does fit time go for large AR(p), SAR/CAR and dense-Gu workloads.
# Usage: R_LIBS=/tmp/somlib Rscript benchmarks/profile_phase0.R [--quick] [--out=/tmp/phase0_profile.csv]
suppressMessages(library(sommer))

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
outArg <- grep("^--out=", args, value = TRUE)
output <- if(length(outArg)) sub("^--out=", "", outArg[1L]) else "/tmp/phase0_profile.csv"

fixtures <- rbind(
  expand.grid(structure = c("ar1", "ar2"), n_levels = if(quick) c(50L, 200L) else c(50L, 200L, 1000L),
              n_ind = 6L, relationship = "none", stringsAsFactors = FALSE),
  expand.grid(structure = c("sar", "car"), n_levels = if(quick) c(200L, 1000L) else c(200L, 1000L, 2000L),
              n_ind = 4L, relationship = "none", stringsAsFactors = FALSE),
  expand.grid(structure = c("ar1", "us"), n_levels = 5L,
              n_ind = if(quick) c(100L, 300L) else c(100L, 300L, 600L),
              relationship = "dense", stringsAsFactors = FALSE)
)
solvers <- c("ldlt", "cholmod")
nIters <- 4L

rows <- list()
for(i in seq_len(nrow(fixtures))){
  fx <- fixtures[i, ]
  wl <- generate_mme_workload(structure = fx$structure, n_levels = fx$n_levels,
                              n_ind = fx$n_ind, relationship = fx$relationship, reml = TRUE)
  for(solver in solvers){
    cat(sprintf("%-32s %-8s ... ", wl$name, solver))
    p <- tryCatch(
      profile_mme_workload(wl, mmes_control(solve=list(method=solver)
), nIters = nIters),
      error = function(e) conditionMessage(e))
    if(is.character(p)){ cat("ERROR:", p, "\n"); next }
    secs <- stats::setNames(p$table$seconds, p$table$phase)
    top <- p$table$phase[which.max(p$table$seconds)]
    cat(sprintf("%.2fs fit, dominant %s (%.0f%%)\n", sum(secs), top, 100 * max(p$table$share)))
    rows[[length(rows) + 1L]] <- data.frame(
      workload = wl$name, solver = solver, structure = fx$structure,
      n_levels = fx$n_levels, n_ind = fx$n_ind, relationship = fx$relationship,
      nEffects = p$counters$nEffects, nnzC = p$counters$nnzC,
      freeParams = p$counters$freeCovarianceParameters,
      numericFactorSeconds = p$counters$numericFactorSeconds,
      fitSeconds = sum(secs), elapsed = p$elapsed_sec, dominantPhase = top,
      as.data.frame(as.list(secs)), stringsAsFactors = FALSE)
  }
}

result <- do.call(rbind, rows)
utils::write.csv(result, output, row.names = FALSE)
print(result[, c("workload", "solver", "nEffects", "fitSeconds", "dominantPhase")], row.names = FALSE)
cat("Wrote", output, "\n")
