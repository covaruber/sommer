# Crossed (full factorial) calibration design for the block-Schur dispatch rule.
#
# Purpose: choose Schur-vs-sparse-CHOLMOD thresholds (group size, block density,
# number of effects) from a design where the knobs vary independently instead of
# from unbalanced fixtures.
#
# Model: response ~ b1 + ... + b<border>
#        random   ~ vsm(ar1m(environment, rho), ism(id), Gu = precision)
#   group_size  -> nlevels(id) = dimension of the user precision Gu (one coupled
#                  block of Gu-size effects per environment level)
#   n_groups    -> nlevels(environment) = number of AR1-chained blocks
#                  (nEffects = group_size * n_groups)
#   density     -> fraction of nonzero off-diagonal entries of Gu (random sparse
#                  SPD precision, diagonally dominant; 1 = fully dense)
#   border      -> number of numeric fixed-effect covariates (engine reports
#                  blockSchurBorder = border + 1, the intercept included)
# Remaining confounding: (1) the AR1 coupling and replicates are held fixed, so
# conclusions hold for that structure only; (2) nEffects = group_size*n_groups
# and the number of records (2 replicates per id x environment cell) scale with
# the two size knobs; (3) with a random pattern, the effective Cholesky fill of
# the sparse path depends on the pattern, not only on the density; (4) the
# (4) the Schur contrast uses experimental$blockSchur force vs off; with "force"
# the engine may still be inactive when it is structurally or memory ineligible
# (blockSchurActive is recorded).
#
# Usage (from the package root):
#   OMP_NUM_THREADS=1 Rscript benchmarks/calibrate_paths_design.R [--quick] [--out=file.csv]
#     [--reps=5] [--stat=min|median]   (hidden: --smoke)

suppressPackageStartupMessages({
  library(sommer)
  library(Matrix)
})

# The Schur family is switched with experimental$blockSchur ("auto" is NOT used
# here: "force" vs "off" gives the cleanest on/off contrast; admission rules are
# calibrated from this contrast).
engine_on <- function(){
  mmes_control(solve=list(method="cholmod", blockSchur="force"),
  execution=list(denseMemoryMB=1000, memoryPolicy="limit"),
  diagnostics=list(profile="force")
)
}
engine_off <- function(){
  mmes_control(solve=list(method="cholmod", blockSchur="off"),
  execution=list(denseMemoryMB=1000, memoryPolicy="limit"),
  diagnostics=list(profile="force")
)
}

# pattern = "random": random sparse pattern (heavy Cholesky fill); "banded": band of
# half-width round(density * size / 2) (almost no fill) with matched nonzero fraction.
make_precision <- function(size, density, pattern="random", seed=1L){
  set.seed(seed)
  P <- matrix(0, size, size)
  upper <- which(upper.tri(P))
  edges <- max(size - 1L, as.integer(round(density * length(upper))))
  adjacent <- which(col(P) - row(P) == 1L)
  if(pattern == "banded"){
    keep <- head(upper[order(col(P)[upper] - row(P)[upper])], edges)
  } else {
    remaining <- setdiff(upper, adjacent)
    keep <- c(adjacent, if(edges > length(adjacent))
      sample(remaining, edges - length(adjacent)) else integer())
  }
  P[keep] <- runif(length(keep), -1, 1)
  P <- P + t(P)
  diag(P) <- rowSums(abs(P)) + 1
  dimnames(P) <- list(paste0("g", seq_len(size)), paste0("g", seq_len(size)))
  attr(P, "inverse") <- TRUE
  P
}

make_workload <- function(group_size, n_groups, density, border, pattern="random", rho=0.4, replicates=2L, seed=1L){
  set.seed(seed)
  precision <- make_precision(group_size, density, pattern, seed)
  set.seed(seed + 10000L)
  dataset <- expand.grid(id=factor(rownames(precision), levels=rownames(precision)),
                         environment=factor(seq_len(n_groups)),
                         replicate=seq_len(replicates))
  n <- nrow(dataset)
  covariates <- NULL
  if(border > 0L){
    covariates <- matrix(rnorm(n * border), n, border,
                         dimnames=list(NULL, paste0("b", seq_len(border))))
    dataset <- cbind(dataset, as.data.frame(covariates))
  }
  dataset$response <- rnorm(n) + if(border > 0L) covariates[, 1L] * 0.3 else 0
  fixed <- if(border > 0L) reformulate(paste0("b", seq_len(border)), "response") else response ~ 1
  list(data=dataset, precision=precision, fixed=fixed, rho=rho)
}

fit_once <- function(workload, engine){
  shape <- ar1m(workload$data$environment, rho=workload$rho)
  precision <- workload$precision
  start <- proc.time()[["elapsed"]]
  fit <- mmes(workload$fixed,
  random=~vsm(shape, ism(id), Gu=precision),
  data=workload$data,
  control=do.call(mmes_control, utils::modifyList(unclass(engine), list(solve=list(method="cholmod"), convergence=list(maxIters=4L), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
  list(time=proc.time()[["elapsed"]] - start, fit=fit)
}

run_cell <- function(cell, reps, stat){
  workload <- make_workload(cell$group_size, cell$n_groups, cell$density, cell$border,
                            pattern=as.character(cell$pattern), seed=cell$seed)
  summarize <- if(stat == "min") min else median
  results <- list(on=list(times=numeric()), off=list(times=numeric()))
  max_difference <- 0.0
  check_pair <- function(first, second){
    components <- c("llik", "covPar", "bu", "InfMat")
    differences <- vapply(components, function(component){
      first_value <- unlist(first[[component]])
      second_value <- unlist(second[[component]])
      if(!length(first_value) || length(first_value) != length(second_value) ||
         any(!is.finite(first_value)) || any(!is.finite(second_value))) return(Inf)
      max(abs(first_value - second_value) / (1 + abs(second_value)))
    }, numeric(1))
    if(any(differences > 1e-6) ||
       first$engineDiagnostics$coefficientEvaluations != second$engineDiagnostics$coefficientEvaluations){
      stop("Paired fits failed numerical or evaluation-count parity")
    }
    max(differences)
  }
  for(name in names(results)){
    invisible(fit_once(workload, if(name == "on") engine_on() else engine_off()))
  }
  for(repetition in seq_len(reps)){
    order <- if((repetition + cell$seed) %% 2L) c("on", "off") else c("off", "on")
    for(name in order){
      run <- fit_once(workload, if(name == "on") engine_on() else engine_off())
      results[[name]]$times[repetition] <- run$time
      results[[name]]$fit <- run$fit
    }
    max_difference <- max(max_difference, check_pair(results$on$fit, results$off$fit))
  }
  for(name in names(results)){
    results[[name]]$time <- summarize(results[[name]]$times)
  }
  diagnostics <- results$on$fit$engineDiagnostics
  pick <- function(field){
    value <- diagnostics[[field]]
    if(is.null(value) || length(value) != 1L) NA_real_ else as.numeric(value)
  }
  counters <- function(fit, field){
    value <- fit$engineDiagnostics$profile$counters[[field]]
    if(is.null(value) || length(value) != 1L) NA_real_ else as.numeric(value)
  }
  data.frame(cell,
    engine_subtype=if(isTRUE(diagnostics$blockChainActive)) "chain" else
      if(isTRUE(diagnostics$factorSchurActive)) "latent" else
      if(isTRUE(diagnostics$blockBandedActive)) "banded" else
      if(isTRUE(diagnostics$blockSchurActive)) "dense" else "inactive",
    parity_relative=max_difference,
    contrast_valid=isTRUE(diagnostics$blockSchurActive) &&
      !isTRUE(results$off$fit$engineDiagnostics$blockSchurActive),
    actual_input_fraction=(Matrix::nnzero(workload$precision) - cell$group_size) /
      (cell$group_size * (cell$group_size - 1)),
    time_on=results$on$time, time_off=results$off$time,
    sparseNnzL=counters(results$off$fit, "sparseFactorNonzeros"),
    sparseFlops=counters(results$off$fit, "sparseFactorFlops"),
    schurCost=counters(results$on$fit, "schurCostFlops"),
    blockSchurActive_on=pick("blockSchurActive"),
    blockSchurActive_off=as.numeric(isTRUE(results$off$fit$engineDiagnostics$blockSchurActive)),
    blockSchurGroups=pick("blockSchurGroups"), blockSchurBorder=pick("blockSchurBorder"),
    nEffects=nrow(results$off$fit$C),
    paired_ratios=paste(signif(results$off$times / results$on$times, 8), collapse=";"),
    llik_diff=abs(tail(results$on$fit$llik[1L, ], 1L) - tail(results$off$fit$llik[1L, ], 1L)))
}

args <- commandArgs(trailingOnly=TRUE)
option <- function(name, default){
  hit <- grep(paste0("^--", name, "="), args, value=TRUE)
  if(length(hit)) sub(paste0("^--", name, "="), "", hit[1L]) else default
}
smoke <- "--smoke" %in% args
quick <- "--quick" %in% args
out <- option("out", "/tmp/calibrate_paths_design.csv")
reps <- as.integer(option("reps", if(smoke) 1L else 5L))
stat <- option("stat", "median")
workers <- as.integer(option("workers", "1"))
stopifnot(is.finite(workers), workers >= 1L, reps >= 1L, stat %in% c("min", "median"))
thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
if(workers > 1L && any(Sys.getenv(thread_variables) != "1")){
  stop("Parallel runs require OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and BLIS_NUM_THREADS all set to 1 before R starts.")
}

levels <- list(group_size=c(8, 16, 32, 64, 128, 256), n_groups=c(10, 50, 200),
               density=c(0.02, 0.1, 0.5, 1), border=c(5, 50, 200))
if(quick) levels <- list(group_size=c(16, 64, 128), n_groups=c(10, 50),
                         density=c(0.1, 1), border=c(5, 50))
# medium: adds the fill-free banded pattern and density 0.02; (256 x 50) is dropped for run time
if("--medium" %in% args){
  levels <- list(group_size=c(16, 32, 64, 128, 256), n_groups=c(10, 50),
                 density=c(0.02, 0.1, 0.5), border=c(5, 50), pattern=c("random", "banded"))
}
if(smoke) levels <- list(group_size=c(4, 8), n_groups=c(3, 4), density=c(0.5, 1), border=c(2, 3))
if(is.null(levels$pattern)) levels$pattern <- c("random", "banded")
levels$seed <- if(smoke) 1L else c(101L, 202L, 303L)
design <- do.call(expand.grid, c(levels, list(KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)))
if(smoke) design <- design[c(1L, 2L, 17L, 18L), ]
cat(sprintf("cells=%d reps=%d stat=%s fits=%d\n", nrow(design), reps, stat, 2L * reps * nrow(design)))

checkpoint <- paste0(out, ".cells")
dir.create(checkpoint, recursive=TRUE, showWarnings=FALSE)
manifest <- list(design=design, reps=reps, stat=stat,
                 workers=workers, script_md5=unname(tools::md5sum("benchmarks/calibrate_paths_design.R")),
                 binary_md5=unname(tools::md5sum(list.files(file.path(find.package("sommer"), "libs"),
                                                          pattern="\\.so$", full.names=TRUE))),
                 package=as.character(packageVersion("sommer")),
                 library=find.package("sommer"), session=capture.output(sessionInfo()),
                 threads=Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")))
manifest_path <- file.path(checkpoint, "manifest.rds")
if(file.exists(manifest_path)){
  previous <- readRDS(manifest_path)
    if(!identical(previous$design, design) || previous$reps != reps || previous$stat != stat ||
      !identical(previous$workers, workers) || !identical(previous$binary_md5, manifest$binary_md5) ||
      !identical(previous$script_md5, manifest$script_md5) || !identical(previous$threads, manifest$threads)){
    stop("Checkpoint design/settings differ; select a new --out path.")
  }
} else saveRDS(manifest, manifest_path)
if("--audit" %in% args){
  print(lapply(levels, identity))
  stopifnot(nrow(design) == prod(lengths(levels)), !anyDuplicated(design))
  stopifnot(all(2 * design$group_size * design$n_groups > design$border + 1L))
  for(size in unique(design$group_size)) for(fraction in unique(design$density)){
    random <- make_precision(size, fraction, "random", 101L)
    banded <- make_precision(size, fraction, "banded", 101L)
    stopifnot(Matrix::nnzero(random) == Matrix::nnzero(banded))
    invisible(chol(random))
    invisible(chol(banded))
  }
  cat("Matched edge counts and positive definiteness verified. Fixed-design dimensions feasible.\n")
  cat("Complete factorial verified. No fits run.\n")
  quit(status=0)
}
execute <- function(index){
  path <- file.path(checkpoint, sprintf("cell-%05d.rds", index))
  if(file.exists(path)) return(readRDS(path))
  result <- run_cell(design[index, ], reps, stat)
  temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(result, temporary)
  if(!file.rename(temporary, path)) stop("Failed to save cell checkpoint")
  cat(sprintf("[%d/%d] %s ratio=%.3f\n", index, nrow(design),
              design$pattern[index], result$time_off / result$time_on))
  result
}
set.seed(817L)
execution_order <- sample(seq_len(nrow(design)))
rows <- parallel::mclapply(execution_order, execute,
                          mc.cores=workers, mc.preschedule=FALSE, mc.set.seed=FALSE)
if(any(vapply(rows, inherits, logical(1), "try-error"))) stop("Cell failed; completed cells remain checkpointed.")
result <- do.call(rbind, rows)
write.csv(result, out, row.names=FALSE)
cat("Engine subtype counts:\n")
print(table(result$engine_subtype))
cat("Valid contrasts:", sum(result$contrast_valid), "of", nrow(result), "\n")
cat("max |llik diff|:", format(max(result$llik_diff, na.rm=TRUE)), "\n")

result$log_ratio <- log(result$time_off / result$time_on)
if("pattern" %in% names(result)){
  model <- lm(log_ratio ~ pattern * log2(group_size) + density + log(n_groups) + log(border), data=result)
} else {
  model <- lm(log_ratio ~ log2(group_size) * density + log(n_groups) + log(border), data=result)
}
print(summary(model)$coefficients)
ok <- is.finite(result$sparseFlops) & is.finite(result$schurCost) & result$schurCost > 0
if(sum(ok) >= 4L){
  result$log_predicted <- log(result$sparseFlops / result$schurCost)
  fill_model <- lm(log_ratio ~ log_predicted, data=result[ok, ])
  cat(sprintf("predicted fill ratio vs measured: R2 = %.2f, slope = %.2f, break-even predicted ratio = %.2e (cells=%d)\n",
              summary(fill_model)$r.squared, coef(fill_model)[2L],
              exp(-coef(fill_model)[1L] / coef(fill_model)[2L]), sum(ok)))
}
cat(sprintf("ratio > 1 (Schur faster) in %d of %d cells\n", sum(result$log_ratio > 0), nrow(result)))
for(d in sort(unique(result$density))){
  sub <- result[result$density == d, ]
  if(nrow(sub) < 2L) next
  slope <- tryCatch(coef(lm(log_ratio ~ log2(group_size), data=sub)), error=function(e) c(NA, NA))
  crossing <- if(is.finite(slope[2L]) && slope[2L] != 0) 2^(-slope[1L] / slope[2L]) else NA_real_
  cat(sprintf("density %g: ratio crosses 1 near group_size %s\n", d,
              if(is.finite(crossing)) format(round(crossing, 1)) else "n/a (no crossing)"))
}
cat("Results written to", out, "\n")
