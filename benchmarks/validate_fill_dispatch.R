suppressPackageStartupMessages(library(sommer))
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2L) stop("usage: validate_fill_dispatch.R label output.rds")
label <- args[1L]
output <- args[2L]

make_graph_workload <- function(pattern, seed=303L){
  size <- 256L
  groups <- 10L
  set.seed(seed)
  precision <- matrix(0, size, size)
  upper <- which(upper.tri(precision))
  edges <- round(0.02 * length(upper))
  adjacent <- which(col(precision) - row(precision) == 1L)
  keep <- if(pattern == "random"){
    c(adjacent, sample(setdiff(upper, adjacent), edges - length(adjacent)))
  } else head(upper[order(col(precision)[upper] - row(precision)[upper])], edges)
  precision[keep] <- runif(length(keep), -1, 1)
  precision <- precision + t(precision)
  diag(precision) <- rowSums(abs(precision)) + 1
  dimnames(precision) <- rep(list(paste0("g", seq_len(size))), 2L)
  attr(precision, "inverse") <- TRUE
  set.seed(seed + 10000L)
  dataset <- expand.grid(id=factor(rownames(precision), levels=rownames(precision)),
                         environment=factor(seq_len(groups)), replicate=1:2)
  covariates <- matrix(rnorm(nrow(dataset) * 5L), ncol=5L,
                       dimnames=list(NULL, paste0("b", 1:5)))
  dataset <- cbind(dataset, as.data.frame(covariates))
  dataset$response <- rnorm(nrow(dataset)) + 0.3 * covariates[, 1L]
  shape <- ar1m(dataset$environment, rho=0.4)
  list(name=paste0("graph_", pattern), data=dataset,
       fixed=response ~ b1 + b2 + b3 + b4 + b5,
       random=~vsm(shape, ism(id), Gu=precision), rcov=~units)
}

workloads <- list(make_graph_workload("random"), make_graph_workload("banded"))
for(configuration in list(
  list("us", 5L, 100L, "sparse"), list("us", 5L, 100L, "dense"),
  list("ar1", 5L, 100L, "dense"), list("ar1", 50L, 64L, "none"),
  list("ar2", 30L, 8L, "none"), list("sar", 100L, 4L, "none"),
  list("car", 100L, 4L, "none"))){
  workloads[[length(workloads) + 1L]] <- generate_mme_workload(
    structure=configuration[[1L]], n_levels=configuration[[2L]],
    n_ind=configuration[[3L]], relationship=configuration[[4L]], seed=303L)
}

fit_model <- function(workload, mode, inverse_mode=0L){
  mmes(workload$fixed, random=workload$random, rcov=workload$rcov,
       data=workload$data, solver="cholmod", nIters=4L, computeCi=inverse_mode,
       getPEV=inverse_mode > 0L, verbose=FALSE, dateWarning=FALSE,
       engine=mmes_control(solver="cholmod", memoryPolicy="limit",
         experimental=list(blockSchur=mode, profile="force")))
}

signature <- function(fit){
  fields <- c("llik", "covPar", "bu", "InfMat", "dLu", "theta_se", "uPevList")
  result <- lapply(fields, function(field) unlist(fit[[field]]))
  names(result) <- fields
  if(isTRUE(fit$CiMode == 2L)) result$Ci <- as.numeric(as.matrix(fit$Ci))
  result
}

compare <- function(first, second){
  errors <- vapply(names(first), function(field){
    reference <- second[[field]]
    candidate <- first[[field]]
    if(length(reference) != length(candidate)) return(Inf)
    if(!length(reference)) return(0)
    if(any(!is.finite(reference)) || any(!is.finite(candidate))) return(Inf)
    max(abs(reference - candidate) / (1 + abs(reference)))
  }, numeric(1))
  if(any(errors > 1e-6)) stop("Parity failure: ", paste(names(errors)[errors > 1e-6], collapse=", "))
  max(errors)
}

rows <- list()
signatures <- list()
for(workload in workloads){
  latest <- list()
  times <- list(auto=numeric(3L), off=numeric(3L))
  for(mode in names(times)) invisible(fit_model(workload, mode))
  for(repetition in 1:3){
    modes <- if(repetition %% 2L) c("auto", "off") else c("off", "auto")
    for(mode in modes){
      started <- proc.time()[["elapsed"]]
      latest[[mode]] <- fit_model(workload, mode)
      times[[mode]][repetition] <- proc.time()[["elapsed"]] - started
    }
    difference <- compare(signature(latest$auto), signature(latest$off))
  }
  diagnostics <- latest$auto$engineDiagnostics
  counters <- diagnostics$profile$counters
  rows[[length(rows) + 1L]] <- data.frame(
    label=label, workload=workload$name, auto_seconds=median(times$auto),
    off_seconds=median(times$off), off_over_auto=median(times$off) / median(times$auto),
    block_active=isTRUE(diagnostics$blockSchurActive), chain=isTRUE(diagnostics$blockChainActive),
    sparse_flops=counters$sparseFactorFlops, schur_cost=counters$schurCostFlops,
    max_relative_difference=difference)
  signatures[[workload$name]] <- signature(latest$auto)
  cat(label, workload$name, "auto", median(times$auto), "off", median(times$off),
      "active", diagnostics$blockSchurActive, "chain", diagnostics$blockChainActive, "\n")
}

small <- generate_mme_workload(structure="us", n_levels=3L, n_ind=12L,
                               relationship="sparse", seed=303L)
for(inverse_mode in 0:2){
  auto <- fit_model(small, "auto", inverse_mode)
  off <- fit_model(small, "off", inverse_mode)
  compare(signature(auto), signature(off))
  signatures[[paste0("inverse_", inverse_mode)]] <- signature(auto)
}
result <- list(results=do.call(rbind, rows), signatures=signatures,
               compare=compare, library=find.package("sommer"), session=sessionInfo())
saveRDS(result, output)
write.csv(result$results, sub("\\.rds$", ".csv", output), row.names=FALSE)
cat("VALIDATION_OK\n")