# Run against a deliberately installed package: Rscript benchmarks/run.R /path/to/output
suppressPackageStartupMessages(library(sommer))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1L || !nzchar(args)){
  stop("Usage: Rscript benchmarks/run.R OUTPUT_DIRECTORY", call.=FALSE)
}
output <- args[[1L]]
if(file.exists(output)){
  stop("Use a new output directory to avoid overwriting benchmark evidence.", call.=FALSE)
}
if(!dir.create(output, recursive=TRUE)){
  stop("Cannot create output directory: ", output, call.=FALSE)
}

workloads <- list(
  generate_mme_workload(name="serial", structure="ar2", n_ind=12L,
    n_levels=5L, seed=814L),
  generate_mme_workload(name="dense-relationship", structure="us", n_ind=16L,
    n_levels=3L, relationship="dense", seed=815L),
  generate_mme_workload(name="incomplete", structure="ar1", n_ind=16L,
    n_levels=4L, missing_rate=0.15, seed=816L),
  generate_mme_workload(name="factor-analytic", structure="fa", n_ind=16L,
    n_levels=4L, seed=817L)
)
names(workloads) <- vapply(workloads, function(workload) workload$name, character(1))
engines <- lapply(c("ldlt", "cholmod"), function(method)
  mmes_control(solve=list(method=method), diagnostics=list(profile="force")))
names(engines) <- c("ldlt", "cholmod")

metadata <- list(
  session=sessionInfo(), package=packageDescription("sommer"),
  library=find.package("sommer"),
  threads=Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "SLURM_CPUS_PER_TASK")),
  command=commandArgs(), system=Sys.info(),
  revision=if(nzchar(Sys.which("git")))
    system2("git", c("rev-parse", "HEAD"), stdout=TRUE) else NA_character_,
  worktree=if(nzchar(Sys.which("git")))
    system2("git", c("status", "--porcelain"), stdout=TRUE) else NA_character_
)
saveRDS(metadata, file.path(output, "environment.rds"))
result <- run_mme_benchmark(workloads, engines=engines, replications=3L,
  nIters=8L, tolerance=1e-5, verbose=TRUE)
saveRDS(result, file.path(output, "benchmark.rds"))
utils::write.csv(result$results, file.path(output, "results.csv"), row.names=FALSE)
print(result)
rows <- result$results
if(nrow(rows) != length(workloads) * length(engines) ||
   anyDuplicated(rows[c("workload", "engine")]) ||
   !setequal(rows$workload, names(workloads)) ||
   !setequal(rows$engine, names(engines))){
  stop("Benchmark did not return exactly one row per workload and engine.",
    call.=FALSE)
}
if(any(!is.na(rows$error)) || anyNA(rows$parity_ok) || any(!rows$parity_ok)){
  stop("Benchmark fitting or numerical parity failed; inspect the saved results.",
    call.=FALSE)
}
