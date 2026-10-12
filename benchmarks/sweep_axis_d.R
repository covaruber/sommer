# Sweep Axis D: Panel Full Inverse (panelCi2) & Matrix-Free PCG
.libPaths(c("/tmp/sommer-review-implementation-lib", .libPaths()))
library(sommer)

cat("=== Axis D Sweep: panelCi2 Full Inverse Batch Sizes ===\n")

q_sizes <- c(300L, 800L)
panel_sizes <- c(16L, 32L, 64L, 128L)

panel_results <- list()

for(q in q_sizes){
  set.seed(42L)
  df <- data.frame(
    id = factor(paste0("ID", seq_len(q))),
    response = stats::rnorm(q)
  )
  M <- matrix(stats::rnorm(q * 20L), q, 20L)
  G <- stats::cor(t(M)) + diag(0.1, q)
  dimnames(G) <- list(levels(df$id), levels(df$id))
  precision <- solve(G)
  attr(precision, "inverse") <- TRUE
  
  for(ps in panel_sizes){
    cat(sprintf("Running computeCi=2, N=%d, panel=%-3d ... ", q + 1L, ps))
    t0 <- proc.time()
    m_fit <- mmes(
      response ~ 1,
      random = ~ vsm(ism(id), Gu = precision),
      data = df,
      solver = "cholmod",
      nIters = 2L,
      computeCi = 2L,
      getPEV = FALSE,
      verbose = FALSE,
      dateWarning = FALSE,
      engine = mmes_control(solver = "cholmod", experimental = list(panelCi2 = ps))
    )
    t1 <- proc.time()
    elapsed <- (t1 - t0)[3L]
    cat(sprintf("Elapsed: %.4fs | Ci dim: %dx%d\n", elapsed, nrow(m_fit$Ci), ncol(m_fit$Ci)))
    panel_results[[length(panel_results) + 1L]] <- data.frame(
      N = q + 1L,
      panel = ps,
      elapsed_sec = elapsed
    )
  }
}

df_panel <- do.call(rbind, panel_results)
print(df_panel)

cat("\n=== Axis D Sweep: Matrix-Free PCG vs Assembled PCG ===\n")
mf_results <- list()
pcg_sizes <- c(400L, 1200L, 2500L)

for(q in pcg_sizes){
  cat(sprintf("Running PCG N=%-4d ... ", q))
  set.seed(42L)
  df <- data.frame(
    id = factor(paste0("ID", seq_len(q))),
    response = stats::rnorm(q)
  )
  M <- matrix(stats::rnorm(q * 15L), q, 15L)
  G <- stats::cor(t(M)) + diag(0.1, q)
  dimnames(G) <- list(levels(df$id), levels(df$id))
  precision <- solve(G)
  attr(precision, "inverse") <- TRUE
  
  # Assembled PCG
  t0 <- proc.time()
  m_asm <- mmes(
    response ~ 1,
    random = ~ vsm(ism(id), Gu = precision),
    data = df,
    solver = "pcg",
    nIters = 3L,
    computeCi = 0L,
    getPEV = FALSE,
    verbose = FALSE,
    dateWarning = FALSE,
    engine = mmes_control(solver = "pcg", experimental = list(matrixFreePcg = "off"))
  )
  t_asm <- (proc.time() - t0)[3L]
  
  # Matrix-Free PCG
  t0 <- proc.time()
  m_mf <- mmes(
    response ~ 1,
    random = ~ vsm(ism(id), Gu = precision),
    data = df,
    solver = "pcg",
    nIters = 3L,
    computeCi = 0L,
    getPEV = FALSE,
    verbose = FALSE,
    dateWarning = FALSE,
    engine = mmes_control(solver = "pcg", experimental = list(matrixFreePcg = "force"))
  )
  t_mf <- (proc.time() - t0)[3L]
  
  speedup <- t_asm / t_mf
  rel_llik <- abs(tail(m_mf$llik, 1L) - tail(m_asm$llik, 1L)) / (abs(tail(m_asm$llik, 1L)) + 1.0)
  cat(sprintf("Asm: %.4fs | MatFree: %.4fs | Speedup: %.2fx | Parity: %s\n",
              t_asm, t_mf, speedup, (rel_llik < 1e-4)))
  mf_results[[length(mf_results) + 1L]] <- data.frame(
    N = q,
    t_asm = t_asm,
    t_mf = t_mf,
    speedup = speedup
  )
}

df_mf <- do.call(rbind, mf_results)
print(df_mf)
saveRDS(list(panel = df_panel, matfree = df_mf), "/tmp/sweep_axis_d.rds")
