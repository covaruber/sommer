# Sweep Axis C: D-First ML Schur vs Baseline Joint C Factorization
.libPaths(c("/tmp/sommer-review-implementation-lib", .libPaths()))
library(sommer)

cat("=== Axis C Sweep: D-First ML Schur vs Joint C ML ===\n")

p_vec <- c(2L, 8L, 25L)
q_vec <- c(60L, 200L, 600L)
solvers <- c("ldlt", "cholmod")

results <- list()

for(solv in solvers){
  for(p in p_vec){
    for(q in q_vec){
      fixture_name <- sprintf("ML_%s_%dp_%dq", solv, p, q)
      cat(sprintf("Running %-25s ... ", fixture_name))
      
      set.seed(42L)
      # Build dataset with p fixed levels and q random levels
      df <- expand.grid(
        fixed_grp = factor(seq_len(p)),
        id = factor(seq_len(q))
      )
      df$response <- stats::rnorm(nrow(df)) + as.integer(df$fixed_grp) * 0.5 + stats::rnorm(q)[df$id]
      
      M <- matrix(stats::rnorm(q * 25L), q, 25L)
      G <- stats::cor(t(M)) + diag(0.1, q)
      dimnames(G) <- list(levels(df$id), levels(df$id))
      precision <- solve(G)
      attr(precision, "inverse") <- TRUE
      
      fit_base <- tryCatch({
        t0 <- proc.time()
        m_base <- mmes(response ~ fixed_grp,
  random=~ vsm(ism(id), Gu = precision),
  data=df,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(solve=list(method=solv, dFirstML="off")
)), list(solve=list(method=solv), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=4L), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
        t1 <- proc.time()
        list(time = (t1 - t0)[3L], fit = m_base)
      }, error = function(e) list(error = e$message))
      
      fit_dfirst <- tryCatch({
        t0 <- proc.time()
        m_df <- mmes(response ~ fixed_grp,
  random=~ vsm(ism(id), Gu = precision),
  data=df,
  control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(solve=list(method=solv, dFirstML="force")
)), list(solve=list(method=solv), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=4L), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))
)
        t1 <- proc.time()
        list(time = (t1 - t0)[3L], fit = m_df)
      }, error = function(e) list(error = e$message))
      
      if(is.null(fit_base$error) && is.null(fit_dfirst$error)){
        t_base <- fit_base$time
        t_df   <- fit_dfirst$time
        speedup <- t_base / t_df
        active <- fit_dfirst$fit$engineDiagnostics$dFirstMLActive
        
        rel_llik <- abs(tail(fit_dfirst$fit$llik, 1L) - tail(fit_base$fit$llik, 1L)) / (abs(tail(fit_base$fit$llik, 1L)) + 1.0)
        parity <- (rel_llik < 1e-4)
        
        cat(sprintf("Base: %.4fs | DFirst: %.4fs | Speedup: %.2fx | Active: %s | Parity: %s\n",
                    t_base, t_df, speedup, active, parity))
        
        results[[length(results) + 1L]] <- data.frame(
          solver = solv,
          p_fixed = p,
          q_random = q,
          t_base = t_base,
          t_df = t_df,
          speedup = speedup,
          active = active,
          parity = parity
        )
      } else {
        cat("ERROR:", fit_base$error, fit_dfirst$error, "\n")
      }
    }
  }
}

df_res_c <- do.call(rbind, results)
print(df_res_c)
saveRDS(df_res_c, "/tmp/sweep_axis_c.rds")
