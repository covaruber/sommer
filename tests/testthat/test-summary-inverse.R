test_that("summary preserves fixed-effect covariance dimensions", {
  set.seed(193)
  dat <- expand.grid(Genotype=factor(seq_len(6)),
    Block=factor(seq_len(3)), Year=factor(seq_len(4)))
  dat$Yield <- 5 + rnorm(nrow(dat))
  dat$x <- rnorm(nrow(dat))
  random <- ~ Year + Genotype + Genotype:Year + Block:Year

  for(fixed in list(Yield ~ 1, Yield ~ x)){
    for(inverse in 0:2){
      fit <- mmes(fixed, random=random, data=dat,
        control=mmes_control(diagnostic=list(verbose=FALSE),
          output=list(inverse=inverse)))
      result <- summary(fit)
      expect_equal(result$varcomp,
        fit$covParNative[, c("term", "parameter", "estimate", "StdError", "Zratio")])
      if(inverse == 2){
        n <- length(fit$b)
        expected <- sqrt(abs(diag(as.matrix(fit$Ci[seq_len(n), seq_len(n),
          drop=FALSE]))))
        expect_equal(unname(result$betas$Std.Error), unname(expected))
      }else{
        expect_true(all(is.na(result$betas$Std.Error)))
      }
    }
  }
})
