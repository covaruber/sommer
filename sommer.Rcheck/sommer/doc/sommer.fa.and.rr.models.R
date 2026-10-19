## -----------------------------------------------------------------------------
library(sommer)
data(DT_h2, package="enhancer")
DT <- DT_h2
DT <- DT[with(DT, order(Env)), ]
length(unique(DT$Env))
length(unique(DT$Name))
head(DT)

## -----------------------------------------------------------------------------
fitFA <- mmes(y ~ Env,
  random = ~ vsm(fam(Env, 2), ism(Name)),
  rcov = ~ units,
  data = DT,
  control = mmes_control(
    # convergence = list(maxIters = 150),
  diagnostics = list(verbose = FALSE)))
summary(fitFA)$varcomp

## -----------------------------------------------------------------------------
fitRR <- mmes(y ~ Env,
  random = ~ vsm(rrm(Env, 2), ism(Name)),
  rcov = ~ units,
  data = DT,
  control = mmes_control(
    # convergence = list(maxIters = 150),
  diagnostics = list(verbose = FALSE)))
summary(fitRR)$varcomp

## -----------------------------------------------------------------------------
c(logLik_FA = tail(fitFA$llik, 1), logLik_RR = tail(fitRR$llik, 1))

## -----------------------------------------------------------------------------
faInfo <- loadings_mmes(fitFA, varianceScale = FALSE, rotation = TRUE)
round(faInfo$loadings, 3)
round(faInfo$specific, 3)
faInfo$sigma2

## -----------------------------------------------------------------------------
rrInfo <- loadings_mmes(fitRR, varianceScale = FALSE, rotation = TRUE)
round(rrInfo$specific, 3)

## -----------------------------------------------------------------------------
faScores <- scores_mmes(fitFA, method = "regression",
                       varianceScale = TRUE, rotation = TRUE)
head(faScores)

## -----------------------------------------------------------------------------
faUnrotated <- loadings_mmes(fitFA, varianceScale = FALSE,
                           rotation = FALSE)
scoresUnrotated <- scores_mmes(fitFA, method = "regression",
                              varianceScale = TRUE, rotation = FALSE)
L <- sqrt(faInfo$sigma2) * faInfo$loadings
LUnrotated <- sqrt(faUnrotated$sigma2) * faUnrotated$loadings
common <- faScores %*% t(L)
commonUnrotated <- scoresUnrotated %*% t(LUnrotated)

stopifnot(
  isTRUE(all.equal(tcrossprod(L), tcrossprod(LUnrotated))),
  isTRUE(all.equal(common, commonUnrotated)),
  isTRUE(all.equal(faInfo$specific, faUnrotated$specific))
)
round(crossprod(L), 3)

## -----------------------------------------------------------------------------
varPerFactor <- colSums(faInfo$loadings^2)
totalVar <- sum(faInfo$loadings^2) + sum(faInfo$specific)
propExplained <- varPerFactor / totalVar
round(100 * propExplained, 1)

## ----fig.show='hold'----------------------------------------------------------
barplot(t(faInfo$loadings), beside = TRUE,
        col = c("steelblue4", "tomato"),
        las = 2, cex.names = 0.7,
        ylab = "Loading",
        main = "Factor-analytic loadings by environment")
legend("topright", legend = colnames(faInfo$loadings),
       fill = c("steelblue4", "tomato"), bty = "n")

## ----fig.show='hold'----------------------------------------------------------
plot(faScores[,2] ~ faScores[,1],
     xlab = "Factor 1 score", ylab = "Factor 2 score",
     main = "Genotype scores")
text(faScores[,2] ~ faScores[,1], labels = rownames(faScores),
     cex = 0.6, pos = 1)
abline(h = 0, v = 0, lty = 3)

## ----fig.width=9, fig.height=5------------------------------------------------
U <- fitFA$uList[[faInfo$term]]
U <- U[rownames(faScores), rownames(L), drop = FALSE]
specificBLUP <- U - common
scoreOrder <- order(faScores[, 1])
selected <- scoreOrder[unique(round(seq(1, length(scoreOrder),
                                       length.out = 4)))]
genotypeColors <- c("steelblue4", "tomato3", "darkgreen", "goldenrod3")

local({
  oldPar <- par(mfrow = c(1, ncol(L)), mar = c(4.5, 4.5, 3, 1))
  on.exit(par(oldPar))
  for (factorIndex in seq_len(ncol(L))) {
    loading <- L[, factorIndex]
    factorResponse <- outer(faScores[selected, factorIndex], loading)
    adjusted <- factorResponse + specificBLUP[selected, , drop = FALSE]
    plot(range(loading), range(c(adjusted, factorResponse, 0)),
         type = "n", xlab = paste("Environment loading:", colnames(L)[factorIndex]),
         ylab = "Adjusted genetic effect", main = paste("Latent regression", factorIndex))
    abline(h = 0, v = 0, lty = 3, col = "grey70")
    for (genotypeIndex in seq_along(selected)) {
      points(loading, adjusted[genotypeIndex, ], pch = 16,
             col = genotypeColors[genotypeIndex])
      abline(a = 0, b = faScores[selected[genotypeIndex], factorIndex],
             col = genotypeColors[genotypeIndex], lwd = 2)
    }
    legend("topleft", legend = rownames(faScores)[selected],
           col = genotypeColors[seq_along(selected)], pch = 16, lty = 1,
           cex = 0.7, bty = "n")
  }
})

## ----fig.show='hold'----------------------------------------------------------
Sigma <- faInfo$sigma2 * (
  faInfo$loadings %*% t(faInfo$loadings) + diag(faInfo$specific)
)
corMat <- cov2cor(Sigma)
stopifnot(isTRUE(all.equal(unname(Sigma),
                          unname(fitFA$theta[[faInfo$term]]))))
heatmap(corMat, symm = TRUE,
        main = "Fitted genetic correlation among environments")

