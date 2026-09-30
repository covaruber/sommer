oatsFit <- function(){
  data(DT_yatesoats, package="enhancer", envir=environment())
  mmes(Y ~ V + N + V:N, random=~ B + B:MP, rcov=~units, data=DT_yatesoats,
       verbose=FALSE)
}

oatsDtable <- function(m){
  Dt <- m$Dtable
  Dt[c(1, 2), "average"] <- TRUE
  Dt[c(3, 4), "include"] <- TRUE
  Dt[4, "average"] <- TRUE
  Dt
}

test_that("include-and-average interactions reproduce emmeans marginal means", {
  skip_if_not_installed("emmeans")
  m <- oatsFit()
  p <- predict(m, Dtable=oatsDtable(m), D="N")
  e <- summary(emmeans::emmeans(m, ~N))
  expect_equal(p$pvals$predicted.value, e$emmean, tolerance=1e-10)
  expect_equal(p$pvals$std.error, e$SE, tolerance=1e-10)
  expect_equal(p$pvals$N, c("0", "0.2", "0.4", "0.6"))
})

test_that("SED and pairwise comparisons follow from the prediction covariance", {
  skip_if_not_installed("emmeans")
  m <- oatsFit()
  p <- predict(m, Dtable=oatsDtable(m), D="N", sed=TRUE, pairwise=TRUE,
               adjust="bonferroni")
  A <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  expect_equal(c(p$sed[1, 2], p$sed[3, 4]),
               sqrt(diag(A %*% as.matrix(p$vcov) %*% t(A))), tolerance=1e-12)
  expect_equal(unname(diag(p$sed)), rep(0, 4))
  expect_equal(nrow(p$pairwise), 6L)
  pr <- summary(pairs(emmeans::emmeans(m, ~N)))
  expect_equal(p$pairwise$difference, pr$estimate, tolerance=1e-10)
  expect_equal(p$pairwise$SED, pr$SE, tolerance=1e-10)
  expect_equal(p$pairwise$p.adjusted, pmin(1, 6 * p$pairwise$p.value))
  expect_equal(p$avsed[["mean"]], mean(p$sed[upper.tri(p$sed)]))

  ref <- predict(m, Dtable=oatsDtable(m), D="N", pairwise="0", df="residual")$pairwise
  expect_equal(ref$level2, rep("0", 3))
  expect_equal(ref$df, rep(nrow(m$W) - 12, 3))
  expect_null(predict(m, Dtable=oatsDtable(m), D="N")$sed)
})

test_that("levels restrict averaging and prediction rows through the Dtable", {
  skip_if_not_installed("emmeans")
  m <- oatsFit()
  keep <- c("Marvellous", "Victory")
  cells <- as.vector(outer(keep, c("0", "0.2", "0.4", "0.6"), paste, sep=":"))
  p <- predict(m, Dtable=oatsDtable(m), D="N", levels=list(V=keep, "V:N"=cells))
  e <- summary(emmeans::emmeans(m, ~N, at=list(V=keep)))
  expect_equal(p$pvals$predicted.value, e$emmean, tolerance=1e-10)
  expect_equal(p$Dtable$levels[[2]], keep)

  rows <- predict(m, Dtable=oatsDtable(m), D="N", levels=list(N=c("0.6", "0.2")))
  expect_equal(rows$pvals$N, c("0.6", "0.2"))

  expect_error(predict(m, D="N", levels=list(foo=1)), "Valid terms")
  expect_error(predict(m, Dtable=oatsDtable(m), D="N", levels=list(V="Nope")),
               "Unknown levels")
})

test_that("levels can request random-effect levels without phenotypes", {
  set.seed(5)
  ids <- paste0("i", 1:12)
  G <- crossprod(matrix(rnorm(144), 12)) / 12 + diag(0.5, 12)
  dimnames(G) <- list(ids, ids)
  Gi <- solve(G); attr(Gi, "inverse") <- TRUE
  d <- data.frame(id=factor(rep(ids[1:10], each=3), levels=ids))
  d$y <- rnorm(30)
  m <- mmes(y~1, random=~vsm(ism(id), Gu=Gi), data=d, verbose=FALSE)
  term <- m$Dtable$term[2]
  p <- predict(m, D="id", levels=setNames(list(c("i1", "i12")), term))
  expect_equal(p$pvals$id, c("i1", "i12"))
  expect_equal(p$pvals$predicted.value[2],
               as.numeric(m$b[1] + m$u["i12", 1]), tolerance=1e-6)
})
