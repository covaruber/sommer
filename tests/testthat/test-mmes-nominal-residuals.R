nominalResidualFixture <- function(){
  nominal <- sommer:::.laplace_nominal_response(
    factor(c("red","blue","green"),levels=c("red","blue","green")),
    reference="blue")
  logits <- matrix(c(.4,-.3,.7,.2,-.5,.1),nrow=3,
    dimnames=list(c("1","2","3"),nominal$predictorLevels))
  structure(list(family=multinomial("blue"),nominal=nominal,
    y=nominal$originalResponse,linear.predictors=logits,
    category.probabilities=sommer:::.laplace_nominal_probabilities(logits,nominal)),
    class=c("mmes.glmm.laplace","mmes.glmm","mmes"))
}

test_that("nominal residuals preserve labels, use conditional probabilities and do not mutate fits", {
  set.seed(27)
  object <- nominalResidualFixture()
  before <- object
  rng <- .Random.seed
  probability <- object$category.probabilities
  indicators <- diag(3)
  dimnames(indicators) <- dimnames(probability)
  expect_equal(residuals(object),indicators-probability,ignore_attr=TRUE)
  expect_equal(residuals(object,type="pearson"),
    (indicators-probability)/sqrt(probability*(1-probability)),ignore_attr=TRUE)
  expect_identical(dimnames(residuals(object)),dimnames(probability))
  expect_identical(dimnames(residuals(object,type="pearson")),dimnames(probability))
  deviance <- residuals(object,type="deviance")
  expect_equal(deviance,sqrt(-2*log(diag(probability))),ignore_attr=TRUE)
  expect_identical(names(deviance),rownames(probability))
  expect_true(all(deviance>=0))
  expect_error(residuals(object,type="working"),"unordered categories")
  expect_error(residuals(object,type="quantile"),"unordered categories")
  expect_identical(object,before)
  expect_identical(.Random.seed,rng)
})

test_that("nominal deviance magnitudes are invariant to reference and category order", {
  object <- nominalResidualFixture()
  expected <- residuals(object,type="deviance")
  logits <- cbind(object$linear.predictors,blue=0)
  recoded <- object
  recoded$y <- factor(as.character(object$y),levels=c("green","red","blue"))
  recoded$nominal <- sommer:::.laplace_nominal_response(recoded$y,"red")
  recoded$linear.predictors <- sweep(
    logits[,recoded$nominal$predictorLevels,drop=FALSE],1L,logits[,"red"],"-")
  expect_equal(residuals(recoded,type="deviance"),expected,tolerance=1e-12)
})

test_that("nominal deviance uses stable log probabilities without clipping tails", {
  object <- nominalResidualFixture()
  object$linear.predictors[,] <- 0
  object$linear.predictors[1L,] <- c(-1000,1000)
  object$category.probabilities <- sommer:::.laplace_nominal_probabilities(
    object$linear.predictors,object$nominal)
  expect_identical(object$category.probabilities[1L,"red"],0)
  expect_equal(unname(residuals(object,type="deviance")[1L]^2),4000,
    tolerance=1e-12)
  expect_error(residuals(object,type="pearson"),"zero or one")
  object$linear.predictors[1L,1L] <- NA_real_
  expect_error(residuals(object,type="deviance"),"finite logits")
  object <- nominalResidualFixture()
  object$y[1L] <- NA
  expect_error(residuals(object,type="deviance"),"valid category")
})
