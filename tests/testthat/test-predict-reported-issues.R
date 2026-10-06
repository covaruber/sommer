prediction_issue_data <- function(){
  data(list="DT_example", package="enhancer", envir=environment())
  data <- droplevels(get("DT_example", envir=environment()))
  data$Env <- factor(data$Env)
  data$Name <- factor(data$Name)
  set.seed(1)
  data$cov <- rnorm(nrow(data), mean=5, sd=0.5)
  data
}

test_that("response scaling preserves no-intercept GLS fits in both engines", {
  data <- prediction_issue_data()
  for(henderson in c(TRUE, FALSE)){
    for(REML in c(TRUE, FALSE)){
      for(fixed in list(Yield~Env-1+cov:Env, Yield~Env+cov:Env, Yield~cov-1)){
        model <- mmes(fixed, random=~Name, data=data, verbose=FALSE,
                      henderson=henderson, REML=REML)
        variance <- summary(model)$varcomp$estimate
        X <- model.matrix(fixed, data)
        Z <- model.matrix(~Name-1, data)
        V <- variance[1] * tcrossprod(Z) + diag(variance[2], nrow(data))
        beta <- solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, data$Yield)))
        expect_equal(as.numeric(model$b), as.numeric(beta), tolerance=1e-7)
        expect_equal(as.numeric(model$bu[seq_len(nrow(model$b)),1]),
                     as.numeric(model$b), tolerance=1e-12)
      }
    }
  }
})

test_that("no-intercept rotation agrees with unrotated fits", {
  set.seed(12)
  ids <- paste0("g", seq_len(6))
  data <- expand.grid(id=ids, env=c("e1", "e2"), KEEP.OUT.ATTRS=FALSE)
  data$cov <- rnorm(nrow(data), 1, 0.4)
  data$y <- 4 + data$cov + rnorm(nrow(data))
  relationship <- outer(seq_along(ids), seq_along(ids), function(first, second) 0.5^abs(first-second))
  dimnames(relationship) <- list(ids, ids)
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  for(henderson in c(TRUE, FALSE)){
    for(fixed in list(y~env-1+cov:env, y~cov-1)){
      ordinary <- mmes(fixed, random=~vsm(ism(id), Gu=precision),
                       data=data, verbose=FALSE, henderson=henderson)
      rotated <- mmes(fixed, random=~vsm(ism(id), Gu=precision, rotation=TRUE),
                      data=data, verbose=FALSE, henderson=henderson)
      expect_equal(as.numeric(rotated$b), as.numeric(ordinary$b), tolerance=1e-6)
      expect_equal(unname(rotated$theta), unname(ordinary$theta), tolerance=1e-6)
      expect_equal(rotated$uList[[1]], ordinary$uList[[1]], tolerance=1e-6)
    }
  }
})

test_that("structured predictions use fitted genetic effects for empty cells", {
  data <- prediction_issue_data()
  counts <- table(data$Name, data$Env)
  for(shape in c("usm", "dsm", "corgm", "fam")){
    random <- as.formula(paste0("~vsm(", shape, "(Env), ism(Name))"))
    model <- mmes(Yield~Env, random=random, data=data, verbose=FALSE)
    prediction <- predict(model, D="Env:Name")
    ids <- prediction$pvals[[1]]
    env <- sub(":.*", "", ids)
    name <- sub("^[^:]*:", "", ids)
    empty <- counts[cbind(name, env)] == 0
    expect_equal(sum(empty), 29L)
    effects <- model$uList[[1]]
    ranges <- model$partitions[[1]]
    columns <- ranges[match(env, colnames(effects)),1] + match(name, rownames(effects)) - 1L
    geneticColumns <- seq.int(nrow(model$b)+1L, ncol(prediction$D))
    expectedD <- prediction$D
    expectedD[,geneticColumns] <- 0
    expectedD[cbind(seq_along(ids), columns)] <- 1
    expected <- predict(model, D=expectedD)
    expect_equal(prediction$D, expectedD, tolerance=1e-12)
    expect_equal(prediction$pvals$predicted.value, expected$pvals$predicted.value)
    expect_equal(prediction$vcov, expected$vcov)
    expect_equal(as.numeric(Matrix::rowSums(prediction$D[empty,geneticColumns,drop=FALSE])),
                 rep(1, sum(empty)))
    requested <- ids[c(which(empty)[1L], which(!empty)[1L])]
    selected <- predict(model, D="Env:Name",
                        levels=setNames(list(requested), names(model$uList)))
    expect_equal(selected$pvals[[1]], requested)
    expect_equal(selected$D, prediction$D[requested,,drop=FALSE])
    reversed <- predict(model, D="Name:Env")
    reversedIds <- paste(name, env, sep=":")
    expect_equal(unname(as.matrix(reversed$D[reversedIds,,drop=FALSE])),
                 unname(as.matrix(prediction$D)))
    marginal <- predict(model, D="Name")
    expect_equal(as.numeric(Matrix::rowSums(marginal$D[,geneticColumns,drop=FALSE])),
                 rep(1, nrow(marginal$D)))
    expect_equal(unique(as.numeric(marginal$D[,geneticColumns])), c(1/3, 0))
  }
})

test_that("structured cell maps include Gu-only levels and remain term local", {
  data <- prediction_issue_data()
  ids <- c(levels(data$Name), "new-genotype")
  precision <- diag(length(ids))
  dimnames(precision) <- list(ids, ids)
  attr(precision, "inverse") <- TRUE
  model <- mmes(Yield~Env, random=~Name+vsm(usm(Env), ism(Name), Gu=precision),
                data=data, verbose=FALSE)
  term <- names(model$uList)[2L]
  requested <- "CA.2011:new-genotype"
  prediction <- predict(model, D="Env:Name", levels=setNames(list(requested), term))
  range <- model$partitions[[2L]][1L,]
  column <- range[1L] + match("new-genotype", rownames(model$uList[[2L]])) - 1L
  expect_equal(as.numeric(prediction$D[,column]), 1)
  otherRange <- model$partitions[[1L]][1L,]
  expect_equal(as.numeric(prediction$D[,seq.int(otherRange[1L], otherRange[2L])]),
               rep(0, diff(otherRange)+1L))
})

test_that("mixed terms use numeric means and preserve supplied numeric values", {
  data <- prediction_issue_data()
  for(fixed in list(Yield~Env+cov:Env, Yield~Env-1+cov:Env)){
    model <- mmes(fixed, random=~Name, data=data, verbose=FALSE)
    term <- names(model$partitionsX)[grepl("cov", names(model$partitionsX))]
    columns <- as.vector(model$partitionsX[[term]])
    prediction <- predict(model, D="Name")
    expect_equal(unname(as.matrix(prediction$D[,columns])),
                 matrix(mean(data$cov)/3, nrow(prediction$D), 3))
    expect_equal(prediction$Dtable$levels[[match(term, prediction$Dtable$term)]]$cov,
                 mean(data$cov))
    for(value in c(-2, 0, 5)){
      evaluated <- predict(model, D="Name", levels=setNames(list(value), term))
      expect_equal(unname(as.matrix(evaluated$D[,columns])),
                   matrix(value/3, nrow(prediction$D), 3))
    }
    selected <- predict(model, D="Name", levels=setNames(list(list(Env="CA.2012", cov=-2)), term))
    expect_equal(unname(as.matrix(selected$D[,columns])),
                 matrix(rep(c(0, -2, 0), each=nrow(selected$D)), nrow(selected$D)))
    expect_equal(predict(model, D="Name", Dtable=selected$Dtable)$D, selected$D)
    expect_error(predict(model, D="Name", levels=setNames(list(NA_real_), term)), "finite")
    expect_error(predict(model, D="Name", levels=setNames(list(list(cov=2, wrong=1)), term)), "Unknown variables")
    Dt <- prediction$Dtable
    Dt$include[Dt$term == term] <- FALSE
    Dt$average[Dt$term == term] <- FALSE
    expect_warning(ignored <- predict(model, D="Name", Dtable=Dt), "ignored")
    expect_equal(as.numeric(ignored$D[,columns]), rep(0, nrow(ignored$D)*length(columns)))
  }
})

test_that("numeric reference grids respect fitted contrasts and included factors", {
  data <- prediction_issue_data()
  contrasts <- list(Env="contr.sum")
  fixed <- Yield~Env*cov
  model <- mmes(fixed, random=~Name, data=data, contrasts=contrasts, verbose=FALSE)
  prediction <- predict(model, D="Name")
  numericTerms <- names(model$partitionsX)[grepl("cov", names(model$partitionsX))]
  columns <- unlist(model$partitionsX[numericTerms], use.names=FALSE)
  grid <- data[seq_len(nlevels(data$Env)),]
  grid$Env <- factor(levels(data$Env), levels=levels(data$Env))
  grid$cov <- mean(data$cov)
  expected <- colMeans(model.matrix(fixed, grid, contrasts.arg=contrasts))
  expect_equal(unname(as.matrix(prediction$D[,columns])),
               matrix(rep(expected[columns], each=nrow(prediction$D)), nrow(prediction$D)), tolerance=1e-12)
  Dt <- model$Dtable
  Dt$include[Dt$type == "fixed"] <- TRUE
  included <- predict(model, D="Env", Dtable=Dt)
  expect_equal(unname(as.matrix(included$D[,columns])),
               unname(model.matrix(fixed, grid, contrasts.arg=contrasts)[,columns]), tolerance=1e-12)
  data$second <- data$cov/2
  product <- mmes(Yield~Env+cov:second, random=~Name, data=data, verbose=FALSE)
  term <- names(product$partitionsX)[grepl("cov", names(product$partitionsX))]
  result <- predict(product, D="Name", levels=setNames(list(list(cov=-2, second=3)), term))
  expect_equal(as.numeric(result$D[,as.vector(product$partitionsX[[term]])]), rep(-6, nrow(result$D)))
})

test_that("levels documentation example is reproducible through the Dtable", {
  data(DT_yatesoats, package="enhancer", envir=environment())
  model <- mmes(Y~V+N+V:N, random=~B+B:MP, data=DT_yatesoats, verbose=FALSE)
  Dt <- model$Dtable
  Dt$average[Dt$term %in% c("1", "V", "V:N")] <- TRUE
  Dt$include[Dt$term %in% c("N", "V:N")] <- TRUE
  prediction <- predict(model, D="N", Dtable=Dt, levels=list(N=c("0.6", "0.2")))
  Dt$levels[[match("N", Dt$term)]] <- c("0.6", "0.2")
  expect_equal(prediction$pvals$N, c("0.6", "0.2"))
  expect_equal(prediction$D, predict(model, D="N", Dtable=Dt)$D)
})