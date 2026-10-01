test_that("sparse non-diagonal weights preserve diagonal-residual fits", {
  set.seed(7)
  n <- 60L
  data <- data.frame(id=factor(rep(seq_len(n / 2L), each=2L)), y=rnorm(n))
  weights <- Matrix::bdiag(replicate(n / 2L,
                                    matrix(c(1.4, -0.2, -0.2, 1.4), 2L),
                                    simplify=FALSE))
  sparse <- suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data,
                                  W=weights, nIters=3, solver="ldlt",
                                  verbose=FALSE, dateWarning=FALSE))
  expect_equal(unname(unlist(sparse$theta)),
               c(0.1532013323059318, 1.0711724582035496), tolerance=1e-8)
  expect_equal(head(as.numeric(sparse$bu), 5L),
               c(0.21835699494394373, 0.09423111073305052,
                 -0.19829464766940430, -0.30130806010703842,
                 0.02750949952965763), tolerance=1e-8)
})