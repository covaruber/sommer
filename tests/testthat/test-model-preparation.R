test_that("formula helpers preserve nested calls and formula environments", {
  expression <- quote(vsm(ar1m(env), ism(id)) + I(x + y))
  pieces <- sommer:::.mmes_split_plus(expression)
  expect_length(pieces, 2L)
  expect_identical(pieces[[2L]], quote(I(x + y)))
  expect_true(sommer:::.mmes_has_call(quote(sommer::vsm(sommer::ism(id))), "ism"))
  expect_identical(sommer:::.mmes_call_name(quote(sommer::vsm(x))), "vsm")
  environment <- new.env(parent=baseenv())
  environment$x <- c(1, NA, 3)
  formula <- stats::as.formula("~x", env=environment)
  expect_identical(sommer:::.mmes_formula_env(formula, globalenv()), environment)
  expect_identical(sommer:::.mmes_observation_ok(formula, NULL, environment, 3L),
    c(TRUE, FALSE, TRUE))
  data <- data.frame(x=c(4, 5, 6))
  expect_identical(sommer:::.mmes_eval_model_expr(quote(x), data, environment), data$x)
})

test_that("weight preparation preserves filtered and cached PQL precision", {
  prepare <- sommer:::.mmes_prepare_weights
  keep <- c(TRUE, FALSE, TRUE, TRUE)
  weights <- diag(c(2, 3, 4, 5))
  identity <- prepare(NULL, TRUE, 4L, keep, 3L, NULL)
  expect_false(identity$useH)
  expect_equal(as.matrix(identity$W), diag(3))
  filtered <- prepare(weights, FALSE, 4L, keep, 3L, NULL)
  retained <- weights[keep, keep]
  expect_equal(as.matrix(filtered$W), retained)
  expect_equal(prepare(retained, FALSE, 4L, keep, 3L, NULL), filtered)

  precision <- c(1, 2, 3, 4)
  ordinary <- prepare(NULL, TRUE, 4L, keep, 3L, NULL,
    pqlPrecision=precision, pqlBaseW=weights)
  cached <- prepare(NULL, TRUE, 4L, keep, 3L, NULL,
    pqlPrecision=precision[keep], pqlBaseFactor=Matrix::chol(Matrix::Matrix(retained)))
  expect_true(ordinary$useH)
  expect_equal(as.matrix(ordinary$W), retained %*% diag(precision[keep]))
  expect_equal(cached, ordinary)
  expect_error(prepare(weights, FALSE, 4L, keep, 3L, c(1, 2)), "retained W")
  expect_error(prepare(NULL, TRUE, 4L, keep, 3L, c(1, 1, 1)), "provide W")
  expect_error(prepare(NULL, TRUE, 4L, keep, 3L, NULL,
    pqlPrecision=c(1, -1, 1)), "finite and positive")
})
