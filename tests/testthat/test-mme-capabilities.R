test_that("capability ownership covers every modular leaf and model formal", {
  map <- sommer:::.mmes_public_catalog()
  leafPaths <- function(value, prefix){
    unlist(lapply(names(value), function(name){
      path <- paste0(prefix, ".", name)
      if(is.list(value[[name]]) && length(value[[name]]))
        leafPaths(value[[name]], path) else path
    }), use.names=FALSE)
  }
  control <- sommer:::.mmes_public_control_defaults(mmes_control())
  paths <- "control.formulation"
  for(module in setdiff(names(control), "formulation")){
    paths <- c(paths, leafPaths(control[[module]], paste0("control.", module)))
  }
  expect_setequal(map$controls$path[startsWith(map$controls$path, "control.")], paths)
  publicModel <- setdiff(names(formals(mmes)), c("control", ".internal", "..."))
  expect_setequal(sub("^mmes\\.", "", map$controls$path[
    startsWith(map$controls$path, "mmes.")]), publicModel)
  expect_false(anyDuplicated(map$controls$path) > 0L)
  expect_true(all(map$controls$module %in% map$modules$id))
  expect_true(all(nzchar(map$controls$requires)))
  expect_true(all(nzchar(map$controls$effect)))
  expect_false(any(map$controls$status == "alias"))
  expect_false(any(startsWith(map$controls$path, "control.task.")))
  expect_setequal(map$controls$path[startsWith(map$controls$path, "method.")],
    paste0("method.", names(formals(mmes_method))))
  expect_true(all(is.na(map$map) | map$map %in% map$controls$path))
  gradient <- map$controls[map$controls$path=="control.optimizer.gradient", ]
  expect_identical(gradient$default, '"auto"')
  expect_identical(gradient$values, "auto, analytic, central.")
  for(value in c("auto","analytic","central"))
    expect_identical(attr(mmes_control(optimizer=list(gradient=value)),
      "publicArguments")$optimizer$gradient, value)
  expect_false("control.likelihood.warmStart" %in% map$controls$path)
  expect_identical(control$mode$warmStart,"auto")
  expect_identical(names(control$likelihood),c("backend","trace","blocks"))
  expect_identical(unname(map$map["control.task.glmm.warmStart"]),
    "control.mode.warmStart")
})
