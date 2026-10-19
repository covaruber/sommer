test_that("capability ownership covers every modular leaf and model formal", {
  map <- sommer:::.mmes_capabilities()
  leafPaths <- function(value, prefix){
    unlist(lapply(names(value), function(name){
      path <- paste0(prefix, ".", name)
      if(is.list(value[[name]]) && name == "pql") leafPaths(value[[name]], path) else path
    }), use.names=FALSE)
  }
  control <- mmes_control()
  paths <- "control.formulation"
  for(module in setdiff(names(control), "formulation")){
    paths <- c(paths, leafPaths(control[[module]], paste0("control.", module)))
  }
  expect_setequal(map$controls$path[startsWith(map$controls$path, "control.")], paths)
  publicModel <- setdiff(names(formals(mmes)), c("control", ".internal"))
  expect_setequal(sub("^mmes\\.", "", map$controls$path[
    startsWith(map$controls$path, "mmes.")]), publicModel)
  expect_false(anyDuplicated(map$controls$path) > 0L)
  expect_true(all(map$controls$module %in% map$modules$id))
  expect_true(all(nzchar(map$controls$requires)))
  expect_true(all(nzchar(map$controls$effect)))
  expect_false(any(map$controls$status == "alias"))
})

test_that("decision graph remains connected after API replacement", {
  map <- sommer:::.mmes_capabilities()
  nodes <- map$decisions$id
  expect_true(all(map$edges$from %in% nodes))
  expect_true(all(map$edges$to %in% nodes))
  reachable <- "task"
  repeat {
    nextNodes <- union(reachable, map$edges$to[map$edges$from %in% reachable])
    if(identical(reachable, nextNodes)) break
    reachable <- nextNodes
  }
  expect_setequal(reachable, nodes)
})

test_that("capability inspection is independent of global option resolution", {
  before <- options()
  map <- sommer:::.mmes_capabilities()
  expect_identical(options(), before)
  old <- options(sommer.mme.gpu=TRUE, sommer.mme.commonScale="force")
  on.exit(options(old), add=TRUE)
  expect_identical(sommer:::.mmes_capabilities(), map)
})
