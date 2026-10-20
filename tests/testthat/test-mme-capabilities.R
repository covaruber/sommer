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
  publicModel <- setdiff(names(formals(mmes)), c("control", ".internal", "..."))
  expect_setequal(sub("^mmes\\.", "", map$controls$path[
    startsWith(map$controls$path, "mmes.")]), publicModel)
  expect_false(anyDuplicated(map$controls$path) > 0L)
  expect_true(all(map$controls$module %in% map$modules$id))
  expect_true(all(nzchar(map$controls$requires)))
  expect_true(all(nzchar(map$controls$effect)))
  expect_false(any(map$controls$status == "alias"))
})
