# Run from the package root after build-mmes-options.R.
stopifnot(requireNamespace("jsonlite", quietly=TRUE))
page <- paste(readLines("pkgdown/assets/mmes-options.html", warn=FALSE), collapse="\n")
json <- regmatches(page, regexpr('(?s)(?<=<script id="catalogData" type="application/json">).*?(?=</script>)',
  page, perl=TRUE))
stopifnot(length(json) == 1L, nzchar(json))
data <- jsonlite::fromJSON(json, simplifyVector=FALSE)
src <- new.env(parent=baseenv())
for(file in c("R/FUN_mme_control.R", "R/FUN_mme_capabilities.R", "R/FUN_mme_modules.R")){
  sys.source(file, envir=src)
}
paths <- vapply(data$controls, `[[`, "", "path")
ids <- vapply(data$nodes, `[[`, "", "id")
placed <- unlist(lapply(data$nodes, `[[`, "controls"), use.names=FALSE)
stopifnot(!anyDuplicated(paths), !anyDuplicated(ids), setequal(paths, placed),
  !"control.output.pev" %in% paths)
for(edge in data$edges) stopifnot(edge$from %in% ids, edge$to %in% ids)
for(row in data$controls){
  stopifnot(nzchar(row$values), nzchar(row$default), nzchar(row$effect), nzchar(row$usage))
}
for(node in data$nodes){
  stopifnot(length(node$controls) > 0L, all(unlist(node$controls) %in% paths))
  if(grepl("^mmes_control\\(", node$example)){
    control <- eval(parse(text=node$example), envir=src)
    stopifnot(inherits(control, "sommer_mmes_control"))
  }
}
for(row in data$structures){
  if(grepl("^mmes_control\\(", row$example)){
    control <- eval(parse(text=row$example), envir=src)
    stopifnot(inherits(control, "sommer_mmes_control"))
  }
}
defaults <- src$mmes_control()
stopifnot(identical(defaults$output$inverse, 0L),
  identical(defaults$formulation, "henderson"),
  identical(defaults$solve$method, "auto"))
cat(sprintf("Guide checks passed: %d options, %d boxes, all constructor examples valid.\n",
  length(paths), length(ids)))
