# Run from the package root after build-mmes-options.R.
stopifnot(requireNamespace("jsonlite", quietly=TRUE))
page <- paste(readLines("pkgdown/assets/mmes-options.html", warn=FALSE), collapse="\n")
json <- regmatches(page, regexpr('(?s)(?<=<script id="catalogData" type="application/json">).*?(?=</script>)',
  page, perl=TRUE))
stopifnot(length(json) == 1L, nzchar(json))
data <- jsonlite::fromJSON(json, simplifyVector=FALSE)
src <- new.env(parent=baseenv())
for(file in c("R/FUN_mme_control.R", "R/FUN_mme_capabilities.R",
    "R/FUN_mme_modules.R", "R/FUN_nonlinear_paths.R", "R/FUN_model_interface.R")){
  sys.source(file, envir=src)
}
paths <- vapply(data$controls, `[[`, "", "path")
ids <- vapply(data$nodes, `[[`, "", "id")
placed <- unlist(lapply(data$nodes, `[[`, "controls"), use.names=FALSE)
stopifnot(!anyDuplicated(paths), !anyDuplicated(ids), setequal(paths, placed),
  !"control.output.pev" %in% paths,
  !any(startsWith(paths,"control.task.")),
  setequal(paths,src$.mmes_public_catalog()$controls$path))
for(edge in data$edges) stopifnot(edge$from %in% ids, edge$to %in% ids)
rule_ids <- vapply(data$rules, `[[`, "", "id")
source_rules <- vapply(src$.mmes_compatibility_rules(), `[[`, "", "id")
stopifnot(!anyDuplicated(rule_ids), all(source_rules %in% rule_ids))
for(rule in data$rules) stopifnot(all(unlist(rule$nodes) %in% ids))
for(node in data$nodes) stopifnot(all(unlist(node$rules) %in% rule_ids))
for(row in data$controls){
  stopifnot(nzchar(row$values), nzchar(row$default), nzchar(row$effect), nzchar(row$usage))
  if(startsWith(row$path, "control.") || startsWith(row$path,"method.")){
    is_control <- startsWith(row$path,"control.")
    defaults <- if(is_control) src$.mmes_public_control_defaults() else
      unclass(src$mmes_method())
    parts <- strsplit(sub("^[^.]+\\.", "", row$path), ".", fixed=TRUE)[[1]]
    value <- Reduce(function(x, part) x[[part]], parts, defaults)
    stopifnot(identical(row$default,paste(deparse(value),collapse="")))
    example <- gsub("...", paste(deparse(value), collapse=""), row$usage, fixed=TRUE)
    control <- eval(parse(text=example), envir=src)
    stopifnot(inherits(control, if(is_control) "sommer_mmes_control" else
      "sommer_mmes_method"))
  }
}
check_constructors <- function(expr){
  if(!is.call(expr)) return(invisible(NULL))
  if(is.symbol(expr[[1L]]) &&
     as.character(expr[[1L]]) %in% c("mmes_control","mmes_method")){
    value <- eval(expr,envir=src)
    stopifnot(inherits(value,c("sommer_mmes_control","sommer_mmes_method")))
    return(invisible(NULL))
  }
  for(child in as.list(expr)[-1L]) check_constructors(child)
}
for(node in data$nodes){
  stopifnot(length(node$controls) > 0L, all(unlist(node$controls) %in% paths))
  for(expr in parse(text=node$example)) check_constructors(expr)
}
for(row in data$structures){
  if(grepl("^mmes_control\\(", row$example)){
    control <- eval(parse(text=row$example), envir=src)
    stopifnot(inherits(control, "sommer_mmes_control"))
  }
}
defaults <- src$mmes_control()
conditional <- src$mmes_control(representation=list(equations="conditional"))
stopifnot(identical(defaults$output$inverse, 0L),
  identical(defaults$formulation, "henderson"),
  identical(defaults$solve$method, "auto"),
  identical(conditional$representation$equations, "conditional"))
method <- src$mmes_method()
resolved <- src$.mmes_resolve_controls(defaults,method,FALSE,NULL,list())
stopifnot(identical(method$approximation,"auto"),
  identical(attr(defaults,"publicArguments")$optimizer$gradient,"auto"),
  identical(attr(defaults,"publicArguments")$mode$warmStart,"auto"),
  identical(names(attr(defaults,"publicArguments")$likelihood),
    c("backend","trace","blocks")),
  identical(resolved$task$glmm$warmStart,"auto"),
  identical(resolved$task$glmm$method,"pql"),
  identical(resolved$task$objective,"reml"),
  !"task" %in% names(attr(defaults,"publicArguments")))
rejects <- function(expr) inherits(tryCatch(force(expr),error=identity),"error")
stopifnot(
  rejects(src$mmes_control(task=list(objective="ml"))),
  rejects(src$mmes_control(likelihood=list(warmStart="off"))),
  rejects(src$.mmes_resolve_controls(src$mmes_control(mode=list(warmStart="on")),
    src$mmes_method(objective="ml"),TRUE,NULL,list())),
  rejects(src$.mmes_resolve_controls(defaults,
    src$mmes_method(approximation="pql",integration="laplace"),FALSE,NULL,list())),
  rejects(src$.mmes_resolve_controls(defaults,
    src$mmes_method(approximation="pql"),TRUE,NULL,list())))
nonlinear <- src$.mmes_resolve_controls(defaults,
  src$mmes_method(objective="ml"),TRUE,NULL,list())
stopifnot(identical(nonlinear$task$nonlinear$mode$warmStart,"off"))
cat(sprintf("Guide checks passed: %d options, %d boxes, all constructor examples valid.\n",
  length(paths), length(ids)))
