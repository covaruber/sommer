.mmes_glmm_path_plan <- function(family, method, REML, rcov=NULL){
  families <- if(inherits(family, "sommer_familym")) family$families else list(family)
  codes <- vapply(families, function(component) component$family, character(1))
  independent <- is.null(rcov) ||
    identical(attr(stats::terms(rcov), "term.labels"), "units")
  gaussianBlocks <- any(codes == "gaussian")
  plan <- sommer_glmm_path_plan(method, !isTRUE(REML),
    any(codes == "ordinal"), any(codes == "multinomial"),
    any(codes == "ziPoisson"), independent,gaussianBlocks)
  plan$mixed <- inherits(family, "sommer_familym")
  plan
}

.mmes_require_glmm_path <- function(plan){
  if(!isTRUE(plan$eligible)){
    stop(structure(list(message=plan$message, call=NULL, glmmPathPlan=plan),
      class=c("sommer_glmm_path_error", "error", "condition")))
  }
  invisible(plan)
}
