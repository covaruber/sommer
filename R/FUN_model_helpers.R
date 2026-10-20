# ---- Helpers ---------------------------------------------------------
.mmes_formula_env <- function(f, fallback){
  e <- if(inherits(f, "formula")) environment(f) else NULL
  if(is.null(e)) fallback else e
}

# Split only top-level additions. '+' inside I(), vsm(), etc. is untouched.
.mmes_split_plus <- function(expr){
  if(is.call(expr) && identical(expr[[1L]], as.name("+"))){
    c(.mmes_split_plus(expr[[2L]]), .mmes_split_plus(expr[[3L]]))
  } else list(expr)
}

.mmes_expr_label <- function(x) paste(deparse(x, width.cutoff=500L), collapse="")

.mmes_call_name <- function(expr){
  if(!is.call(expr)) return(NULL)
  head <- expr[[1L]]
  if(is.symbol(head)) return(as.character(head))
  if(is.call(head) && length(head) == 3L &&
     as.character(head[[1L]]) %in% c("::", ":::")){
    return(as.character(head[[3L]]))
  }
  NULL
}

.mmes_has_call <- function(expr, names){
  if(!is.call(expr)) return(FALSE)
  head <- .mmes_call_name(expr)
  if(length(head) == 1L && head %in% names) return(TRUE)
  any(vapply(as.list(expr)[-1L], .mmes_has_call, logical(1), names=names))
}

.mmes_eval_model_expr <- function(expr, data_full, enclos){
  if(is.null(data_full)) eval(expr, envir=enclos)
  else eval(expr, envir=data_full, enclos=enclos)
}

# Observation-level variables are symbols resolving to vectors of length n,
# or matrices/data.frames with n rows. Objects such as Gu/Ai are therefore
# not mistaken for observation variables unless they actually have n rows.
.mmes_observation_ok <- function(expr, data_full, enclos, n){
  vars <- unique(all.vars(expr))
  if(!length(vars)) return(rep(TRUE, n))
  ok <- rep(TRUE, n)
  for(v in vars){
    val <- tryCatch({
      if(!is.null(data_full) && v %in% names(data_full)) data_full[[v]]
      else get(v, envir=enclos, inherits=TRUE)
    }, error=function(e) NULL)
    if(is.null(val)) next
    if(is.data.frame(val) || is.matrix(val)){
      if(nrow(val) == n) ok <- ok & stats::complete.cases(val)
    } else if(length(val) == n){
      ok <- ok & !is.na(val)
    }
  }
  ok
}

.mmes_method_keep <- function(ok, method, what){
  method <- tolower(method)
  if(method %in% c("exclude", "omit")) return(ok)
  if(method %in% c("include", "pass")) return(rep(TRUE, length(ok)))
  if(method == "fail" && any(!ok))
    stop("Missing values found in ", what, ".", call.=FALSE)
  if(method == "fail") return(rep(TRUE, length(ok)))
  stop("Unknown missing-data method '", method, "' for ", what, ".", call.=FALSE)
}
