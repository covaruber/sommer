#### =========== ####
## VARIANCE-PARAMETER TABLE AND EQUALITY/SCALING CONSTRAINTS (ASReml vcc)
#### =========== ####

.vc_kind <- function(cs, j){
  if(j == 1L) return(list(kind="scale", lower=NA_real_, upper=NA_real_))
  for(f in cs$factors){
    s <- as.integer(f$par_start); e <- as.integer(f$par_end)
    if(e >= s && j >= s && j <= e){
      k <- j - s + 1L
      rep <- f$report
      if(is.null(rep) || !identical(rep$backend, "builtin")){
        return(list(kind="identity", lower=NA_real_, upper=NA_real_))
      }
      return(list(kind=rep$transform[k], lower=rep$lower[k], upper=rep$upper[k]))
    }
  }
  list(kind="identity", lower=NA_real_, upper=NA_real_)
}

.vc_to_natural <- function(w, kind, lower, upper){
  switch(kind, scale=exp(w), exp=exp(w), identity=w, tanh=tanh(w),
         bounded_logit=lower + (upper - lower) * stats::plogis(w), w)
}

.vc_to_working <- function(v, kind, lower, upper){
  switch(kind, scale=log(v), exp=log(v), identity=v, tanh=atanh(v),
         bounded_logit=stats::qlogis((v - lower) / (upper - lower)), v)
}

.vc_param_table <- function(covStruct, termLabels){
  rows <- lapply(seq_along(covStruct), function(i){
    cs <- covStruct[[i]]
    kinds <- lapply(seq_along(cs$par), function(j) .vc_kind(cs, j))
    data.frame(term=termLabels[i], parameter=as.character(cs$par_names),
               kind=vapply(kinds, `[[`, character(1), "kind"),
               start=vapply(seq_along(cs$par), function(j)
                 .vc_to_natural(cs$par[j], kinds[[j]]$kind, kinds[[j]]$lower, kinds[[j]]$upper),
                 numeric(1)),
               free=as.logical(cs$free), structure=i, position=seq_along(cs$par),
               stringsAsFactors=FALSE)
  })
  out <- do.call(rbind, rows)
  out <- cbind(index=seq_len(nrow(out)), out)
  rownames(out) <- NULL
  out
}

.vc_constraint_map <- function(covStruct, vcParams, vcc){
  vcc <- as.data.frame(vcc, stringsAsFactors=FALSE)
  if(!all(c("parameter", "group") %in% names(vcc))){
    stop("vcc must be a data frame with columns 'parameter' and 'group' (and optionally 'scale').",
         call.=FALSE)
  }
  if(is.null(vcc$scale)) vcc$scale <- 1
  fullNames <- paste(vcParams$term, vcParams$parameter, sep=":")
  idx <- if(is.numeric(vcc$parameter)) as.integer(vcc$parameter) else {
    p <- as.character(vcc$parameter)
    m <- match(p, fullNames)
    alone <- is.na(m) & p %in% vcParams$parameter &
      vapply(p, function(z) sum(vcParams$parameter == z) == 1L, logical(1))
    m[alone] <- match(p[alone], vcParams$parameter)
    m
  }
  if(anyNA(idx) || any(idx < 1L | idx > nrow(vcParams))){
    stop("Unknown parameters in vcc: ", paste(vcc$parameter[is.na(idx)], collapse=", "),
         ". See the vcParams table from mmes(..., returnParam=TRUE).", call.=FALSE)
  }
  if(anyDuplicated(idx)) stop("A parameter appears more than once in vcc.", call.=FALSE)
  if(any(!is.finite(vcc$scale)) || any(vcc$scale == 0)){
    stop("vcc scale values must be finite and non-zero.", call.=FALSE)
  }

  nPar <- nrow(vcParams)
  group <- rep(NA_integer_, nPar)
  scale <- rep(1, nPar)
  free <- vcParams$free
  natural <- vcParams$start
  groups <- split(seq_along(idx), vcc$group)
  for(g in names(groups)){
    rowsG <- groups[[g]]
    members <- idx[rowsG]
    if(length(members) < 2L) stop("Constraint group ", g, " has a single parameter.", call.=FALSE)
    kinds <- vcParams$kind[members]
    if(length(unique(kinds)) != 1L){
      stop("Constraint group ", g, " mixes parameters of different kinds (",
           paste(unique(kinds), collapse=", "), "); variance scales (sigma2) can only be ",
           "constrained with other variance scales, and other parameters with parameters of the same transform.",
           call.=FALSE)
    }
    kind <- kinds[1L]
    s <- vcc$scale[rowsG] / vcc$scale[rowsG[1L]]
    if(kind %in% c("scale", "exp") && any(s <= 0)){
      stop("Scaling coefficients of variance parameters must be positive.", call.=FALSE)
    }
    if(kind %in% c("tanh", "bounded_logit") && any(abs(s - 1) > 1e-12)){
      stop("Correlation-type parameters can only be constrained to be equal (scale = 1).",
           call.=FALSE)
    }
    group[members] <- as.integer(match(g, names(groups)))
    scale[members] <- s
    anchor <- members[1L]
    if(any(!free[members])){
      anchor <- members[!free[members]][1L]
      free[members] <- FALSE
    }
    natural[members] <- natural[anchor] / scale[anchor] * s
  }

  # working values consistent with the constraints, written back as starting values
  structures <- vcParams$structure
  positions <- vcParams$position
  for(k in which(!is.na(group))){
    cs <- covStruct[[structures[k]]]
    kd <- .vc_kind(cs, positions[k])
    covStruct[[structures[k]]]$par[positions[k]] <-
      .vc_to_working(natural[k], kd$kind, kd$lower, kd$upper)
    covStruct[[structures[k]]]$free[positions[k]] <- free[k]
    if(positions[k] == 1L) covStruct[[structures[k]]]$sigma2_is_default <- FALSE
  }

  # w = o + T phi over free parameters; fixed rows are zero and never projected
  cols <- list()
  o <- numeric(nPar)
  for(k in which(free & is.na(group))){
    v <- numeric(nPar); v[k] <- 1; cols[[length(cols) + 1L]] <- v
  }
  for(g in sort(unique(group[!is.na(group) & free]))){
    members <- which(group == g)
    v <- numeric(nPar)
    kind <- vcParams$kind[members[1L]]
    if(kind %in% c("scale", "exp")){
      v[members] <- 1
      o[members] <- log(scale[members])
    }else{
      v[members] <- scale[members]
    }
    cols[[length(cols) + 1L]] <- v
  }
  T <- if(length(cols)) do.call(cbind, cols) else matrix(0, nPar, 0L)
  list(covStruct=covStruct, T=T, o=o, group=group)
}
