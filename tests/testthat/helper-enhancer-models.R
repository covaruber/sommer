.enhancer_regression_models <- function(){
  env <- new.env(parent=environment())
  for(name in c("DT_yatesoats", "DT_example", "DT_legendre", "DT_cpdata")){
    utils::data(list=name, package="enhancer", envir=env)
  }
  example <- env$DT_example
  example <- example[order(example$Env),]
  wide <- example
  wide[c("Yield","Weight")] <- lapply(wide[c("Yield","Weight")],
    function(x) as.numeric(scale(x)))
  long <- stackTraits(wide, c("Yield","Weight"), keep=c("Name","Env","Block"))
  cp <- env$DT_cpdata
  cp$Rowf <- factor(cp$Row); cp$Colf <- factor(cp$Col)
  relationship <- .95*A.mat(env$GT_cpdata)+.05*diag(nrow(cp))
  Ai <- solve(relationship)
  attr(Ai,"inverse") <- TRUE
  list(
    split_plot=list(fixed=Y~V+N+V:N, random=~B+B:MP, rcov=~units,
      data=env$DT_yatesoats),
    gxe_cs=list(fixed=Yield~Env, random=~Name+Env:Name, rcov=~units, data=example),
    gxe_diagonal=list(fixed=Yield~Env, random=~Name+vsm(dsm(Env),ism(Name)),
      rcov=~vsm(dsm(Env),ism(units)), data=example),
    multivariate=list(fixed=value~trait, random=~vsm(usm(trait),ism(Name)),
      rcov=~vsm(dsm(trait),ism(record)), data=long),
    random_regression=list(fixed=Y~1+Xf,
      random=~vsm(usm(leg(X,1)),ism(SUBJECT)), rcov=~units, data=env$DT_legendre),
    field_row_column=list(fixed=Yield~1, random=~id+Rowf+Colf,
      rcov=~units, data=cp),
    field_genomic=list(fixed=Yield~1, random=~vsm(ism(id),Gu=Ai)+Rowf+Colf,
      rcov=~units, data=cp)
  )
}
