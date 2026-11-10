test_that("GLMM decision policy distinguishes likelihood layout and estimator", {
  scalar <- familym("trait",count=poisson(),continuous=gaussian())
  ordinalMixed <- familym("trait",rating=ordinal(),count=poisson())
  nominalMixed <- familym("trait",choice=multinomial(),count=poisson())
  both <- familym("trait",rating=ordinal(),choice=multinomial(),count=poisson())
  for(family in list(scalar,ordinalMixed,nominalMixed,both)){
    nominal <- any(vapply(family$families,function(x)x$family=="multinomial",logical(1)))
    plan <- sommer:::.mmes_glmm_path_plan(family,"laplace",FALSE,~units)
    expect_true(plan$eligible)
    expect_identical(plan$requested,plan$selected)
    expect_identical(plan$selected,if(nominal)
      "native-laplace-category-block" else "native-laplace-scalar")
    expect_identical(plan$likelihood,if(nominal) "category-block" else "scalar")
    expect_identical(plan$reason,"eligible")
    expect_true(plan$mixed)
  }
  for(REML in c(TRUE,FALSE)){
    # PQL retains its working-model objective and structured residual support.
    plan <- sommer:::.mmes_glmm_path_plan(scalar,"pql",REML,~vsm(ism(group)))
    expect_true(plan$eligible)
    expect_identical(plan$selected,"pql")
  }
  expect_identical(sommer:::.mmes_glmm_path_plan(ordinal(),"laplace",FALSE)$composition,
    "ordinal")
  expect_identical(sommer:::.mmes_glmm_path_plan(multinomial(),"laplace",FALSE)$selected,
    "native-laplace-category-block")
})

test_that("GLMM policy reports stable reasons without silently falling back", {
  family <- familym("trait",rating=ordinal(),choice=multinomial(),count=poisson())
  cases <- list(
    list(method="pql",REML=FALSE,rcov=NULL,reason="categorical-pql-unsupported"),
    list(method="laplace",REML=TRUE,rcov=NULL,reason="laplace-requires-ml"),
    list(method="laplace",REML=FALSE,rcov=~vsm(ism(group)),
      reason="laplace-residual-dependence"))
  for(case in cases){
    plan <- sommer:::.mmes_glmm_path_plan(family,case$method,case$REML,case$rcov)
    expect_false(plan$eligible)
    expect_identical(plan$selected,"none")
    expect_identical(plan$reason,case$reason)
    error <- tryCatch(sommer:::.mmes_require_glmm_path(plan),error=identity)
    expect_s3_class(error,"sommer_glmm_path_error")
    expect_identical(error$glmmPathPlan,plan)
  }
  zip <- sommer:::.mmes_glmm_path_plan(zipoisson(),"pql",FALSE)
  expect_identical(zip$reason,"mixture-pql-unsupported")
})

test_that("public dispatch rejects ineligible paths before invoking a fitter", {
  local_mocked_bindings(
    .mmes_laplace=function(...)stop("Laplace backend must not run"),
    .mmes_pql=function(...)stop("PQL backend must not run"),.package="sommer")
  data <- data.frame(y=rep(1:3,8),trait=rep(c("rating","count"),each=12),
    group=factor(rep(1:4,6)))
  family <- familym("trait",rating=ordinal(),count=poisson())
  for(case in list(
      list(method="pql",objective="ml",rcov=~units,reason="categorical-pql-unsupported"),
      list(method="laplace",objective="reml",rcov=~units,reason="laplace-requires-ml"),
      list(method="laplace",objective="ml",rcov=~vsm(ism(group)),
        reason="laplace-residual-dependence"))){
    control <- sommer:::.mmes_control_internal(task=list(objective=case$objective,
      glmm=list(method=case$method)),diagnostics=list(verbose=FALSE,dateWarning=FALSE))
    error <- tryCatch(sommer:::.mmes_dispatch(y~trait,random=~group,rcov=case$rcov,data=data,
      family=family,control=control),error=identity)
    expect_s3_class(error,"sommer_glmm_path_error")
    expect_identical(error$glmmPathPlan$reason,case$reason)
    expect_identical(error$glmmPathPlan$selected,"none")
  }
})
