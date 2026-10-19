## ----eval=FALSE---------------------------------------------------------------
# vsm(
#   dsm(Env),
#   ar1m(Row),
#   ism(Name)
# )

## ----eval=FALSE---------------------------------------------------------------
# library(sommer)
# data(DT_example, package="enhancer")
# DT <- DT_example

## ----eval=FALSE---------------------------------------------------------------
# fit_cs <- mmes(Yield ~ Env,
#   random = ~ vsm(csm(Env), ism(Name)),
#   rcov = ~ vsm(ism(units)),
#   data = DT,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))

## ----eval=FALSE---------------------------------------------------------------
# n_env <- nlevels(factor(DT$Env))
# 
# fit_csh <- mmes(Yield ~ Env,
#   random = ~ vsm(
#     csm(Env, variance="heterogeneous", values=rep(1, n_env)),
#     ism(Name)
#   ),
#   rcov = ~ vsm(ism(units)),
#   data = DT,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))

## ----eval=FALSE---------------------------------------------------------------
# DT$EnvOrder <- factor(DT$Env, levels=unique(DT$Env), ordered=TRUE)
# n_env <- nlevels(DT$EnvOrder)
# 
# fit_ar1h <- mmes(Yield ~ Env,
#   random = ~ vsm(
#     ar1m(EnvOrder, rho=0.30, variance="heterogeneous",
#          values=rep(1, n_env)),
#     ism(Name)
#   ),
#   rcov = ~ vsm(ism(units)),
#   data = DT,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))

## ----eval=FALSE---------------------------------------------------------------
# n_env <- nlevels(factor(DT$Env))
# 
# environment_cs <- csm(
#   DT$Env,
#   variance="heterogeneous",
#   values=rep(1, n_env),
#   fixed=c(FALSE, rep(TRUE, n_env - 1L))
# )

## ----eval=FALSE---------------------------------------------------------------
# fit_fixed_residual <- mmes(Yield ~ Env,
#   random = ~ vsm(ism(Name)),
#   rcov = ~ vsm(ism(units), sigma2=1, fixedSigma2=TRUE),
#   data = DT,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))

## ----eval=FALSE---------------------------------------------------------------
# library(sommer)
# data(DT_example, package="enhancer")
# DT <- DT_example
# 
# env_levels <- unique(as.character(DT$Env))
# DT$EnvOrder <- factor(DT$Env, levels=env_levels, ordered=TRUE)
# DT$EnvCoordinate <- as.numeric(DT$EnvOrder)
# n_env <- length(env_levels)
# 
# # DT_example has three environments. AR(3) needs at least four ordered
# # levels, so this demonstration-only partition is used for the AR(3) rows
# # below. Replace it with a scientific time, distance, or ordered factor.
# DT$CatalogOrder <- factor(rep(seq_len(4L), length.out=nrow(DT)), ordered=TRUE)
# n_catalog <- nlevels(DT$CatalogOrder)
# 
# # Named first-neighbour adjacency among ordered environments.
# W_env <- matrix(0, n_env, n_env,
#                 dimnames=list(env_levels, env_levels))
# W_env[cbind(seq_len(n_env - 1L), 2:n_env)] <- 1
# W_env[cbind(2:n_env, seq_len(n_env - 1L))] <- 1
# 
# # A known positive-definite covariance shape for ownm().
# K_env <- 0.40 ^ abs(outer(seq_len(n_env), seq_len(n_env), "-"))
# 
# # Every object below has the same observation layout and can be used as
# # a covariance factor in vsm(shape, ism(Name)). Most use EnvOrder; AR(3)
# # uses CatalogOrder because EnvOrder has too few levels for that structure.
# environment_shapes <- list(
#   identity = ism(DT$EnvOrder),
#   diagonal = dsm(DT$EnvOrder),
#   selected_diagonal = atm(DT$EnvOrder, levs=env_levels[1:3]),
#   compound_symmetry = csm(DT$EnvOrder),
#   compound_symmetry_heterogeneous = csm(
#     DT$EnvOrder, variance="heterogeneous", values=rep(1, n_env)
#   ),
#   ar1 = ar1m(DT$EnvOrder),
#   ar1_heterogeneous = ar1m(
#     DT$EnvOrder, variance="heterogeneous", values=rep(1, n_env)
#   ),
#   ar2 = ar2m(DT$EnvOrder),
#   ar2_heterogeneous = ar2m(
#     DT$EnvOrder, variance="heterogeneous", values=rep(1, n_env)
#   ),
#   ar3 = ar3m(DT$CatalogOrder),
#   ar3_heterogeneous = ar3m(
#     DT$CatalogOrder, variance="heterogeneous", values=rep(1, n_catalog)
#   ),
#   ma1 = mam(DT$EnvOrder, order=1L),
#   ma2 = mam(DT$EnvOrder, order=2L),
#   unstructured = usm(DT$EnvOrder),
#   general_correlation = corgm(DT$EnvOrder),
#   factor_analytic = fam(DT$EnvOrder, k=1L),
#   antedependence = antem(DT$EnvOrder, order=1L),
#   user_defined = ownm(DT$EnvOrder, K=K_env),
#   reduced_rank = rrm(DT$EnvOrder, k=1L),
#   matern = maternm(DT$EnvCoordinate),
#   toeplitz = toeplitzm(DT$EnvOrder),
#   sar = sar(DT$EnvOrder, W=W_env),
#   car = car(DT$EnvOrder, W=W_env)
# )

## ----eval=FALSE---------------------------------------------------------------
# fit_environment_shape <- function(shape){
#   mmes(Yield ~ Env,
#   random = ~ vsm(shape, ism(DT$Name)),
#   rcov = ~ vsm(ism(units)),
#   data = DT,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))
# }
# 
# fit_identity <- fit_environment_shape(environment_shapes$identity)
# fit_ar1 <- fit_environment_shape(environment_shapes$ar1)
# fit_matern <- fit_environment_shape(environment_shapes$matern)
# fit_car <- fit_environment_shape(environment_shapes$car)

## ----eval=FALSE---------------------------------------------------------------
# effect_1 <- vsm(ism(DT$Name))
# effect_2 <- vsm(ism(DT$Name))
# joint_effect <- covm(effect_1, effect_2, labels=c("effect_1", "effect_2"))

## ----eval=FALSE---------------------------------------------------------------
# fit <- mmes(y ~ 1,
#             random = ~ strm(dir = vsm(ism(id)), mat = vsm(ism(dam)),
#                             pe = vsm(ism(pe)), cov = usm, Gu = Ainv),
#             data = animals)
# covparams_mmes(fit, 1)  # variances and covariances among dir, mat and pe

## ----eval=FALSE---------------------------------------------------------------
# p <- mmes(Y ~ V * N,
#   random = ~ B + B:MP,
#   rcov = ~ units,
#   data = DT_yatesoats,
#   control = mmes_control(output = list(setup = TRUE)))
# p$vcParams
# fit <- mmes(Y ~ V * N, random = ~ B + B:MP, rcov = ~ units, data = DT_yatesoats,
#             vcc = data.frame(parameter = c("vsm(ism(B:MP)):sigma2", "vsm(ism(B)):sigma2"),
#                              group = 1, scale = c(1, 2)))

## ----eval=FALSE---------------------------------------------------------------
# structure <- csm(DT$Env, variance="heterogeneous")
# str(structure$covFactor)

## ----eval=FALSE---------------------------------------------------------------
# random_structure <- vsm(csm(DT$Env, variance="heterogeneous"), ism(DT$Name))
# random_structure$covStruct$par_names
# random_structure$covStruct$free

