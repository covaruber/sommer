## ----setup, message=FALSE-----------------------------------------------------
library(sommer)

## ----basic-call, eval=FALSE---------------------------------------------------
# fit <- mmes(fixed = outcome ~ treatment,
#   random = ~ subject,
#   rcov = ~ units,
#   data = dat,
#   family = binomial(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----structured-random, eval=FALSE--------------------------------------------
# fit <- mmes(y ~ environment,
#   random = ~ vsm(usm(environment), ism(genotype)),
#   rcov = ~ units,
#   data = dat,
#   family = poisson(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----offset, eval=FALSE-------------------------------------------------------
# fit <- mmes(events ~ treatment + offset(log(exposure)),
#   random = ~ site,
#   rcov = ~ units,
#   data = dat,
#   family = poisson(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----family-data, eval=FALSE--------------------------------------------------
# set.seed(2026)
# n_group <- 30
# n_per_group <- 8
# dat <- data.frame(
#   group = factor(rep(seq_len(n_group), each = n_per_group)),
#   x = rep(c(0, 1), length.out = n_group * n_per_group),
#   exposure = runif(n_group * n_per_group, 0.5, 2)
# )

## ----gaussian-example, eval=FALSE---------------------------------------------
# dat$y_gaussian <- 2 + 0.7 * dat$x + rnorm(nrow(dat), sd = 1)
# fit_gaussian <- mmes(
#   y_gaussian ~ x, random = ~ group, rcov = ~ units, data = dat,
#   family = gaussian()
# )

## ----binomial-example, eval=FALSE---------------------------------------------
# probability <- plogis(-0.7 + 1.1 * dat$x)
# dat$y_binomial <- rbinom(nrow(dat), size = 1, prob = probability)
# fit_binomial <- mmes(y_binomial ~ x,
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = binomial(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----poisson-example, eval=FALSE----------------------------------------------
# rate <- dat$exposure * exp(0.2 + 0.4 * dat$x)
# dat$y_poisson <- rpois(nrow(dat), lambda = rate)
# fit_poisson <- mmes(y_poisson ~ x + offset(log(exposure)),
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = poisson(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----gamma-example, eval=FALSE------------------------------------------------
# mean_gamma <- exp(0.3 + 0.25 * dat$x)
# dat$y_gamma <- rgamma(nrow(dat), shape = 4, scale = mean_gamma / 4)
# fit_gamma <- mmes(y_gamma ~ x,
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = Gamma(link = "log"),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----inverse-gaussian-example, eval=FALSE-------------------------------------
# dat$y_inverse_gaussian <- exp(0.2 + 0.3 * dat$x + rnorm(nrow(dat), sd = 0.2))
# fit_inverse_gaussian <- mmes(y_inverse_gaussian ~ x,
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = inverse.gaussian(link = "log"),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----quasi-example, eval=FALSE------------------------------------------------
# dat$y_quasi <- 1 + 0.5 * dat$x + rnorm(nrow(dat), sd = 0.5)
# fit_quasi <- mmes(y_quasi ~ x,
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = quasi(link = "identity", variance = "constant"),
#   control = mmes_control(update = list(blendWeight = 0)))
# 
# fit_quasibinomial <- mmes(y_binomial ~ x,
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = quasibinomial(),
#   control = mmes_control(update = list(blendWeight = 0)))
# 
# fit_quasipoisson <- mmes(y_poisson ~ x + offset(log(exposure)),
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = quasipoisson(),
#   control = mmes_control(update = list(blendWeight = 0)))

## ----extract, eval=FALSE------------------------------------------------------
# fitted(fit_poisson)
# fitted(fit_poisson, type = "link")
# residuals(fit_poisson, type = "deviance")
# residuals(fit_poisson, type = "working")
# 
# fit_poisson$pqlMonitor
# fit_poisson$pqlConverged
# fit_poisson$family

## ----pql-control, eval=FALSE--------------------------------------------------
# fit <- mmes(y_poisson ~ x + offset(log(exposure)),
#   random = ~ group,
#   rcov = ~ units,
#   data = dat,
#   family = poisson(),
#   control = mmes_control(
#     update = list(blendWeight = 0),
#     task = list(pql = list(maxit = 30, tol = 1e-6)),
#     convergence = list(maxIters = 20),
#     solve = list(method = "auto")
#   ))

