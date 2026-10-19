## ----setup, include=FALSE-----------------------------------------------------
library(sommer)

## ----eval=FALSE---------------------------------------------------------------
# # Standalone scalar PX-GEM
# fit <- mmes(y ~ 1,
#   random = ~id,
#   data = dat,
#   control = mmes_control(update = list(recovery = "pxem", blendWeight = 1)))
# 
# # PX-GEM/AI with a five-iteration exponential taper
# fit <- mmes(y ~ 1,
#   random = ~id,
#   data = dat,
#   control = mmes_control(update = list(recovery = "pxem", blendWeight = 5)))

