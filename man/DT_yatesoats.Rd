\name{DT_yatesoats}
\alias{DT_yatesoats}
\docType{data}
\title{
  Yield of oats in a split-block experiment
}
\description{
  The yield of oats from a split-plot field trial using three varieties
  and four levels of manurial treatment. The experiment was laid out in
  6 blocks of 3 main plots, each split into 4 sub-plots. The varieties
  were applied to the main plots and the manurial (nitrogen) treatments
  to the sub-plots.
}

\format{
  \describe{
    \item{\code{block}}{block factor with 6 levels}
    \item{\code{nitro}}{nitrogen treatment in hundredweight per acre}
    \item{\code{Variety}}{genotype factor, 3 levels}
    \item{\code{yield}}{yield in 1/4 lbs per sub-plot, each 1/80 acre.}
    \item{\code{row}}{row location}
    \item{\code{column}}{column location}
  }
}

\source{
  Yates, Frank (1935) Complex experiments,
  \emph{Journal of the Royal Statistical Society Suppl}. 2, 181--247.
}

\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
}

\examples{

### ========================== ###
### ========================== ###
data(DT_yatesoats)
DT <- DT_yatesoats
head(DT)
# m3 <- mmer(fixed=Y ~ V + N + V:N,
#            random = ~ B + B:MP,
#            rcov=~units,
#            data = DT)
# summary(m3)$varcomp
# 
# m3b <- mmec(fixed=Y ~ V + N + V:N,
#            random = ~ B + B:MP,
#            rcov=~units,
#            data = DT)
# summary(m3b)$varcomp
}

