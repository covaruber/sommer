\name{dsumm}
\alias{dsumm}
\title{Direct sum of residual covariance structures across sections}
\description{
\code{dsumm} replicates a residual \code{\link{vsm}} structure independently for
every level (section) of a grouping variable, such as a trial, location, or
experiment. Each section receives its own variance and its own copy of every
covariance parameter of the inner structure. This is the sommer counterpart of
the ASReml specification \code{dsum(~ ... | section)}.
}
\usage{
dsumm(x, by, levels = NULL)
}
\arguments{
  \item{x}{A residual \code{vsm()} term, for example
    \code{vsm(ar1m(range), ar1m(row), ism(units))}. Its covariance-shaping
    factors define the structure fitted within each section.}
  \item{by}{Observation-level grouping variable defining the sections. Missing
    values of \code{by} exclude the observation.}
  \item{levels}{Optional character vector of \code{by} levels that receive the
    structured covariance of \code{x}. All other levels receive an independent
    and identically distributed residual with their own variance. The default
    \code{NULL} applies the structure to every level.}
}
\details{
\strong{Direct sum versus Kronecker product.} A \code{vsm()} term is a Kronecker
product with a single overall variance. Including a diagonal section factor,
\deqn{\mathrm{vsm}(\mathrm{dsm}(t), \mathrm{ar1m}(r), \mathrm{ar1m}(c), \mathrm{ism}(units)):\quad
R = \sigma^2\, D_t \otimes K_r(\rho_r) \otimes K_c(\rho_c),}
gives each section its own variance, but a single correlation shared by all
sections, because each factor of a Kronecker product appears only once. In
contrast, \code{dsumm} builds the block-diagonal direct sum
\deqn{R = \bigoplus_{t} \sigma^2_t\, K_r(\rho_{r,t}) \otimes K_c(\rho_{c,t}),}
so every section has its own variance and its own correlations. Sections are
mutually independent, so their parameters are estimated only from their own
observations; with section-specific fixed effects the joint fit reproduces
separate per-section fits.

\strong{Choosing between the two.} Use the \code{vsm()} form when trials are
small or similar and a common spatial correlation is adequate (fewer
parameters). Use \code{dsumm} when correlations are expected to differ between
trials. With \eqn{T} sections and \eqn{p} inner covariance parameters,
\code{dsumm} estimates \eqn{T(p+1)} residual parameters.

\strong{Unlisted levels.} Levels omitted from \code{levels} are modelled as
\eqn{\sigma^2_t I}. Their observations are retained even when the variables of
the inner structure (for example row or range) are missing, so trials without
spatial coordinates can be analysed together with spatially designed trials.

\strong{Computation.} Each section is embedded in the rectangle spanned by its
observed levels of the inner factors. Covariance applications use the small
factor inverses, and missing cells (e.g. plots with a missing response) are
absorbed exactly by a Schur complement, so the dense observation-level
covariance matrix is never formed. This requires sections that are close to
complete rectangles and no \code{W} weights.

\strong{Reporting.} \code{\link{covparams_mmes}} returns one row per parameter
and section, with the section in the \code{section} column. The variance rows
belong to the \code{by} factor; correlation rows are labelled by their inner
factor, e.g. \code{ar1m(range)}.
}
\value{
A residual covariance structure, in the same form as a residual \code{vsm()}
term, for use in the \code{rcov} argument of \code{\link{mmes}}.
}
\section{Limitations}{
\code{dsumm} is currently available only for residual (\code{rcov}) structures,
without \code{W} weights, and cannot be nested.
}
\examples{
\dontrun{
# ASReml: residual = ~dsum(~ar1(range):ar1(row) | trial)
mix <- mmes(yield ~ trial,
            rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by = trial),
            data = dat)
covparams_mmes(mix)

# Spatial structure only for two trials; remaining trials are iid.
mix2 <- mmes(yield ~ trial,
             rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)),
                            by = trial, levels = c("T1", "T2")),
             data = dat)

# Site-specific unstructured covariance among environments.
mix3 <- mmes(y ~ site,
             rcov = ~ dsumm(vsm(usm(env), ism(units)), by = site),
             data = dat2)
}
}
\seealso{
\code{\link{vsm}}, \code{\link{ar1m}}, \code{\link{dsm}}, \code{\link{mmes}},
\code{\link{covparams_mmes}}.
}
