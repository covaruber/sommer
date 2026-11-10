\name{graph_term}
\alias{graph_term}
\title{Reference a named graph predictor in a model formula}
\usage{
graph_term(name)
}
\arguments{
\item{name}{A literal, non-empty character name in the fitting call's \code{predictors} list.}
}
\description{
A formula-only marker referencing a \code{nonlinear_graph} object in
\code{predictors}. It must not be evaluated outside a model formula.
}
\value{No standalone value; the fitting layer interprets the formula call.}
\details{
Graph coordinates have namespaced parameter names \code{name::coordinate}.
They default to intercept-only identity-linked declarations; named entries in
\code{parameters} can override these defaults with \code{parameter()}.
Natural starts use the same namespaced names. Graph-internal random maps and
joint random declarations are not supported by this route; use namespaced
parameter random regressions instead. Unused predictors are errors.
Graph fits require a Gaussian identity family and an explicit identity output
node. The existing graph evaluator supplies graph values to the same expression
and likelihood route used by explicit parameter expressions.
}
\examples{
f <- y ~ Env + graph_term("net")
f
}
\seealso{\code{\link{parameter}}, \code{\link{nonlinear_graph}}}
