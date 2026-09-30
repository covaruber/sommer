\name{stackTraits}
\alias{stackTraits}
\title{Stack several trait columns into the long format used by multi-trait models}
\description{
Converts a wide data frame (one column per trait) into the long format used by
multi-trait \code{\link{mmes}} models: one row per record and trait, a factor
identifying the trait, a single response column and a record key that pairs
the traits measured on the same experimental unit.
}
\usage{
stackTraits(data, traits, keep = NULL, trait = "trait", value = "value",
            record = "record")
}
\arguments{
  \item{data}{A data frame in wide format.}
  \item{traits}{Character vector with the names of the trait columns to stack.
  Their order defines the levels of the trait factor.}
  \item{keep}{Columns of \code{data} copied to every stacked row. By default,
  all columns that are not in \code{traits}.}
  \item{trait}{Name of the new trait factor column.}
  \item{value}{Name of the new response column.}
  \item{record}{Name of the new record-key factor column.}
}
\details{
Rows are stacked trait by trait, so row \eqn{i} of trait \eqn{t} is row
\eqn{(t-1)n+i} of the output. The record key is taken from the row names of
\code{data} when these are unique and not the default \code{1..n}; otherwise
it is \code{"r1", "r2", ...}. Missing trait values are kept as \code{NA} (they
are handled by the \code{naMethodY} rules of \code{mmes}, so incomplete records
are allowed).

Using the record key in the residual,
\code{rcov = ~ vsm(usm(trait), ism(record))}, pairs the residuals of the traits
of each unit explicitly. With \code{ism(units)} the pairing is implied by the
order of the rows within each trait instead.
}
\value{A data frame with the \code{keep} columns and the \code{record},
\code{trait} and \code{value} columns, with \code{length(traits) * nrow(data)}
rows.}
\seealso{\code{\link{mmes}}, \code{\link{covmatrix_mmes}}, \code{\link{familym}}}
\examples{
set.seed(17)
n <- 60
DT <- data.frame(id = factor(rep(seq_len(n), each = 3)))
u <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, 0.35, 0.35, 1.4), 2))
DT$trait_a <- u[DT$id, 1] + rnorm(nrow(DT), sd = 0.7)
DT$trait_b <- u[DT$id, 2] + rnorm(nrow(DT), sd = 1)
DTL <- stackTraits(DT, traits = c("trait_a", "trait_b"), keep = "id")
head(DTL)
table(DTL$trait)

# bivariate model with the residual paired through the record key
fit <- mmes(value ~ trait,
            random = ~ vsm(usm(trait), ism(id)),
            rcov = ~ vsm(usm(trait), ism(record)),
            data = DTL, verbose = FALSE)
covmatrix_mmes(fit, 1)$correlation   # genetic correlation
}
\keyword{manip}
