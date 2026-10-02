\name{APY}
\alias{APY}
\title{Construct an APY Relationship Precision Matrix}
\description{
Constructs the Algorithm for Proven and Young (APY) approximation for a
positive-definite relationship matrix and a user-specified core set. In the
genomic evaluation literature, the input is the genomic relationship matrix
\eqn{G} for genotyped animals, usually calculated from SNP marker data. The
core and non-core partition is within that genotyped set; it is not a
pedigree-versus-genomic partition. The non-core conditional covariance is
approximated by its diagonal. The returned matrix is a sparse precision matrix
suitable for the \code{Gu} argument of \code{mmes()} in a genomic BLUP model.
}
\usage{
APY(G, core, tol = 1e-10, return.details = FALSE)
}
\arguments{
  \item{G}{A finite, symmetric, positive-definite numeric relationship matrix.
  For the standard APY genetic evaluation use, this is the genomic relationship
  matrix among genotyped animals, commonly obtained from marker genotypes with
  \code{A.mat()}. Row and column names, when present, must agree exactly. This
  implementation receives \code{G} as a dense matrix.}
  \item{core}{Unique one-based row indices or row names selecting the core
  individuals.}
  \item{tol}{Positive tolerance below one, used for symmetry validation and
  the minimum acceptable conditional variance relative to the relationship
  scale. Small or non-positive conditional variances cause an error; they are
  not silently regularized.}
  \item{return.details}{If \code{TRUE}, return a list containing the sparse
  precision matrix as \code{Gu}, core and non-core indices, regression matrix
  \code{B}, and non-core conditional variances.}
}
\value{
A sparse symmetric APY precision matrix with individual dimnames and
\code{attr(x, "inverse") = TRUE}, ready for \code{Gu}. Metadata describing
the core and conditional variances is stored in \code{attr(x, "APY")}. If
\code{return.details = TRUE}, returns the diagnostic list described above.
}
\details{
Partition \eqn{G} into core and non-core blocks. With
\eqn{B = G_{nc}G_{cc}^{-1}} and
\eqn{D = \operatorname{diag}(G_{nn} - BG_{cn})}, APY uses the precision
\eqn{\left[\begin{smallmatrix}
G_{cc}^{-1} + B' D^{-1} B & -B' D^{-1} \\
-D^{-1} B & D^{-1}
\end{smallmatrix}\right]}. It preserves the core block and core-to-non-core
relationships while replacing non-core conditional correlations with zero.
Thus APY is exact when the non-core conditional covariance is diagonal;
otherwise it is an approximation to the supplied relationship model.

APY was developed to provide a computationally economical approximation to
the inverse genomic relationship matrix for large genomic evaluations
(Misztal, 2016; Fragomeni et al., 2015). The algorithm partitions genotyped
animals into core and non-core individuals. The core contains the relationships
that are modeled jointly; each non-core animal is regressed on the core, with
its remaining conditional variance treated as independent of other non-core
animals. Thus non-core animals are not independent of the core: their
core-to-non-core relationships are retained.

Pedigree relationships are a distinct input. For pedigree-only models, the
pedigree relationship matrix \eqn{A} and its sparse inverse \eqn{A^{-1}} are
normally used directly; APY is not required. In single-step genomic BLUP,
pedigree and genomic information are combined, conventionally through
\eqn{H^{-1} = A^{-1} + \left[\begin{smallmatrix}0&0\\0&G^{-1}-A_{22}^{-1}\end{smallmatrix}\right]},
where genotyped animals are ordered last and \eqn{A_{22}} is their pedigree
relationship block. APY approximates the genomic term \eqn{G^{-1}} in that
correction; it does not replace the pedigree term or construct \eqn{H^{-1}}.
This function returns only the APY genomic precision. Do not use it as the
complete single-step precision for a population that also contains
non-genotyped animals.

This function constructs APY from an already available relationship matrix. It
does not avoid the cost of constructing a full dense relationship matrix from
markers, and it does not change the solver algorithm. Marker-based block
construction and single-step matrix assembly are not implemented here.
}
\examples{
## Marker-derived GBLUP for genotyped animals.
set.seed(2)
ids <- paste0("animal", seq_len(12))
markers <- matrix(sample(c(-1, 0, 1), 12 * 80, replace = TRUE), nrow = 12,
                  dimnames = list(ids, NULL))
G <- A.mat(markers)

## Blend slightly with the identity so this small illustrative G is positive
## definite, then construct APY precision for a six-animal core.
G <- 0.95 * G + 0.05 * diag(nrow(G))
dimnames(G) <- list(ids, ids)
Ginv.apy <- APY(G, core = ids[seq_len(6)])

phenotypes <- data.frame(
  id = factor(rep(ids, each = 3), levels = ids),
  y = rnorm(length(ids) * 3)
)
fit <- mmes(y ~ 1, random = ~vsm(ism(id), Gu = Ginv.apy),
            rcov = ~units, data = phenotypes, verbose = FALSE)

## In a single-step evaluation, this APY approximation supplies G^{-1} for
## the genotyped block correction; it is not itself the complete H^{-1}.
}
\references{
Misztal I (2016). Inexpensive Computation of the Inverse of the Genomic
Relationship Matrix in Populations with Small Effective Population Size.
\emph{Genetics}, 202(2), 401--409. doi:10.1534/genetics.115.182089.

Fragomeni BO, Lourenco DA, Tsuruta S, Masuda Y, Aguilar I, Legarra A, Lawlor
TJ, Misztal I (2015). Use of genomic recursions and algorithm for proven and
young animals for single-step genomic BLUP analyses--a simulation study.
\emph{Journal of Animal Breeding and Genetics}, 132(5), 340--345.
doi:10.1111/jbg.12161.
}
\seealso{\code{\link{mmes}}, \code{\link{A.mat}}, \code{\link{H.mat}}}