\name{postVarU}
\alias{postVarU}
\title{Post-fit Sampling Covariance of BLUPs}
\description{
Computes the plug-in sampling variance or covariance of fitted random-effect
predictors, \eqn{\hat u}, using Henderson's mixed-model coefficient system.
This is distinct from the prior variance of the random effects and their
prediction error variance (PEV).
}
\usage{
postVarU(object, mode=1L)
}
\arguments{
\item{object}{A Gaussian \code{mmes} Henderson fit with stored \code{C} and
\code{Cscale}. Direct-engine, PQL and matrix-free \code{solveOnly} fits are
not supported.}
\item{mode}{\code{1} (default) returns diagonal BLUP sampling variances in
\code{uVarList}; \code{2} additionally returns the joint sampling covariance
matrix \code{VarU}; \code{0} clears only VarU-derived outputs.}
}
\details{
For \eqn{y=X\beta+Zu+e}, with \eqn{\operatorname{Var}(u)=G} and
\eqn{\operatorname{Var}(e)=R}, the unscaled Henderson matrix is
\deqn{C_0=\left[\begin{array}{cc}
X'R^{-1}X & X'R^{-1}Z \\
Z'R^{-1}X & Z'R^{-1}Z+G^{-1}
\end{array}\right].}
Its random-effect inverse block is
\eqn{Q=(C_0^{-1})_{uu}=\operatorname{Var}(u-\hat u)}. Consequently,
\deqn{\operatorname{Var}(\hat u)=G-Q.}
For known covariance parameters this identity is exact. At estimated
parameters the result is a plug-in covariance, omitting covariance-parameter
estimation uncertainty. The same identity is used for ML and REML fits.

The implementation does not construct or invert the observation covariance
\eqn{V=ZGZ'+R}. Error covariance is computed by solving the stored Henderson
system and multiplying by \code{Cscale} to recover the original response
scale. For each independent random term,
\eqn{G_k=\Sigma_k\otimes A_k}, where \eqn{\Sigma_k} is the fitted
covariance-coordinate matrix and \eqn{A_k} the relationship covariance.
Products with \eqn{A_k} are evaluated by solving its retained sparse
precision, not by explicitly forming its inverse or the full joint \eqn{G}.

\code{mode=1} processes identity contrasts in batches of at most 128 rows,
so neither a full coefficient inverse nor a full random-effect covariance
matrix is stored. It requires solves for all requested effect diagonals and
can be expensive for large models. \code{mode=2} requires quadratic storage
in the number of random effects. Although independent terms have zero prior
cross-covariance, their BLUPs can have nonzero sampling cross-covariance; the
joint \code{VarU} includes these cross-term entries.

Rotation and factor-score augmentation are mapped to the public
\code{uList} coordinates before computing error covariance. New fits retain
\code{randomPrecision} in original level order, so external \code{Gu}
objects are not needed after fitting. Older ordinary fits attempt to rebuild
relationship precisions from saved model inputs, requiring the original
formula objects; dimension checks cannot detect all changes to those objects.
Older rotated or augmented fits must be refitted. No clipping is performed:
small negative variances from numerical cancellation remain visible.

Existing \code{Ci} and \code{uPevList} are left unchanged. The variance here
is not the empirical variance across BLUP values and is not a posterior
variance of \eqn{u}. It is the repeated-sampling covariance of \eqn{\hat u}.
}
\value{
The original fitted object with \code{uVarList}, a named list matching the
dimensions and labels of \code{uList}, and \code{uVarMode}. With
\code{mode=2}, \code{VarU} is the joint covariance in term order, with each
\code{uList} matrix stacked column by column. Its labels contain term,
covariance-coordinate and main-effect labels. Modes 0 and 1 remove any
previous full \code{VarU}; mode 0 also removes \code{uVarList}.
}
\seealso{\code{\link{postPEV}}, \code{\link{predict.mmes}},
\code{\link{covmatrix_mmes}}}
\examples{
data(DT_example, package="enhancer")
model <- mmes(Yield~Env, random=~Name, data=DT_example,
              computeCi=0, verbose=FALSE)
model <- postVarU(model)
head(model$uVarList[[1]])
full <- postVarU(model, mode=2)
head(diag(full$VarU))

prediction <- predict(model, D="Name", PEV=TRUE, VarU=TRUE)
head(prediction$pvals)
prediction$PEV[1:3, 1:3]
prediction$VarU[1:3, 1:3]
prediction$sampling.vcov[1:3, 1:3]
}