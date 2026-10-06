# Post-fit Sampling Covariance of BLUPs

`postVarU()` operates on a Gaussian `mmes` Henderson model object. It
returns the model with sampling variances of its BLUPs, not prediction
error variances and not the empirical variance across reported effects.

```r
postVarU(object, mode=1L)
```

## Modes and Outputs

- `mode=1`: `uVarList` contains effect-level sampling variances with the
  same dimensions and labels as `uList`. No full covariance is stored.
- `mode=2`: additionally returns joint covariance `VarU`, including
  cross-term entries. Effects are ordered by term and each `uList` matrix
  stacked column by column. Labels include term, coordinate and main-effect
  level. Storage is quadratic in the number of random effects.
- `mode=0`: removes only `uVarList` and `VarU`.

`uVarMode` records the mode. Existing `Ci` and `uPevList` remain unchanged.

## Henderson Calculation

For $y=X\beta+Zu+e$, with $\operatorname{Var}(u)=G$ and
$\operatorname{Var}(e)=R$, define

$$
C_0=\begin{bmatrix}
X'R^{-1}X & X'R^{-1}Z\\
Z'R^{-1}X & Z'R^{-1}Z+G^{-1}
\end{bmatrix}.
$$

Then $Q=(C_0^{-1})_{uu}=\operatorname{Var}(u-\hat u)$ and

$$\operatorname{Var}(\hat u)=G-Q.$$

The implementation solves the stored Henderson system, applies `Cscale`,
and computes prior covariance products using sparse relationship-precision
solves and $G_k=\Sigma_k\otimes A_k$. It never constructs or inverts the
observation-level covariance $V$, nor forms the full joint $G$.
Mode 1 uses batches of at most 128 identity contrasts; obtaining all
diagonals may still be expensive for large models.

Independent prior terms can have correlated BLUPs, so mode 2 includes
their sampling cross-covariances. Rotation and factor-score effects are
mapped to public coordinates. New fits retain original relationship
precisions. Older ordinary fits attempt input reconstruction, requiring
unchanged formula objects; dimension checks cannot detect all changes.
Older transformed fits must be refitted.

These covariances are exact at known parameters and plug-in approximations
at estimated parameters, for ML and REML, excluding variance-component
estimation uncertainty. Tiny negative variances from cancellation are not
clipped. Direct, PQL and matrix-free `solveOnly` fits are unsupported.

## Example

```r
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
```

`postPEV()` and `postVarU()` attach model-effect outputs. In contrast,
`predict()` projects covariance through the final `D` specified by `D`,
`Dtable` and `levels`. Prediction `VarU` is for the random contribution;
`sampling.vcov` additionally includes fixed-effect sampling uncertainty.
