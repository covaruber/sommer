\name{mmes_control}
\alias{mmes_control}
\title{Modular Numerical Controls for Mixed Models}
\description{
Creates a strictly validated \code{sommer_mmes_control} object for the
numerical methods used by \code{\link{mmes}}. Model formulas remain separate.
Every module accepts only its documented fields; unknown module or field names
are errors. Constructing controls does not fit a model or establish that
requested methods are eligible for a particular model.
}
\usage{
mmes_control(formulation = "henderson", task = list(),
             representation = list(), solve = list(), logdet = list(),
             score = list(), curvature = list(), update = list(),
             proposal = list(), convergence = list(), output = list(),
             execution = list(), diagnostics = list())
}
\arguments{
\item{formulation}{Equation domain: \code{"henderson"} (default) or
\code{"direct"}. Henderson forms coefficient-space mixed-model equations;
direct uses observation-space inversion.}
\item{task}{Named list with \code{mode}, \code{objective}, and \code{pql}.
Defaults: \code{mode="fit"}, \code{objective="reml"}, and
\code{pql=list(maxit=20L, tol=1e-5, estimateTheta=FALSE, warmStart=TRUE)}.
\code{mode} is \code{"fit"} or \code{"solve"}; \code{objective} is
\code{"reml"} or \code{"ml"}. The nested \code{pql} list accepts only
\code{maxit}, \code{tol}, \code{estimateTheta}, and \code{warmStart}, with
the defaults shown. These fields own task selection and outer PQL iteration,
not inner covariance iterations.}
\item{representation}{Named list with \code{compression="off"},
\code{factorScores="none"}, \code{equations="normal"},
\code{sparseFactors="off"}, and \code{matrixFree="auto"} by default,
unless corresponding \code{sommer.mme.*} options override experimental modes.
If \code{sommer.mme.saddlePoint} resolves to \code{"force"}, the default
\code{equations} becomes \code{"saddle-point"}.
\code{compression} is \code{"off"}, \code{"groups"}, or \code{"classes"};
\code{factorScores} is \code{"none"}, \code{"fixed-shape"}, or
\code{"profile"}; \code{equations} is \code{"normal"} or
\code{"saddle-point"}. The remaining fields use experimental mode values
\code{"auto"}, \code{"off"}, or \code{"force"} (logical TRUE/FALSE are
accepted). This module owns equation, compression, covariance-factor, and
factor-score representations.}
\item{solve}{Named list. \code{method="auto"} selects among
\code{"auto"}, \code{"ldlt"}, \code{"cholmod"}, and \code{"pcg"}.
PCG fields default to \code{tol=1e-8}, \code{maxIters=0L},
\code{preconditioner="diagonal"}, \code{nystromRank=32L},
\code{warmStart=FALSE}, \code{recycleRank=0L},
\code{recycleMethod="coarse"}, \code{recycleState=NULL},
\code{exportRecycle=FALSE}, \code{rhsMethod="independent"},
\code{blockSize=16L}, and \code{residualRefresh=0L}.
Also accepts the experimental solve requests \code{blockSchur="auto"},
\code{blockBanded="auto"}, \code{factorSchurMultiFactor="auto"}, and
\code{dFirstML="auto"}; global \code{sommer.mme.*} options may override
these mode defaults. The method owns checked linear solves,
preconditioning, reuse, and eligible elimination. PCG settings apply only
when an eligible PCG path uses them.}
\item{logdet}{Named list with \code{scaling="none"} and
\code{lanczosSteps=20L} by default. Scaling is \code{"none"} or
\code{"diagonal"}. This module owns determinant approximation and its
correction, not linear-solve accuracy.}
\item{score}{Named list with \code{traceProbes=8L},
\code{estimator="hutchinson"}, \code{rank=16L},
\code{covarianceChecks="auto"}, and \code{traceBatch="off"} by default.
The estimator is \code{"hutchinson"} or \code{"control-variate"}.
\code{covarianceChecks} and \code{traceBatch} accept experimental mode
values. This module owns covariance derivatives, trace estimation, and
covariance checks.}
\item{curvature}{Named list with \code{method="exact-ai"} and
\code{reduced="off"} by default. The method is \code{"exact-ai"} or
\code{"hybrid-ai"}; \code{reduced} accepts experimental mode values.
This module owns the AI metric and reduced-direction assembly.}
\item{update}{Named list with \code{recovery="em"},
\code{blendWeight=NULL}, and \code{commonScale="off"} by default.
\code{recovery} is \code{"em"} or \code{"pxem"}; \code{commonScale}
accepts experimental mode values. This module owns covariance-coordinate
updates, EM/PXEM recovery, blending, and common-scale profiling.}
\item{proposal}{Named list with \code{stepWeight=NULL} and
\code{extrapolation="none"} by default. Extrapolation is \code{"none"} or
\code{"aitken"}. This module owns proposal scaling and extrapolation.
Aitken extrapolation cannot be combined with
\code{curvature$method="hybrid-ai"}.}
\item{convergence}{Named list with \code{maxIters=30L},
\code{likelihoodTol=1e-4}, \code{parameterTol=1e-4},
\code{accuracy="fixed"}, \code{maxTraceProbes=NULL},
\code{maxLanczosSteps=NULL}, \code{logDetTol=1e-3}, and
\code{scoreTol=1e-3}. Accuracy is \code{"fixed"} or \code{"adaptive"};
NULL stochastic limits resolve to the initial \code{traceProbes} and
\code{lanczosSteps}. This module owns inner iteration convergence and
stochastic accuracy requests.}
\item{output}{Named list with \code{inverse=0L}, \code{covariance="full"},
\code{coefficientMatrix=TRUE},
\code{setup=FALSE}, \code{selectedInverse="auto"}, and
\code{panelSize=64L} by default. \code{inverse} is 0 (no final coefficient
inverse or fit-time PEV computation),
1 (selected inverse and random-effect PEVs), or 2 (full inverse and
random-effect PEVs). Modes 1 and 2 compute PEVs during fitting; mode 1
does not form the complete coefficient-matrix inverse. With mode 0, use
\code{\link{postPEV}} after fitting to compute PEVs or the full inverse,
provided \code{coefficientMatrix=TRUE};
\code{covariance} is \code{"full"} or \code{"none"}; \code{selectedInverse}
accepts experimental mode values. This module owns coefficient uncertainty,
covariance-parameter uncertainty, retained coefficient-matrix output, and
setup output.}
\item{execution}{Named list with \code{denseMemoryMB=1000},
\code{memoryPolicy="speed"}, \code{gpu=FALSE}, \code{gpuMinSize=256L},
\code{blasThreads=NULL}, \code{inverseTol=1e-6}, and \code{seed=0L} by
default. A non-NULL \code{sommer.mme.denseMemoryMB},
\code{sommer.mme.memoryPolicy}, \code{sommer.mme.gpu}, or
\code{sommer.mme.gpuMinSize} option supplies its corresponding default.
This module owns resource policy, hardware, threading, numerical inverse
tolerance, and stochastic seed.}
\item{diagnostics}{Named list with \code{verbose=TRUE},
\code{dateWarning=TRUE}, and \code{profile="off"} by default.
\code{profile} accepts experimental mode values. This module owns ordinary
fit messages, package-age warnings, and phase profiling.}
}
\details{
Module ownership describes mathematical responsibility, not a promise that
every cross-module combination is supported. For example, formulation and
solver choices are model-dependent; PCG has REML and output restrictions;
compression, saddle-point equations, factor-score augmentation, common-scale
profiling, and reduced AI each have additional prerequisites. Consult
\code{vignette("sommer.mmes.options", package="sommer")} and the internal
\code{.mmes_capabilities()} catalog for the current conceptual decision graph
and capability inventory. The catalog is descriptive; fitting validation
remains authoritative. The modular API does not imply a kernel rewrite or
arbitrary method composition.

\code{task$mode="solve"} selects a known-covariance path rather than ordinary
variance-component fitting; it requires an eligible Gaussian identity model
and has narrower PCG/output options. \code{task$pql$estimateTheta=TRUE}
requires a negative-binomial family. \code{representation$compression} needs
eligible matching designs and independent residuals. Saddle-point equations,
factor-score representation, sparse factors, and matrix-free PCG each have
separate model and backend gates. PCG fitting is a Henderson REML path;
selected-inverse output mode 1 is not supported by PCG. The vignette and
capability catalog describe these prerequisites, but the fit-time validator
decides eligibility.

Within PCG, \code{solve$preconditioner="nystrom"} uses
\code{solve$nystromRank}; \code{solve$recycleState} and
\code{solve$exportRecycle} require \code{solve$recycleMethod="krylov"} and a
positive \code{solve$recycleRank}. \code{logdet$scaling} and Lanczos depth
apply to stochastic determinant evaluation. \code{score$rank} is used by the
control-variate estimator. Adaptive trace and Lanczos budgets apply to
stochastic PCG accuracy. \code{output$panelSize} applies to full inverse
panels, not selected-inverse mode.

Experimental mode fields accept \code{"auto"}, \code{"off"}, or
\code{"force"}; logical TRUE and FALSE mean force and off. Package options
may supply defaults for these requests. \code{"force"} requests a path but
does not bypass mathematical, structural, or memory eligibility checks.
Nested lists reject unknown names and duplicate or empty names.

Defaults retain the established numerical bundle. Aitken belongs to
\code{proposal$extrapolation}; hybrid AI belongs to
\code{curvature$method}; they are mutually exclusive.

\strong{Migration from removed numerical arguments}

\code{mmes()} now takes numerical settings through \code{control}. The old
formals are not compatibility aliases and must not be passed to \code{mmes()}.
Map former settings to canonical module paths:
\describe{
\item{\code{henderson}}{\code{formulation="henderson"} or
\code{formulation="direct"}.}
\item{\code{solver}, \code{engine}}{\code{solve$method}; pass the complete
control object as \code{control=...}.}
\item{\code{nIters}, \code{tolParConvLL}, \code{tolParConvNorm},
\code{tolParInv}}{\code{convergence$maxIters},
\code{convergence$likelihoodTol}, \code{convergence$parameterTol}, and
\code{execution$inverseTol}, respectively.}
\item{\code{REML}, \code{solveOnly}, \code{pqlControl}}{\code{task$objective}
(\code{"reml"} or \code{"ml"}), \code{task$mode}
(\code{"fit"} or \code{"solve"}), and \code{task$pql}.}
\item{\code{computeCi}, \code{getPEV}, \code{returnParam}}{\code{output$inverse}
and \code{output$setup}. PEV output is controlled by \code{output$inverse};
the former \code{getPEV} flag has no separate replacement.}
\item{\code{blend}, \code{blendWeight}, \code{stepWeight}}{\code{update$recovery},
\code{update$blendWeight}, and \code{proposal$stepWeight}.}
\item{\code{acceleration="aitken"}}{\code{proposal$extrapolation="aitken"}.}
\item{\code{acceleration="hybrid-ai"}}{\code{curvature$method="hybrid-ai"}.
Do not combine with Aitken.}
\item{\code{compression}, \code{factorScoreAugmentation},
\code{covarianceInformation}, \code{retainCoefficientMatrix}}{
\code{representation$compression}, \code{representation$factorScores},
\code{output$covariance}, and \code{output$coefficientMatrix}.}
\item{Former \code{pcg} settings}{Map solve tolerances, iteration limits,
preconditioning, and reuse to \code{solve}; determinant scaling and Lanczos
depth to \code{logdet}; trace probes, estimator, rank, and covariance checks
to \code{score}; stochastic accuracy budgets to \code{convergence}; seed to
\code{execution}.}
\item{Former \code{experimental} settings}{Map structured representation
requests to \code{representation}, elimination requests to \code{solve},
trace batching to \code{score}, reduced AI to \code{curvature},
common-scale profiling to \code{update}, selected/full inverse requests to
\code{output}, and profiling to \code{diagnostics}.}
\item{\code{denseMemoryMB}, \code{memoryPolicy}, \code{gpu},
\code{gpuMinSize}, \code{blasThreads}}{\code{execution} fields with the same
names.}
\item{\code{verbose}, \code{dateWarning}}{\code{diagnostics} fields with the
same names.}
\item{\code{recordDiagnostics}, \code{verboseDiagnostics}}{Removed flags with
no direct replacement. Use \code{diagnostics$verbose} for fit messages and
\code{diagnostics$profile} for implemented phase profiling.}
}
With \code{diagnostics$verbose=TRUE}, the existing Henderson \code{MME path}
line includes the coefficient solver, active/inactive latent reduced-rank and
Schur variants, and the selected residual engine (including inactive grid and
repeated-Kronecker paths). An inactive automatic Schur path is reported as
``not selected'', not as proof that the model lacks exploitable structure;
disabled controls and unsupported coefficient backends are distinguished.
}
\section{Compatibility}{
Settings are checked against one compatibility table when the control is
built and again in \code{\link{mmes}}. A supplied setting that conflicts with
another, or that would have no effect on the selected path, stops the call with
a message that names every offending setting and how to fix it. A setting counts
as supplied only when it differs from its default, so restating defaults is
always accepted. The main rules are:
\itemize{
\item \code{formulation="direct"} accepts only \code{task}, \code{update$recovery},
\code{update$blendWeight},
\code{proposal$stepWeight}, convergence iterations and tolerances,
\code{output$setup}, \code{output$inverse} 0 or 2, memory
and threading \code{execution} settings, and \code{diagnostics}.
\item \code{task$mode="solve"} accepts only \code{solve$tol}, \code{solve$maxIters},
\code{solve$warmStart}, \code{solve$residualRefresh}, a non-Nystrom
\code{solve$preconditioner}, \code{output$setup}, CPU \code{execution} settings and
\code{diagnostics}.
\item PCG-only settings (solve tolerances, preconditioner, recycling,
\code{logdet}, \code{score} estimators, adaptive accuracy,
\code{representation$matrixFree}, \code{execution$seed}) require
\code{solve$method="pcg"}; \code{"auto"} never selects PCG.
\item \code{solve$method="pcg"} rejects EM/PXEM recovery, blending,
common-scale profiling, reduced AI, Aitken, hybrid AI, saddle-point
equations, factor scores and \code{task$objective="ml"}.
\item Aitken extrapolation, hybrid AI, saddle-point equations, factor scores,
blending and \code{output$coefficientMatrix=FALSE} each reject the modules
they cannot combine with.
}
Rules that depend on the data, family, \code{vcc} or rotation are checked by
\code{\link{mmes}} at fit time; for example a Gaussian identity fit stops when
\code{task$pql} is changed.
}
\value{
An object of class \code{"sommer_mmes_control"}.
}
\seealso{
\code{\link{mmes}}; \code{vignette("sommer.mmes.options", package="sommer")}
}
\examples{
ctrl <- mmes_control(
  formulation = "henderson",
  solve = list(method = "pcg", tol = 1e-8),
  score = list(traceProbes = 16L),
  convergence = list(accuracy = "adaptive")
)
print(ctrl)
}
