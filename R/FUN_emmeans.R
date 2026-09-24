# ##############################################################################
# #    Copyright (c) 2012-2019 Russell V. Lenth                                #
# #                                                                            #
# #    This file is part of the emmeans package for R (*emmeans*)              #
# #                                                                            #
# #    *emmeans* is free software: you can redistribute it and/or modify       #
# #    it under the terms of the GNU General Public License as published by    #
# #    the Free Software Foundation, either version 2 of the License, or       #
# #    (at your option) any later version.                                     #
# #                                                                            #
# #    *emmeans* is distributed in the hope that it will be useful,            #
# #    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# #    GNU General Public License for more details.                            #
# #                                                                            #
# #    You should have received a copy of the GNU General Public License       #
# #    along with R and *emmeans*.  If not, see                                #
# #    <https://www.r-project.org/Licenses/> and/or                            #
# #    <http://www.gnu.org/licenses/>.                                         #
# ##############################################################################
# 
# # sommer package support
# 
# recover_data.mmer = function(object, data, ...) {
#     if (is.null(data))
#         data = object$data
#     fcall = call("mmer", formula = object$call$fixed, data = data)
#     emmeans::recover_data(fcall, delete.response(terms(object$call$fixed)), 
#                           object$call$na.method.V, ...)
# }
recover_data.mmes <- function(object, data=NULL, ...) {
	fixed <- object$args$fixed
	if(is.null(fixed))
		return("The fitted mmes object does not retain its fixed-effects formula.")

	if(is.null(data)) data <- object$dataOriginal
	if(is.null(data))
		return("The fitted mmes object does not retain the original data; supply it via data =.")

	fitCall <- call("mmes", fixed=fixed, data=data)
	emmeans::recover_data(fitCall, stats::delete.response(stats::terms(fixed)),
												na.action=NULL, ...)
}

emm_basis.mmes <- function(object, trms, xlev, grid, ...) {
	if(ncol(object$y) != 1L)
		stop("emmeans support for mmes currently requires a univariate response.", call.=FALSE)

	betaNames <- rownames(object$b)
	if(is.null(betaNames))
		stop("The fitted mmes object does not retain fixed-effect names.", call.=FALSE)

	modelFrame <- stats::model.frame(trms, grid, na.action=stats::na.pass, xlev=xlev)
	X <- stats::model.matrix(trms, modelFrame)
	colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"

	missingColumns <- setdiff(betaNames, colnames(X))
	if(length(missingColumns)) {
		X <- cbind(X, matrix(0, nrow=nrow(X), ncol=length(missingColumns),
												 dimnames=list(NULL, missingColumns)))
	}
	unexpectedColumns <- setdiff(colnames(X), betaNames)
	if(length(unexpectedColumns))
		stop("The emmeans reference grid produced fixed-effect columns absent from the fitted mmes model.",
				 call.=FALSE)
	X <- X[, betaNames, drop=FALSE]

	nCoefficients <- nrow(object$bu)
	if(is.null(nCoefficients) || nCoefficients < ncol(X))
		stop("The fitted mmes object has incompatible coefficient dimensions.", call.=FALSE)
	D <- cbind(Matrix::Diagonal(ncol(X)),
						 Matrix::Matrix(0, nrow=ncol(X), ncol=nCoefficients - ncol(X),
														sparse=TRUE))

	list(X=X, bhat=as.numeric(object$b),
			 V=as.matrix(predict_mmes_vcov_cpp(object, D)),
			 nbasis=estimability::all.estble,
			 dffun=function(k, dfargs) Inf, dfargs=list())
}
# 
# emm_basis.mmer = function(object, trms, xlev, grid, ...) {
#     cf = object$Beta
#     bhat = cf$Estimate
#     m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
#     # if we can get contrasts from the object, fix next line
#     X = model.matrix(trms, m, contrasts.arg = NULL)
#     V = .my.vcov(object, vcov. = function(., ...) .$VarBeta)
#     
#     nbasis = estimability::all.estble   # soup this up if can have rank deficiencies
#     misc = list()
#     # soup-up following if (1) glms allowed or (2) d.f. available
#     dfargs = list(df = object$df.residual)
#     dffun = function(k, dfargs) Inf
#     bas = list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
#                dffun = dffun, dfargs = dfargs, misc = misc)
#     # check for multiv resp
#     k = length(levels(cf$Trait))
#     if (k > 1) {
#         bas$misc$ylevs = list(Trait = levels(cf$Trait))
#         bas$X = kronecker(diag(rep(1, k)), bas$X)
#         # reorder coefs to go one trait at a time
#         ord = as.integer(matrix(seq_along(bas$bhat), ncol = k, byrow = TRUE))
#         bas$bhat = bas$bhat[ord]
#         bas$V = bas$V[ord, ord, drop = FALSE]
#     }
#     bas
# }